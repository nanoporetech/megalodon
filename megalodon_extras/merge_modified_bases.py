import queue
from time import sleep
import multiprocessing as mp

from tqdm import tqdm

from megalodon import logging, mods, megalodon_helper as mh
from ._extras_parsers import get_parser_merge_modified_bases


LOGGER = logging.get_logger()

QUEUE_SIZE_LIMIT = 100


########################
# data table functions #
########################

def get_data_dbids(out_mods_db, chrm, strand, pos, mod_data, uuid):
    # extract output database ids
    pos_dbid = out_mods_db.get_pos_dbid_or_insert(
        out_mods_db.get_chrm_dbid(chrm), strand, pos)
    mod_base_dbid = out_mods_db.get_mod_base_dbid_or_insert(*mod_data)
    read_dbid = out_mods_db.get_read_dbid_or_insert(uuid)
    return pos_dbid, mod_base_dbid, read_dbid


def extract_data_worker(
        in_db_fns_q, data_q, out_mods_db_fn, batch_size, force_uint32,
        db_safety):
    # load output database with all in memory indices
    out_mods_db = mods.ModsDb(
        out_mods_db_fn, read_only=True,
        in_mem_chrm_to_dbid=True, in_mem_mod_to_dbid=True,
        in_mem_uuid_to_dbid=True, in_mem_pos_to_dbid=True,
        force_uint32_pos_to_dbid=force_uint32, db_safety=db_safety)
    while True:
        try:
            in_mod_db_fn = in_db_fns_q.get(block=True, timeout=1)
        except queue.Empty:
            sleep(0.001)
            continue
        if in_mod_db_fn is None:
            break

        mods_db = mods.ModsDb(in_mod_db_fn)
        batch_data = []
        for (score, uuid, mod_base, motif, motif_pos, raw_motif, strand, pos,
             chrm, chrm_len) in mods_db.iter_data():
            batch_data.append((score, *get_data_dbids(
                out_mods_db, chrm, strand, pos,
                (mod_base, motif, motif_pos, raw_motif), uuid)))
            if len(batch_data) >= batch_size:
                data_q.put(batch_data)
                batch_data = []
        if len(batch_data) > 0:
            data_q.put(batch_data)
        mods_db.close()
        out_mods_db.db.commit()
    out_mods_db.close()


def insert_data_mp(
        in_mod_db_fns, out_mods_db, out_mods_db_fn, batch_size, max_proc,
        force_uint32, db_safety):
    LOGGER.info('Merging modified base data using multiprocessing')
    num_proc = min(max_proc, len(in_mod_db_fns))
    in_db_fns_q = mp.Queue()
    for in_mod_db_fn in in_mod_db_fns:
        in_db_fns_q.put(in_mod_db_fn)
    for _ in range(num_proc):
        in_db_fns_q.put(None)
    data_q = mp.Queue(maxsize=QUEUE_SIZE_LIMIT)
    data_ps = []
    for _ in range(num_proc):
        p = mp.Process(
            target=extract_data_worker,
            args=(in_db_fns_q, data_q, out_mods_db_fn, batch_size,
                  force_uint32, db_safety), daemon=True)
        p.start()
        data_ps.append(p)

    total_batches = 0
    for in_mod_db_fn in in_mod_db_fns:
        mods_db = mods.ModsDb(in_mod_db_fn)
        total_batches += (mods_db.get_num_uniq_stats() // batch_size) + 1
        mods_db.close()
    bar = tqdm(desc='Statistics Batches', total=total_batches,
               smoothing=0, dynamic_ncols=True)
    while any(p.is_alive() for p in data_ps):
        try:
            batch_data = data_q.get(block=True, timeout=1)
        except queue.Empty:
            sleep(0.001)
            continue
        out_mods_db.insert_read_data(batch_data)
        bar.update()
    while not data_q.empty():
        batch_data = data_q.get(block=False)
        out_mods_db.insert_read_data(batch_data)
        bar.update()
    bar.close()


def insert_data(in_mod_db_fns, out_mods_db, batch_size):
    LOGGER.info('Inserting modified base data')
    for in_mod_db_fn in in_mod_db_fns:
        mods_db = mods.ModsDb(in_mod_db_fn)
        bar = tqdm(
            desc=in_mod_db_fn, total=mods_db.get_num_uniq_stats(), smoothing=0,
            dynamic_ncols=True)
        batch_data = []
        for (score, uuid, mod_base, motif, motif_pos, raw_motif, strand,
             pos, chrm, chrm_len) in mods_db.iter_data():
            batch_data.append((score, *get_data_dbids(
                out_mods_db, chrm, strand, pos,
                (mod_base, motif, motif_pos, raw_motif), uuid)))
            if len(batch_data) >= batch_size:
                out_mods_db.insert_read_data(batch_data)
                batch_data = []
            bar.update()
        if len(batch_data) > 0:
            out_mods_db.insert_read_data(batch_data)
        mods_db.close()
        bar.close()


#######################
# pos table functions #
#######################

def init_pos_dict(mods_db):
    return dict(
        ((chrm, strand), [])
        for _, chrm, _ in mods_db.iter_chrms()
        for strand in (1, -1))


def insert_pos_data(dir_pos, out_mods_db):
    for (chrm, strand), cs_pos in dir_pos.items():
        out_mods_db.get_pos_dbids_or_insert(
            cs_pos, out_mods_db.get_chrm_dbid(chrm), strand)


def extract_pos_worker(in_mod_db_fn, batch_size, pos_q):
    mods_db = mods.ModsDb(in_mod_db_fn)
    pos_batch = init_pos_dict(mods_db)
    num_pos = 0
    for _, chrm_dbid, strand, pos in mods_db.iter_pos():
        pos_batch[(mods_db.get_chrm(chrm_dbid)[0], strand)].append(pos)
        num_pos += 1
        if num_pos >= batch_size:
            pos_q.put(pos_batch)
            pos_batch = init_pos_dict(mods_db)
            num_pos = 0
    if num_pos >= 0:
        pos_q.put(pos_batch)
    mods_db.close()


def insert_pos_mp(in_mod_db_fns, out_mods_db, batch_size):
    LOGGER.info('Merging pos tables using multiprocessing')
    total_batches = 0
    pos_q = mp.Queue(maxsize=QUEUE_SIZE_LIMIT)
    pos_ps = []
    for in_mod_db_fn in in_mod_db_fns:
        mods_db = mods.ModsDb(in_mod_db_fn)
        total_batches += (mods_db.get_num_uniq_mod_pos() // batch_size) + 1
        mods_db.close()
        p = mp.Process(
            target=extract_pos_worker,
            args=(in_mod_db_fn, batch_size, pos_q), daemon=True)
        p.start()
        pos_ps.append(p)

    bar = tqdm(desc='Position Batches', total=total_batches,
               smoothing=0, dynamic_ncols=True)
    while any(p.is_alive() for p in pos_ps):
        try:
            pos_batch = pos_q.get(block=True, timeout=1)
        except queue.Empty:
            sleep(0.001)
            continue
        insert_pos_data(pos_batch, out_mods_db)
        bar.update()
    while not pos_q.empty():
        pos_batch = pos_q.get(block=False)
        insert_pos_data(pos_batch, out_mods_db)
        bar.update()
    bar.close()


def insert_pos(in_mod_db_fns, out_mods_db, batch_size):
    LOGGER.info('Merging pos tables')
    for in_mod_db_fn in in_mod_db_fns:
        mods_db = mods.ModsDb(in_mod_db_fn)
        num_pos = 0
        pos_batch = init_pos_dict(mods_db)
        bar = tqdm(desc=in_mod_db_fn, total=mods_db.get_num_uniq_mod_pos(),
                   smoothing=0, dynamic_ncols=True)
        for _, chrm_dbid, strand, pos in mods_db.iter_pos():
            pos_batch[(mods_db.get_chrm(chrm_dbid)[0], strand)].append(pos)
            num_pos += 1
            if num_pos >= batch_size:
                insert_pos_data(pos_batch, out_mods_db)
                num_pos = 0
                pos_batch = init_pos_dict(mods_db)
            bar.update()
        if num_pos > 0:
            insert_pos_data(pos_batch, out_mods_db)
        mods_db.close()
        bar.close()


########################
# read table functions #
########################

def extract_reads_worker(in_mod_db_fn, batch_size, uuids_q):
    mods_db = mods.ModsDb(in_mod_db_fn)
    uuids_batch = []
    for read_dbid, uuid in mods_db.iter_uuids():
        uuids_batch.append(uuid)
        if len(uuids_batch) >= batch_size:
            uuids_q.put(uuids_batch)
            uuids_batch = []
    if len(uuids_batch) >= 0:
        uuids_q.put(uuids_batch)
    mods_db.close()


def insert_reads_mp(in_mod_db_fns, out_mods_db, batch_size):
    LOGGER.info('Merging read uuid tables using multiprocessing')
    total_batches = 0
    uuids_q = mp.Queue()
    uuids_ps = []
    for in_mod_db_fn in in_mod_db_fns:
        mods_db = mods.ModsDb(in_mod_db_fn)
        total_batches += (mods_db.get_num_uniq_reads() // batch_size) + 1
        mods_db.close()
        p = mp.Process(
            target=extract_reads_worker,
            args=(in_mod_db_fn, batch_size, uuids_q), daemon=True)
        p.start()
        uuids_ps.append(p)

    bar = tqdm(desc='UUID Batches', total=total_batches,
               smoothing=0, dynamic_ncols=True)
    while any(p.is_alive() for p in uuids_ps):
        try:
            uuids_batch = uuids_q.get(block=True, timeout=1)
        except queue.Empty:
            sleep(0.001)
            continue
        out_mods_db.get_read_dbids_or_insert(uuids_batch)
        bar.update()
    while not uuids_q.empty():
        uuids_batch = uuids_q.get(block=False)
        out_mods_db.get_read_dbids_or_insert(uuids_batch)
        bar.update()
    bar.close()


def insert_reads(in_mod_db_fns, out_mods_db):
    LOGGER.info('Merging read uuid tables')
    for in_mod_db_fn in in_mod_db_fns:
        mods_db = mods.ModsDb(in_mod_db_fn)
        bar = tqdm(desc=in_mod_db_fn, total=mods_db.get_num_uniq_reads(),
                   smoothing=0, dynamic_ncols=True)
        for read_dbid, uuid in mods_db.iter_uuids():
            out_mods_db.get_read_dbid_or_insert(uuid)
            bar.update()
        mods_db.close()
        bar.close()


############################
# mod/chrm table functions #
############################

def insert_mods(in_mod_db_fns, out_mods_db):
    LOGGER.info('Merging mod tables')
    all_mod_long_names = set()
    for in_mod_db_fn in in_mod_db_fns:
        mods_db = mods.ModsDb(in_mod_db_fn)
        all_mod_long_names.update(mods_db.get_mod_long_names())
        bar = tqdm(desc=in_mod_db_fn, total=mods_db.get_num_uniq_mods(),
                   smoothing=0, dynamic_ncols=True)
        for (_, mod_base, motif, motif_pos,
             raw_motif) in mods_db.iter_mod_bases():
            out_mods_db.get_mod_base_dbid_or_insert(
                mod_base, motif, motif_pos, raw_motif)
            bar.update()
        mods_db.close()
        bar.close()
    out_mods_db.insert_mod_long_names(list(all_mod_long_names))


def insert_chrms(in_mod_db_fns, out_mods_db):
    LOGGER.info('Merging chrm tables')
    ref_names_and_lens = [[], []]
    for in_mod_db_fn in in_mod_db_fns:
        mods_db = mods.ModsDb(in_mod_db_fn)
        bar = tqdm(desc=in_mod_db_fn, total=mods_db.get_num_uniq_chrms(),
                   smoothing=0, dynamic_ncols=True)
        for _, chrm, chrm_len in mods_db.iter_chrms():
            if chrm not in ref_names_and_lens[0]:
                ref_names_and_lens[0].append(chrm)
                ref_names_and_lens[1].append(chrm_len)
            bar.update()
        mods_db.close()
        bar.close()
    # insert chrms at the end to avoid errors for in memory position datasets
    out_mods_db.insert_chrms(ref_names_and_lens)
    out_mods_db.create_chrm_index()


########
# main #
########

def _main(args):
    mh.mkdir(args.output_megalodon_results_dir, args.overwrite)
    logging.init_logger(args.output_megalodon_results_dir)

    LOGGER.info('Opening new modified base statistics database')
    out_mods_db_fn = mh.get_megalodon_fn(args.output_megalodon_results_dir,
                                         mh.PR_MOD_NAME)
    out_mods_db = mods.ModsDb(
        out_mods_db_fn, read_only=False, init_db_tables=True,
        in_mem_chrm_to_dbid=True, in_mem_mod_to_dbid=True,
        in_mem_uuid_to_dbid=True, in_mem_pos_to_dbid=True,
        force_uint32_pos_to_dbid=args.force_uint32_pos_index,
        db_safety=args.database_safety)

    in_mod_db_fns = [mh.get_megalodon_fn(mega_dir, mh.PR_MOD_NAME)
                     for mega_dir in args.megalodon_results_dirs]

    LOGGER.info(
        'Merging will proceed in five stages:\n\t1) chrmosomes\n\t2) ' +
        'modified base definitions\n\t3) read identifiers\n\t4) reference ' +
        'positions\n\t5) modified base statistics')
    insert_chrms(in_mod_db_fns, out_mods_db)
    insert_mods(in_mod_db_fns, out_mods_db)
    if args.single_process:
        insert_reads(in_mod_db_fns, out_mods_db)
    else:
        insert_reads_mp(in_mod_db_fns, out_mods_db, args.data_batch_size)
    if args.single_process:
        insert_pos(in_mod_db_fns, out_mods_db, args.data_batch_size)
    else:
        insert_pos_mp(in_mod_db_fns, out_mods_db, args.data_batch_size)
    out_mods_db.db.commit()
    if args.single_process:
        insert_data(in_mod_db_fns, out_mods_db, args.data_batch_size)
    else:
        insert_data_mp(
            in_mod_db_fns, out_mods_db, out_mods_db_fn, args.data_batch_size,
            args.max_processes, args.force_uint32_pos_index,
            db_safety=args.database_safety)
    out_mods_db.db.commit()

    LOGGER.info(
        'Creating data covering index for efficient searching by position')
    out_mods_db.create_data_covering_index()
    out_mods_db.db.commit()
    out_mods_db.close()


if __name__ == '__main__':
    _main(get_parser_merge_modified_bases().parse_args())
