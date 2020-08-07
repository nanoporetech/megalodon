import queue
from time import sleep
import multiprocessing as mp

from tqdm import tqdm

from megalodon import (
    backends, logging, mods, megalodon_helper as mh,
    megalodon_multiprocessing as mega_mp)
from ._extras_parsers import get_parser_merge_modified_bases


LOGGER = logging.get_logger()

QUEUE_SIZE_LIMIT = 100


########################
# data table functions #
########################

def extract_data_worker(in_db_fns_q, data_conn, out_mods_db_fn, batch_size):
    # load output database with uuid in-memory indices
    out_mods_db = mods.ModsDb(out_mods_db_fn, in_mem_uuid_to_dbid=True)
    while True:
        try:
            in_mod_db_fn = in_db_fns_q.get(block=True, timeout=0.1)
        except queue.Empty:
            sleep(0.001)
            continue
        if in_mod_db_fn is None:
            break

        in_mods_db = mods.ModsDb(in_mod_db_fn)
        batch_data = []
        for score, uuid, mod_base, in_pos_dbid in in_mods_db.iter_data():
            out_pos_dbid = out_mods_db.get_pos_dbid(*in_mods_db.get_pos(
                in_pos_dbid))
            batch_data.append((
                score, out_pos_dbid, out_mods_db.get_mod_base_dbid(mod_base),
                out_mods_db.get_read_dbid(uuid)))
            if len(batch_data) >= batch_size:
                data_conn.put(batch_data)
                batch_data = []
        if len(batch_data) > 0:
            data_conn.put(batch_data)
            batch_data = []
        in_mods_db.close()
    out_mods_db.close()


def insert_data_mp(
        in_mod_db_fns, out_mods_db, out_mods_db_fn, batch_size, max_proc):
    LOGGER.info('Merging modified base data using multiprocessing')
    num_proc = min(max_proc, len(in_mod_db_fns))
    in_db_fns_q = mp.Queue()
    for in_mod_db_fn in in_mod_db_fns:
        in_db_fns_q.put(in_mod_db_fn)
    for _ in range(num_proc):
        in_db_fns_q.put(None)
    data_q = mega_mp.SimplexManyToOneQueue(max_size=QUEUE_SIZE_LIMIT)
    data_ps = []
    for _ in range(num_proc):
        data_conn = data_q.get_conn()
        p = mp.Process(
            target=extract_data_worker,
            args=(in_db_fns_q, data_conn, out_mods_db_fn, batch_size),
            daemon=True)
        p.start()
        data_conn.close()
        del data_conn
        data_ps.append(p)

    total_batches = 0
    for in_mod_db_fn in in_mod_db_fns:
        in_mods_db = mods.ModsDb(in_mod_db_fn)
        total_batches += (in_mods_db.get_num_uniq_stats() // batch_size) + 1
        in_mods_db.close()
    bar = tqdm(desc='Data Batches', unit='Batches', total=total_batches,
               smoothing=0, dynamic_ncols=True)

    while data_q.has_valid_conns:
        for batch_data in data_q.wait_recv():
            out_mods_db.insert_batch_data(batch_data)
            bar.update()
    bar.close()


def insert_data(in_mod_db_fns, out_mods_db, batch_size):
    LOGGER.info('Inserting modified base data')
    total_batches = 0
    for in_mod_db_fn in in_mod_db_fns:
        in_mods_db = mods.ModsDb(in_mod_db_fn)
        total_batches += (in_mods_db.get_num_uniq_stats() // batch_size) + 1
        in_mods_db.close()
    bar = tqdm(desc='Data Batches', unit='Batches', total=total_batches,
               smoothing=0, dynamic_ncols=True)
    for in_mod_db_fn in in_mod_db_fns:
        in_mods_db = mods.ModsDb(in_mod_db_fn)
        batch_data = []
        for score, uuid, mod_base, in_pos_dbid in in_mods_db.iter_data():
            out_pos_dbid = out_mods_db.get_pos_dbid(*in_mods_db.get_pos(
                in_pos_dbid))
            batch_data.append((
                score, out_pos_dbid, out_mods_db.get_mod_base_dbid(mod_base),
                out_mods_db.get_read_dbid(uuid)))
            if len(batch_data) >= batch_size:
                out_mods_db.insert_batch_data(batch_data)
                batch_data = []
                bar.update()
        if len(batch_data) > 0:
            out_mods_db.insert_batch_data(batch_data)
            bar.update()
        in_mods_db.close()
    bar.close()


########################
# read table functions #
########################

def extract_reads_worker(in_mod_db_fn, uuids_q):
    in_mods_db = mods.ModsDb(in_mod_db_fn)
    in_uuids = set(uuid for _, uuid in in_mods_db.iter_uuids())
    in_mods_db.close()
    uuids_q.put(in_uuids)


def insert_reads_mp(in_mod_db_fns, out_mods_db):
    LOGGER.info('Extracting read uuid tables using multiprocessing')
    uuids_q = mp.Queue()
    uuids_ps = []
    for in_mod_db_fn in in_mod_db_fns:
        p = mp.Process(
            target=extract_reads_worker,
            args=(in_mod_db_fn, uuids_q), daemon=True)
        p.start()
        uuids_ps.append(p)

    bar = tqdm(desc='Databases', unit='DBs', total=len(in_mod_db_fns),
               smoothing=0, dynamic_ncols=True)
    in_uuids = set()
    while any(p.is_alive() for p in uuids_ps):
        try:
            in_db_uuids = uuids_q.get(block=True, timeout=0.1)
        except queue.Empty:
            continue
        in_uuids.update(in_db_uuids)
        bar.update()
    while not uuids_q.empty():
        in_db_uuids = uuids_q.get(block=False)
        in_uuids.update(in_db_uuids)
        bar.update()
    bar.close()
    LOGGER.info('Inserting read uuid tables into output DB')
    out_mods_db.insert_uuids(in_uuids)


def insert_reads(in_mod_db_fns, out_mods_db):
    LOGGER.info('Merging read uuid tables')
    in_uuids = set()
    for in_mod_db_fn in tqdm(in_mod_db_fns, desc='Databases', unit='DBs',
                             smoothing=0, dynamic_ncols=True):
        in_mods_db = mods.ModsDb(in_mod_db_fn)
        in_uuids.update(uuid for _, uuid in in_mods_db.iter_uuids())
        in_mods_db.close()
    out_mods_db.insert_uuids(in_uuids)


############################
# mod/chrm table functions #
############################

def extract_mods(in_mod_db_fns):
    LOGGER.info('Merging mod tables')
    # collect modified base definitions from input databases
    mod_base_to_can = dict()
    for in_mod_db_fn in tqdm(in_mod_db_fns, desc='Databases', unit='DBs',
                             smoothing=0, dynamic_ncols=True):
        mods_db = mods.ModsDb(in_mod_db_fn)
        for mod_base, can_base, mln in mods_db.get_full_mod_data():
            if mod_base in mod_base_to_can and \
               (can_base, mln) != mod_base_to_can[mod_base]:
                raise mh.MegaError(
                    'Modified base associated with mutliple canonical bases ' +
                    'or long names in different databases. {} != {}'.format(
                        str((can_base, mln)),
                        str(mod_base_to_can[mod_base])))
            mod_base_to_can[mod_base] = (can_base, mln)
    # check that mod long names are unique
    mlns = [mln for _, mln in mod_base_to_can.values()]
    if len(mlns) != len(set(mlns)):
        raise mh.MegaError(
            'Modified base long name assigned to more than one modified ' +
            'base single letter code.')

    # extract canonical bases associated with modified base
    can_bases = set(can_base for can_base, _ in mod_base_to_can.values())
    # determine first valid canonical alphabet compatible with database
    can_alphabet = None
    for v_alphabet in mh.VALID_ALPHABETS:
        if len(can_bases.difference(v_alphabet)) == 0:
            can_alphabet = v_alphabet
            break
    if can_alphabet is None:
        LOGGER.error(
            'Mods database does not contain valid canonical bases ({})'.format(
                ''.join(sorted(can_bases))))
        raise mh.MegaError('Invalid alphabet.')

    # compute full output alphabet and ordered modified base long names
    can_base_to_mods = dict(
        (can_base, [(mod_base, mln)
                    for mod_base, (mcan_base, mln) in mod_base_to_can.items()
                    if mcan_base == can_base])
        for can_base in can_alphabet)
    alphabet = ''
    mod_long_names = []
    for can_base in can_alphabet:
        alphabet += can_base
        for mod_base, mln in can_base_to_mods[can_base]:
            alphabet += mod_base
            mod_long_names.append(mln)

    return alphabet, mod_long_names


def extract_chrms(in_mod_db_fns):
    LOGGER.info('Merging chrm tables')
    ref_names_and_lens = [[], []]
    for in_mod_db_fn in tqdm(in_mod_db_fns, desc='Databases', unit='DBs',
                             smoothing=0, dynamic_ncols=True):
        mods_db = mods.ModsDb(in_mod_db_fn)
        for _, chrm, chrm_len in mods_db.iter_chrms():
            if chrm in ref_names_and_lens[0]:
                prev_chrm_len = ref_names_and_lens[1][
                    ref_names_and_lens[0].index(chrm)]
                if prev_chrm_len != chrm_len:
                    raise mh.MegaError((
                        'Chromosome lengths from databases do not agree ' +
                        'for {}: {} != {}').format(
                            chrm, prev_chrm_len, chrm_len))
            else:
                ref_names_and_lens[0].append(chrm)
                ref_names_and_lens[1].append(chrm_len)
        mods_db.close()
    return ref_names_and_lens


########
# main #
########

def _main(args):
    mh.mkdir(args.output_megalodon_results_dir, args.overwrite)
    logging.init_logger(args.output_megalodon_results_dir)

    LOGGER.info('Extracting mods and chrms from input databases')
    in_mod_db_fns = [mh.get_megalodon_fn(mega_dir, mh.PR_MOD_NAME)
                     for mega_dir in args.megalodon_results_dirs]
    alphabet, mod_long_names = extract_mods(in_mod_db_fns)
    ref_names_and_lens = extract_chrms(in_mod_db_fns)

    LOGGER.info('Opening new per-read modified base statistics database')
    model_info = backends.DetachedModelInfo(
        alphabet=alphabet, mod_long_names=mod_long_names)
    mods_info = mods.ModInfo(
        model_info, out_dir=args.output_megalodon_results_dir)
    mods.init_mods_db(mods_info, ref_names_and_lens)

    # load uuids in memory in main out db only in single process mode.
    # else worker threads only have to load uuid lookup tables
    out_mods_db = mods.ModsDb(
        mods_info.mods_db_fn, read_only=False,
        in_mem_uuid_to_dbid=args.single_process)

    LOGGER.info('Inserting read UUIDs from input databases')
    if args.single_process:
        insert_reads(in_mod_db_fns, out_mods_db)
    else:
        insert_reads_mp(in_mod_db_fns, out_mods_db)
    # commit so read uuids are available to worker processes
    out_mods_db.commit()
    LOGGER.info('Inserting per-read calls from input databases')
    if args.single_process:
        insert_data(in_mod_db_fns, out_mods_db, args.data_batch_size)
    else:
        insert_data_mp(
            in_mod_db_fns, out_mods_db, mods_info.mods_db_fn,
            args.data_batch_size,
            args.max_processes)
    out_mods_db.commit()

    LOGGER.info(
        'Creating data covering index for efficient iteration by position')
    out_mods_db.create_data_covering_index()
    out_mods_db.commit()
    out_mods_db.close()


if __name__ == '__main__':
    _main(get_parser_merge_modified_bases().parse_args())
