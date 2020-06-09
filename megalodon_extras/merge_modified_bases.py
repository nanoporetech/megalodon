from tqdm import tqdm

from megalodon import logging, mods, megalodon_helper as mh
from ._extras_parsers import get_parser_merge_modified_bases


LOGGER = logging.get_logger()


def _main(args):
    mh.mkdir(args.output_megalodon_results_dir, args.overwrite)
    logging.init_logger(args.output_megalodon_results_dir)

    LOGGER.info('Opening new modified base statistics database')
    out_mods_db = mods.ModsDb(
        mh.get_megalodon_fn(args.output_megalodon_results_dir, mh.PR_MOD_NAME),
        read_only=False, init_db_tables=True, in_mem_chrm_to_dbid=True,
        in_mem_mod_to_dbid=True, in_mem_uuid_to_dbid=True,
        in_mem_pos_to_dbid=not args.mod_positions_on_disk)

    LOGGER.info(
        'Merging will proceed in five stages:\n\t1) chrmosomes\n\t2) ' +
        'modified base definitions\n\t3) read identifiers\n\t4)reference ' +
        'positions\n\t5) modified base statistics')
    LOGGER.info('Merging chrm tables')
    ref_names_and_lens = [[], []]
    for mega_dir in args.megalodon_results_dirs:
        mods_db = mods.ModsDb(mh.get_megalodon_fn(mega_dir, mh.PR_MOD_NAME))
        bar = tqdm(desc=mega_dir, total=mods_db.get_num_uniq_chrms(),
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

    LOGGER.info('Merging mod tables')
    for mega_dir in args.megalodon_results_dirs:
        mods_db = mods.ModsDb(mh.get_megalodon_fn(mega_dir, mh.PR_MOD_NAME))
        bar = tqdm(desc=mega_dir, total=mods_db.get_num_uniq_mods(),
                   smoothing=0, dynamic_ncols=True)
        for (_, mod_base, motif, motif_pos,
             raw_motif) in mods_db.iter_mod_bases():
            out_mods_db.get_mod_base_dbid_or_insert(
                mod_base, motif, motif_pos, raw_motif)
            bar.update()
        mods_db.close()
        bar.close()

    LOGGER.info('Merging read uuid tables')
    for mega_dir in args.megalodon_results_dirs:
        mods_db = mods.ModsDb(mh.get_megalodon_fn(mega_dir, mh.PR_MOD_NAME))
        bar = tqdm(desc=mega_dir, total=mods_db.get_num_uniq_reads(),
                   smoothing=0, dynamic_ncols=True)
        for read_dbid, uuid in mods_db.iter_uuids():
            out_mods_db.get_read_dbid_or_insert(uuid)
            bar.update()
        mods_db.close()
        bar.close()

    LOGGER.info('Merging pos tables')
    for mega_dir in args.megalodon_results_dirs:
        mods_db = mods.ModsDb(mh.get_megalodon_fn(mega_dir, mh.PR_MOD_NAME))
        bar = tqdm(desc=mega_dir, total=mods_db.get_num_uniq_mod_pos(),
                   smoothing=0, dynamic_ncols=True)
        for _, chrm_dbid, strand, pos in mods_db.iter_pos():
            out_mods_db.get_pos_dbid_or_insert(
                out_mods_db.get_chrm_dbid(out_mods_db.get_chrm(chrm_dbid)[0]),
                strand, pos)
            bar.update()
        mods_db.close()
        bar.close()

    LOGGER.info('Inserting modified base data')
    for mega_dir in args.megalodon_results_dirs:
        mods_db = mods.ModsDb(mh.get_megalodon_fn(mega_dir, mh.PR_MOD_NAME))
        bar = tqdm(
            desc=mega_dir, total=mods_db.get_num_uniq_stats(), smoothing=0,
            dynamic_ncols=True)
        # TODO could possibly be more efficient by grouping inserts by uuid
        for (score, uuid, mod_base, motif, motif_pos, raw_motif, strand,
             pos, chrm, chrm_len) in mods_db.iter_data():
            chrm_id = out_mods_db.get_chrm_dbid(chrm)
            pos_id = out_mods_db.get_pos_dbid_or_insert(chrm_id, strand, pos)
            mod_base_id = out_mods_db.get_mod_base_dbid_or_insert(
                mod_base, motif, motif_pos, raw_motif)
            read_id = out_mods_db.get_read_dbid_or_insert(uuid)
            out_mods_db.insert_data(score, pos_id, mod_base_id, read_id)
            bar.update()
        bar.close()

    LOGGER.info(
        'Creating data covering index for efficient searching by position')
    out_mods_db.create_data_covering_index()
    out_mods_db.close()


if __name__ == '__main__':
    _main(get_parser_merge_modified_bases().parse_args())
