from tqdm import tqdm

from megalodon import logging, mods, megalodon_helper as mh
from ._extras_parsers import get_parser_merge_modified_bases


def _main(args):
    mh.mkdir(args.output_megalodon_results_dir, args.overwrite)
    logging.init_logger(args.output_megalodon_results_dir)
    logger = logging.get_logger()

    logger.info('Opening new modified base statistics database')
    out_mods_db = mods.ModsDb(
        mh.get_megalodon_fn(args.output_megalodon_results_dir, mh.PR_MOD_NAME),
        read_only=False, in_mem_chrm_to_dbid=True, in_mem_mod_to_dbid=True,
        in_mem_uuid_to_dbid=True,
        in_mem_pos_to_dbid=not args.mod_positions_on_disk)

    # TODO fill chrm table, then mods, then uuid, then pos
    for mega_dir in args.megalodon_results_dirs:
        logger.info('Adding modified base statistics from {}'.format(
            mega_dir))
        # full read only mode with no indices read into memory
        mods_db = mods.ModsDb(
            mh.get_megalodon_fn(mega_dir, mh.PR_MOD_NAME),
            read_only=True, chrm_index_in_memory=False,
            mod_index_in_memory=False, uuid_index_in_memory=False)
        bar = tqdm(
            desc=mega_dir, total=mods_db.get_num_uniq_stats(), smoothing=0,
            dynamic_ncols=True)
        for (score, uuid, mod_base, motif, motif_pos, raw_motif, strand,
             pos, chrm, chrm_len) in mods_db.iter_data():
            chrm_id = out_mods_db.get_chrm_dbid_or_insert(chrm, chrm_len)
            pos_id = out_mods_db.get_pos_dbid_or_insert(chrm_id, strand, pos)
            mod_base_id = out_mods_db.get_mod_base_dbid_or_insert(
                mod_base, motif, motif_pos, raw_motif)
            read_id = out_mods_db.get_read_dbid_or_insert(uuid)
            out_mods_db.insert_data(score, pos_id, mod_base_id, read_id)
            bar.update()
        bar.close()

    logger.info('Creating indices and closing database')
    out_mods_db.create_chrm_index()
    out_mods_db.create_data_covering_index()
    out_mods_db.close()


if __name__ == '__main__':
    _main(get_parser_merge_modified_bases().parse_args())
