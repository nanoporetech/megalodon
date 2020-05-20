from tqdm import tqdm

from megalodon import logging, megalodon_helper as mh, variants
from._extras_parsers import get_parser_merge_variants


def _main(args):
    mh.mkdir(args.output_megalodon_results_dir, args.overwrite)
    logging.init_logger(args.output_megalodon_results_dir)
    logger = logging.get_logger()

    logger.info('Opening new sequence variant statistics database')
    out_vars_db = variants.VarsDb(
        mh.get_megalodon_fn(args.output_megalodon_results_dir, mh.PR_VAR_NAME),
        read_only=False, loc_index_in_memory=not args.var_locations_on_disk,
        uuid_index_in_memory=True)

    for mega_dir in args.megalodon_results_dirs:
        logger.info('Adding sequence variant statistics from {}'.format(
            mega_dir))
        # full read only mode with no indices read into memory
        vars_db = variants.VarsDb(
            mh.get_megalodon_fn(mega_dir, mh.PR_VAR_NAME),
            read_only=True, chrm_index_in_memory=False,
            alt_index_in_memory=False, uuid_index_in_memory=False)
        bar = tqdm(
            desc=mega_dir, total=vars_db.get_num_uniq_stats(), smoothing=0,
            dynamic_ncols=True)
        for (score, uuid, strand, alt_seq, ref_seq, pos, var_name,
             test_end, test_start, chrm, chrm_len) in vars_db.iter_data():
            chrm_id = out_vars_db.get_chrm_id_or_insert(chrm, chrm_len)
            loc_id = out_vars_db.get_loc_id_or_insert(
                chrm_id, test_start, test_end, pos, ref_seq, var_name)
            alt_id = out_vars_db.get_alt_id_or_insert(alt_seq)
            read_id = out_vars_db.get_read_id_or_insert(uuid)
            out_vars_db.insert_data(score, loc_id, alt_id, read_id)
            bar.update()
        bar.close()

    logger.info('Creating indices and closing database')
    if out_vars_db.chrm_idx_in_mem:
        out_vars_db.create_chrm_index()
    if out_vars_db.loc_idx_in_mem:
        out_vars_db.create_loc_index()
    if out_vars_db.alt_idx_in_mem:
        out_vars_db.create_alt_index()
    out_vars_db.create_data_covering_index()
    out_vars_db.close()


if __name__ == '__main__':
    _main(get_parser_merge_variants().parse_args())
