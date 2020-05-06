from megalodon import megalodon, variants, megalodon_helper as mh
from._extras_parsers import get_parser_merge_variants


def _main(args):
    megalodon.mkdir(args.output_megalodon_results_dir, False)
    out_vars_db = variants.VarsDb(
        mh.get_megalodon_fn(args.output_megalodon_results_dir, mh.PR_VAR_NAME),
        read_only=False, loc_index_in_memory=not args.var_locations_on_disk)

    for mega_dir in args.megalodon_results_dirs:
        # full read only mode with no indices read into memory
        vars_db = variants.VarsDb(
            mh.get_megalodon_fn(mega_dir, mh.PR_VAR_NAME),
            read_only=True, chrm_index_in_memory=False,
            alt_index_in_memory=False, uuid_index_in_memory=False)
        for (score, uuid, strand, alt_seq, ref_seq, pos, var_name,
             test_end, test_start, chrm, chrm_len) in vars_db.iter_data():
            chrm_id = out_vars_db.get_chrm_id_or_insert(chrm, chrm_len)
            loc_id = out_vars_db.get_loc_id_or_insert(
                chrm_id, test_start, test_end, pos, ref_seq, var_name)
            alt_id = out_vars_db.get_alt_id_or_insert(alt_seq)
            read_id = out_vars_db.get_read_id_or_insert(uuid)
            out_vars_db.insert_data(score, loc_id, alt_id, read_id)

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
