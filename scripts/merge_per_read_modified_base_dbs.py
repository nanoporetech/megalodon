import argparse

from megalodon import megalodon, mods, megalodon_helper as mh


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'megalodon_results_dirs', nargs='+',
        help='Output megalodon directories with per_read_mods in output.')
    parser.add_argument(
        '--output-megalodon-results-dir',
        default='megalodon_merge_mods_results',
        help='Output directory. Cannot exist before this command. ' +
        'Default: %(default)s')
    parser.add_argument(
        '--mod-positions-on-disk', action='store_true',
        help='Force modified base positions to be stored only within on ' +
        'disk database table. This option will reduce the RAM memory ' +
        'requirement, but may slow processing. Default: ' +
        'Store positions in memory.')

    return parser

def main():
    args = get_parser().parse_args()

    megalodon.mkdir(args.output_megalodon_results_dir, False)
    out_mods_db = mods.ModsDb(
        mh.get_megalodon_fn(args.output_megalodon_results_dir, mh.PR_MOD_NAME),
        read_only=False, pos_index_in_memory=not args.mod_positions_on_disk)

    for mega_dir in args.megalodon_results_dirs:
        # full read only mode with no indices read into memory
        mods_db = mods.ModsDb(
            mh.get_megalodon_fn(mega_dir, mh.PR_MOD_NAME),
            read_only=True, chrm_index_in_memory=False,
            mod_index_in_memory=False, uuid_index_in_memory=False)
        for (score, uuid, mod_base, motif, motif_pos, raw_motif, strand,
             pos, chrm, chrm_len) in mods_db.iter_data():
            chrm_id = out_mods_db.get_chrm_id_or_insert(chrm, chrm_len)
            pos_id = out_mods_db.get_pos_id_or_insert(chrm_id, strand, pos)
            mod_base_id = out_mods_db.get_mod_base_id_or_insert(
                mod_base, motif, motif_pos, raw_motif)
            read_id =  out_mods_db.get_read_id_or_insert(uuid)
            out_mods_db.insert_data(score, pos_id, mod_base_id, read_id)

    if out_mods_db.chrm_idx_in_mem:
        out_mods_db.create_chrm_index()
    if out_mods_db.pos_idx_in_mem:
        out_mods_db.create_pos_index()
    if out_mods_db.mod_idx_in_mem:
        out_mods_db.create_mod_index()
    out_mods_db.create_data_covering_index()
    out_mods_db.close()

    return

if __name__ == '__main__':
    main()
