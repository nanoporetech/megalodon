import argparse
from collections import defaultdict

import numpy as np

from megalodon import mods, megalodon_helper as mh


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'megalodon_results_dir',
        help='Output directory from megalodon with per_read_mods in output.')
    parser.add_argument(
        '--out-filename',
        help='Output filename for text summary. Default: output into ' +
        'megalodon results directory')

    return parser

def main():
    args = get_parser().parse_args()

    mods_db = mods.ModsDb(
        mh.get_megalodon_fn(args.megalodon_results_dir, mh.PR_MOD_NAME))
    mods_txt_fp = open(
        mh.get_megalodon_fn(args.megalodon_results_dir, mh.PR_MOD_TXT_NAME)
        if args.out_filename is None else args.out_filename, 'w')
    mods_txt_fp.write('\t'.join(mods_db.text_field_names) + '\n')
    for pos_id, pos_chrm, strand, pos in mods_db.iter_pos():
        pr_mod_stats = mods_db.get_pos_stats(
            (pos_id, pos_chrm, strand, pos), return_uuids=True)
        mod_type_stats = defaultdict(dict)
        for r_stats in pr_mod_stats:
            mod_type_stats[r_stats.read_id][r_stats.mod_base] = (
                r_stats.score, r_stats.raw_motif, r_stats.motif_pos,
                r_stats.chrm)

        mod_out_text = ''
        for read_id, r_mod_stats in mod_type_stats.items():
            mod_lps = np.array(list(zip(*r_mod_stats.values()))[0])
            with np.errstate(divide='ignore'):
                can_lp = np.log1p(-np.exp(mod_lps).sum())
            mod_out_text += '\n'.join((
                ('\t'.join('{}' for _ in mods_db.text_field_names)).format(
                    read_id, chrm, strand, pos, mod_lp,
                    can_lp, mod_base, '{}:{}'.format(raw_motif, motif_pos))
                for mod_base, (mod_lp, raw_motif, motif_pos, chrm) in
                r_mod_stats.items())) + '\n'
        mods_txt_fp.write(mod_out_text)

    return

if __name__ == '__main__':
    main()
