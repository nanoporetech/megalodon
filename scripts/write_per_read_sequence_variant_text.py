import argparse
from collections import defaultdict

import numpy as np

from megalodon import variants, megalodon_helper as mh


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'megalodon_results_dir',
        help='Output directory from megalodon with per_read_variants ' +
        'in output.')
    parser.add_argument(
        '--out-filename',
        help='Output filename for text summary. Default: output into ' +
        'megalodon results directory')

    return parser

def main():
    args = get_parser().parse_args()

    vars_db = variants.VarsDb(
        mh.get_megalodon_fn(args.megalodon_results_dir, mh.PR_VAR_NAME),
        uuid_strand_index_in_memory=True)
    vars_txt_fp = open(
        mh.get_megalodon_fn(args.megalodon_results_dir, mh.PR_VAR_TXT_NAME)
        if args.out_filename is None else args.out_filename, 'w')
    vars_txt_fp.write('\t'.join(vars_db.text_field_names) + '\n')
    for loc_id, loc_chrm, pos, ref_seq, var_name in vars_db.iter_locs():
        pr_var_stats = vars_db.get_loc_stats(
            (loc_id, loc_chrm, pos, ref_seq, var_name))
        alt_type_stats = defaultdict(dict)
        for r_stats in pr_var_stats:
            alt_type_stats[r_stats.read_id][r_stats.alt_seq] = (
                r_stats.score, r_stats.chrm)

        var_out_text = ''
        for read_id, r_var_stats in alt_type_stats.items():
            uuid, strand = vars_db.get_uuid_strand(read_id)
            alt_lps = np.array(list(zip(*r_var_stats.values()))[0])
            with np.errstate(divide='ignore'):
                ref_lp = np.log1p(-np.exp(alt_lps).sum())
            var_out_text += '\n'.join((
                ('\t'.join('{}' for _ in vars_db.text_field_names)).format(
                    uuid, chrm, strand, pos, ref_lp, alt_lp,
                    ref_seq, alt_seq, var_name)
                for alt_seq, (alt_lp, chrm) in r_var_stats.items())) + '\n'
        vars_txt_fp.write(var_out_text)

    return

if __name__ == '__main__':
    main()
