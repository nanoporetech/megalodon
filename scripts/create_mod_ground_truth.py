import argparse

from megalodon import megalodon_helper as mh


STRAND_CONV = {1: '+', -1: '-', None: '.'}


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--bed-methyl-files', nargs='+', required=True,
        help='Bed methyl files.')
    parser.add_argument(
        '--coverage-threshold', type=int, default=1,
        help='Only include sites with sufficient coverage. ' +
        'Default: 1 (= All sites)')
    parser.add_argument(
        '--pct-mod-thresholds', type=float, nargs=2, default=[10.0, 90.0],
        help='Lower and upper percent modified thresholds for ground truth ' +
        'modified positions. Default: %(default)s')
    parser.add_argument(
        '--out-csv', default='ground_truth_modifications.csv',
        help='Output filename for ground truth calls. Default: %(default)s')
    parser.add_argument(
        '--strand-offset', type=int,
        help='Offset to combine stranded results. Positive value indicates ' +
        'reverse strand sites have higher position values. Default treat ' +
        'strands independently.')

    return parser


def main():
    args = get_parser().parse_args()
    samp_cov, samp_mod_cov = mh.parse_bed_methyls(
        args.bed_methyl_files, strand_offset=args.strand_offset)
    with open(args.out_csv, 'w') as gt_fp:
        for (chrom, strand), ctg_cov in samp_cov.items():
            for pos, cov in ctg_cov.items():
                if cov < args.coverage_threshold:
                    continue
                pct_mod = 100 * samp_mod_cov[(chrom, strand)][pos] / cov
                if pct_mod <= args.pct_mod_thresholds[0]:
                    gt_fp.write(','.join(map(str, (
                        chrom, STRAND_CONV[strand], pos, 'False'))) + '\n')
                elif pct_mod >= args.pct_mod_thresholds[1]:
                    gt_fp.write(','.join(map(str, (
                        chrom, STRAND_CONV[strand], pos, 'True'))) + '\n')


if __name__ == '__main__':
    main()
