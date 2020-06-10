from megalodon import megalodon_helper as mh
from ._extras_parsers import get_parser_modified_bases_create_ground_truth


def _main(args):
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
                        chrom, mh.int_strand_to_str(strand), pos,
                        'False'))) + '\n')
                    if args.strand_offset is not None:
                        gt_fp.write(','.join(map(str, (
                            chrom, mh.int_strand_to_str(strand),
                            pos + args.strand_offset, 'False'))) + '\n')
                elif pct_mod >= args.pct_mod_thresholds[1]:
                    gt_fp.write(','.join(map(str, (
                        chrom, mh.int_strand_to_str(strand), pos,
                        'True'))) + '\n')
                    if args.strand_offset is not None:
                        gt_fp.write(','.join(map(str, (
                            chrom, mh.int_strand_to_str(strand),
                            pos + args.strand_offset, 'True'))) + '\n')


if __name__ == '__main__':
    _main(get_parser_modified_bases_create_ground_truth().parse_args())
