import pysam

from ._extras_parsers import get_parser_phase_variants_extract_haplotype_reads


def _main(args):
    out_fps = {}
    for rec in pysam.AlignmentFile(args.alignment_filename):
        try:
            hp = dict(rec.tags)['HP']
        except KeyError:
            # skip un-tagged reads
            continue
        if hp not in out_fps:
            out_fps[hp] = open('{}.haplotype_{}_read_ids.txt'.format(
                args.out_basename, hp), 'w')
        out_fps[hp].write(rec.qname + '\n')

    for fp in out_fps.values():
        fp.close()

    return


if __name__ == '__main__':
    _main(get_parser_phase_variants_extract_haplotype_reads().parse_args())
