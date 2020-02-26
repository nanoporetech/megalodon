import argparse

import pysam


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'alignment_filename', help='Alignment filename.')
    parser.add_argument(
        'out_basename', help='Basename for read ids output.')

    return parser


def main():
    args = get_parser().parse_args()

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
    main()
