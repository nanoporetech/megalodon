import sys
import argparse
from collections import defaultdict

import pysam
import numpy as np


HOM_REF_TXT = 'hom_ref'
HET_TXT = 'het'
HOM_ALT_TXT = 'hom_alt'

SNP_TXT = 'SNP'
DEL_TXT = 'DEL'
INS_TXT = 'INS'


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'ground_truth_variants',
        help='VCF file containing ground truth diploid variant calls.')
    parser.add_argument(
        'megalodon_variants', default='megalodon_results/variants.vcf',
        help='VCF file containing diploid variant calls from megalodon.')

    return parser


def main():
    def conv_call_str(gt_vals):
        gt_set = set(gt_vals)
        if gt_set == set([0]):
            return HOM_REF_TXT
        elif gt_set == set([0, 1]):
            return HET_TXT
        return HOM_ALT_TXT

    args = get_parser().parse_args()

    gt_calls = defaultdict(dict)
    for variant in pysam.VariantFile(args.ground_truth_variants).fetch():
        # skip mutli-allelic sites
        if variant.alts is None or len(variant.alts) > 1: continue
        if len(variant.ref) == len(variant.alts[0]):
            gt_calls[SNP_TXT][(variant.contig, variant.pos, variant.ref,
                               variant.alts[0])] = conv_call_str(
                                   variant.samples.values()[0]['GT'])
        elif len(variant.ref) > len(variant.alts[0]):
            gt_calls[DEL_TXT][(variant.contig, variant.pos, variant.ref,
                               variant.alts[0])] = conv_call_str(
                                   variant.samples.values()[0]['GT'])
        else:
            gt_calls[INS_TXT][(variant.contig, variant.pos, variant.ref,
                               variant.alts[0])] = conv_call_str(
                                   variant.samples.values()[0]['GT'])
    mega_calls = defaultdict(dict)
    for variant in pysam.VariantFile(args.megalodon_variants).fetch():
        # skip mutli-allelic sites
        if len(variant.alts) > 1: continue
        # TODO remove minus one from previous off-by-one bug
        if len(variant.ref) == len(variant.alts[0]):
            mega_calls[SNP_TXT][(variant.contig, variant.pos - 1, variant.ref,
                                 variant.alts[0])] = conv_call_str(
                                     variant.samples.values()[0]['GT'])
        elif len(variant.ref) > len(variant.alts[0]):
            mega_calls[DEL_TXT][(variant.contig, variant.pos - 1, variant.ref,
                                 variant.alts[0])] = conv_call_str(
                                     variant.samples.values()[0]['GT'])
        else:
            mega_calls[INS_TXT][(variant.contig, variant.pos - 1, variant.ref,
                                 variant.alts[0])] = conv_call_str(
                                     variant.samples.values()[0]['GT'])

    for var_type in (SNP_TXT, DEL_TXT, INS_TXT):
        counts = defaultdict(int)
        for chrm_pos_ref_alt in set(gt_calls[var_type]).intersection(
                mega_calls[var_type]):
            counts[(gt_calls[var_type][chrm_pos_ref_alt],
                    mega_calls[var_type][chrm_pos_ref_alt])] += 1
        sys.stdout.write(var_type + '\n')
        sys.stdout.write('{:25}{:<10}{:<10}{:<10}\n'.format(
            'Truth\Calls', 'HomRef', 'Het', 'HomAlt'))
        for truth in (HOM_REF_TXT, HET_TXT, HOM_ALT_TXT):
            sys.stdout.write('{:25}{:<10}{:<10}{:<10}\n'.format(
                truth, *map(str, (
                    counts[(truth, mega_call)]
                    for mega_call in (HOM_REF_TXT, HET_TXT, HOM_ALT_TXT)))))
        sys.stdout.write('\n\n')

    return

if __name__ == '__main__':
    main()
