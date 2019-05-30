import sys
from collections import defaultdict

import pysam
import numpy as np


GT_SNPS_FN = 'ground_truth_variants.vcf'
MEGA_SNPS_FN = 'megalodon_results/variants.vcf'

HOM_REF_TXT = 'hom_ref'
HET_TXT = 'het'
HOM_ALT_TXT = 'hom_alt'

SNP_TXT = 'SNP'
DEL_TXT = 'DEL'
INS_TXT = 'INS'


def conv_call_str(gt_vals):
    gt_set = set(gt_vals)
    if gt_set == set([0]):
        return HOM_REF_TXT
    elif gt_set == set([0, 1]):
        return HET_TXT
    return HOM_ALT_TXT

gt_calls = defaultdict(dict)
for variant in pysam.VariantsFile(GT_SNPS_FN).fetch():
    # skip mutli-allelic sites
    if len(variants.alts) > 1: continue
    gt_calls[(variant.contig, variant.pos, variant.ref,
              variant.alts[0])] = conv_call_str(
                  variant.samples.values()[0]['GT'])
mega_calls = defaultdict(dict)
for variant in pysam.VariantsFile(MEGA_SNPS_FN).fetch():
    # skip mutli-allelic sites
    if len(variants.alts) > 1: continue
    mega_calls[(variant.contig, variant.pos, variant.ref,
                variant.alts[0])] = conv_call_str(
                    variant.samples.values()[0]['GT'])

calls = defaultdict(int)
for chrm_pos_ref_alt in set(gt_calls).intersection(mega_calls):
    calls[(gt_calls[chrm_pos_ref_alt], mega_calls[chrm_pos_ref_alt])] += 1

sys.stdout.write('{:20}{:<10}{:<10}{:<10}\n'.format(
    'Truth\Calls', 'Hom. Ref.', 'Het.', 'Hom Alt.'))
for truth in (HOM_REF_TXT, HET_TXT, HOM_ALT_TXT):
    sys.stdout.write('{:20}{:<10}{:<10}{:<10}\n'.format(
        truth, map(str, (
            count_calls[(truth, mega_call)]
            for mega_call in (HOM_REF_TXT, HET_TXT, HOM_ALT_TXT)))))
