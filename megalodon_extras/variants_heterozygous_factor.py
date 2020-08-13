from collections import defaultdict

import pysam
import numpy as np

from megalodon import logging
from ._extras_parsers import get_parser_variants_heterozygous_factor


LOGGER = logging.get_logger()

HOM_REF_TXT = 'hom_ref'
HET_TXT = 'het'
HOM_ALT_TXT = 'hom_alt'

SNP_TXT = 'SNP'
DEL_TXT = 'DEL'
INS_TXT = 'INS'

STAT_WIDTH = 12
STATS_FMT_STR = '{:<' + str(STAT_WIDTH) + '}'
FLOAT_FMT_STR = '{:<' + str(STAT_WIDTH) + '.4f}'
STAT_NAMES = ('HomRef', 'Het', 'HomAlt', 'F1', 'Precision', 'Recall')
N_STATS = len(STAT_NAMES)
N_INT_STATS = 3
N_FLOAT_STATS = N_STATS - N_INT_STATS
HEADER_TMPLT = STATS_FMT_STR * (N_STATS + 1)
STATS_TMPLT = STATS_FMT_STR * (N_INT_STATS + 1) + \
    FLOAT_FMT_STR * N_FLOAT_STATS


def _main(args):
    def conv_call_str(gt_vals):
        gt_set = set(gt_vals)
        if gt_set == set([0]):
            return HOM_REF_TXT
        elif gt_set == set([0, 1]):
            return HET_TXT
        return HOM_ALT_TXT

    logging.init_logger()
    gt_calls = defaultdict(dict)
    for variant in pysam.VariantFile(args.ground_truth_variants).fetch():
        # skip mutli-allelic sites
        if variant.alts is None or len(variant.alts) > 1:
            continue
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
        if len(variant.alts) > 1:
            continue
        if len(variant.ref) == len(variant.alts[0]):
            mega_calls[SNP_TXT][(variant.contig, variant.pos, variant.ref,
                                 variant.alts[0])] = conv_call_str(
                                     variant.samples.values()[0]['GT'])
        elif len(variant.ref) > len(variant.alts[0]):
            mega_calls[DEL_TXT][(variant.contig, variant.pos, variant.ref,
                                 variant.alts[0])] = conv_call_str(
                                     variant.samples.values()[0]['GT'])
        else:
            mega_calls[INS_TXT][(variant.contig, variant.pos, variant.ref,
                                 variant.alts[0])] = conv_call_str(
                                     variant.samples.values()[0]['GT'])

    for var_type in (SNP_TXT, DEL_TXT, INS_TXT):
        counts = defaultdict(int)
        for chrm_pos_ref_alt in set(gt_calls[var_type]).intersection(
                mega_calls[var_type]):
            counts[(gt_calls[var_type][chrm_pos_ref_alt],
                    mega_calls[var_type][chrm_pos_ref_alt])] += 1

        # compute F1 stat
        vt_stats = []
        for truth_type in (HOM_REF_TXT, HET_TXT, HOM_ALT_TXT):
            gt_count = sum(
                counts[(truth_type, mega_call)]
                for mega_call in (HOM_REF_TXT, HET_TXT, HOM_ALT_TXT))
            mega_count = sum(
                counts[(gt_call, truth_type)]
                for gt_call in (HOM_REF_TXT, HET_TXT, HOM_ALT_TXT))
            if gt_count == 0 or mega_count == 0:
                vt_stats.append((np.NAN, np.NAN, np.NAN))
            else:
                prec = counts[(truth_type, truth_type)] / mega_count
                recall = counts[(truth_type, truth_type)] / gt_count
                vt_stats.append((
                    2 * (prec * recall) / (prec + recall), prec, recall))

        # print output
        LOGGER.info(var_type)
        LOGGER.info(HEADER_TMPLT.format('Truth\tCalls', *STAT_NAMES))
        for truth, (f1, prec, recall) in zip(
                (HOM_REF_TXT, HET_TXT, HOM_ALT_TXT),
                vt_stats):
            LOGGER.info(STATS_TMPLT.format(
                truth, counts[(truth, HOM_REF_TXT)], counts[(truth, HET_TXT)],
                counts[(truth, HOM_ALT_TXT)], f1, prec, recall))
        mean_f1_fmt = ('{:>' + str(STAT_WIDTH * (N_STATS - 2)) + '}' +
                       FLOAT_FMT_STR * N_FLOAT_STATS + '\n')
        mean_stats = map(np.nanmean, zip(*vt_stats))
        LOGGER.info(mean_f1_fmt.format('Mean Stats:   ', *mean_stats))


if __name__ == '__main__':
    _main(get_parser_variants_heterozygous_factor().parse_args())
