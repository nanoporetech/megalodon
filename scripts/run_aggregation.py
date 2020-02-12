#!/usr/bin/env python3
import os
import sys
# set blas library environment variables (without these the cblas calls
# can completely halt processing)
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import argparse
from time import sleep

from megalodon import (
    aggregate, backends, logging, mapping, mods, variants,
    megalodon_helper as mh)


def get_parser():
    parser = argparse.ArgumentParser(
        'Aggregate per-read, per-site statistics from previous megalodon call.')

    out_grp = parser.add_argument_group('Output Arguments')
    out_grp.add_argument(
        '--outputs', nargs='+',
        default=[mh.VAR_NAME, mh.MOD_NAME],
        choices=[mh.VAR_NAME, mh.MOD_NAME],
        help='Output type(s) to produce. Default: %(default)s')
    out_grp.add_argument(
        '--megalodon-directory',
        default='megalodon_results',
        help='Megalodon output directory containing per-read database(s) ' +
        'where aggregated results will be added. Default: %(default)s')
    out_grp.add_argument(
        '--output-suffix', default='re_aggregated',
        help='Suffix to apply to aggregated results, to avoid ' +
        'overwriting results. Default: %(default)s')
    out_grp.add_argument(
        '--read-ids-filename',
        help='File containing read ids to process (one per ' +
        'line). Default: All reads')

    var_grp = parser.add_argument_group('Sequence Variant Arguments')
    var_grp.add_argument(
        '--haploid', action='store_true',
        help='Compute sequence variant aggregation for haploid genotypes. ' +
        'Default: diploid')
    var_grp.add_argument(
        '--heterozygous-factors', type=float, nargs=2,
        default=[mh.DEFAULT_SNV_HET_FACTOR, mh.DEFAULT_INDEL_HET_FACTOR],
        help='Bayesian prior factor for snv and indel heterozygous calls ' +
        '(compared to 1.0 for hom ref/alt). Default: %(default)s')
    var_grp.add_argument(
        '--write-vcf-log-probs', action='store_true',
        help='Write alt log prbabilities out in non-standard VCF field.')

    mod_grp = parser.add_argument_group('Modified Base Arguments')
    mod_grp.add_argument(
        '--mod-aggregate-method', choices=list(mods.AGG_METHOD_NAMES),
        default=mods.BIN_THRESH_NAME,
        help='Modified base aggregation method. Default: %(default)s')
    mod_grp.add_argument(
        '--mod-binary-threshold', type=float, nargs=1,
        default=mods.DEFAULT_BINARY_THRESH,
        help='Threshold for modified base aggregation (probability of ' +
        'modified/canonical base). Only applicable for ' +
        '"--mod-aggregate-method binary_threshold". Default: %(default)s')
    mod_grp.add_argument(
        '--mod-output-formats', nargs='+',
        default=[mh.MOD_BEDMETHYL_NAME,],
        choices=tuple(mh.MOD_OUTPUT_FMTS.keys()),
        help='Modified base aggregated output format(s). Default: %(default)s')
    mod_grp.add_argument(
        '--write-mod-log-probs', action='store_true',
        help='Write per-read modified base log probabilities ' +
        'out in non-standard modVCF field.')

    misc_grp = parser.add_argument_group('Miscellaneous Arguments')
    misc_grp.add_argument(
        '--processes', type=int, default=1,
        help='Number of parallel processes. Default: %(default)d')
    misc_grp.add_argument(
        '--suppress-progress', action='store_true',
        help='Suppress progress bar output.')

    return parser

def main():
    args = get_parser().parse_args()
    log_suffix = ('aggregation' if args.output_suffix is None else
                  'aggregation.' + args.output_suffix)
    logging.init_logger(args.megalodon_directory, out_suffix=log_suffix)
    logger = logging.get_logger()
    logger.debug('Command: """' + ' '.join(sys.argv) + '"""')

    if args.mod_aggregate_method == mods.EM_NAME:
        mod_agg_info = mods.AGG_INFO(mods.EM_NAME, None)
    elif args.mod_aggregate_method == mods.BIN_THRESH_NAME:
        mod_agg_info = mods.AGG_INFO(
            mods.BIN_THRESH_NAME, args.mod_binary_threshold)
    valid_read_ids = None
    if args.read_ids_filename is not None:
        with open(args.read_ids_filename) as read_ids_fp:
            valid_read_ids = set(line.strip() for line in read_ids_fp)
    aggregate.aggregate_stats(
        args.outputs, args.megalodon_directory, args.processes,
        args.write_vcf_log_probs, args.heterozygous_factors,
        variants.HAPLIOD_MODE if args.haploid else variants.DIPLOID_MODE,
        mod_agg_info, args.write_mod_log_probs, args.mod_output_formats,
        args.suppress_progress, valid_read_ids, args.output_suffix)

    if mh.VAR_NAME in args.outputs:
        logger.info('Sorting output variant file')
        variant_fn = mh.add_fn_suffix(
            mh.get_megalodon_fn(args.megalodon_directory, mh.VAR_NAME),
            args.output_suffix)
        sort_variant_fn = mh.add_fn_suffix(variant_fn, 'sorted')
        variants.sort_variants(variant_fn, sort_variant_fn)
        logger.info('Indexing output variant file')
        index_var_fn = variants.index_variants(sort_variant_fn)

    return


if __name__ == '__main__':
    main()
