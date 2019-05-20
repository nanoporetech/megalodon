#!/usr/bin/env python3
import os
# set blas library environment variables (without these the cblas calls
# can completely halt processing)
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import argparse

from megalodon import aggregate, backends, mods, snps, megalodon_helper as mh


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--taiyaki-model-filename',
        help='Taiyaki model checkpoint file.')
    parser.add_argument(
        '--flappie-model-name',
        help='Flappie model name.')
    parser.add_argument(
        '--outputs', nargs='+',
        default=[mh.SNP_NAME, mh.MOD_NAME],
        choices=[mh.SNP_NAME, mh.MOD_NAME],
        help='Output type(s) to produce. Default: %(default)s')
    parser.add_argument(
        '--haploid', action='store_true',
        help='Compute SNP aggregation for haploid genotypes. Default: diploid')
    parser.add_argument(
        '--heterozygous-factors', type=float, nargs=2,
        default=[mh.DEFAULT_SNV_HET_FACTOR, mh.DEFAULT_INDEL_HET_FACTOR],
        help='Bayesian prior factor for snv and indel heterozygous calls ' +
        '(compared to 1.0 for hom ref/alt). Default: %(default)s')
    parser.add_argument(
        '--mod-binary-threshold', type=float, nargs=2,
        default=mods.DEFAULT_AGG_INFO.binary_threshold,
        help='Thresholds for modified base aggregation. Default: %(default)s')
    parser.add_argument(
        '--output-directory',
        default='megalodon_results',
        help='Directory to store output results. Default: %(default)s')
    parser.add_argument(
        '--processes', type=int, default=1,
        help='Number of parallel processes. Default: %(default)d')
    parser.add_argument(
        '--write-vcf-llr', action='store_true',
        help='Write log-likelihood ratios out in non-standard VCF field.')
    parser.add_argument(
        '--suppress-progress', action='store_true',
        help='Suppress progress bar output.')

    return parser

def main():
    args = get_parser().parse_args()
    model_info = backends.ModelInfo(
        args.flappie_model_name, args.taiyaki_model_filename)
    mod_names = (model_info.mod_long_names
                 if mh.MOD_NAME in args.outputs else [])
    mod_agg_info = mods.AGG_INFO(
        mods.BIN_THRESH_NAME, args.mod_binary_threshold)
    aggregate.aggregate_stats(
        args.outputs, args.output_directory, args.processes,
        args.write_vcf_llr, args.heterozygous_factors,
        snps.HAPLIOD_MODE if args.haploid else snps.DIPLOID_MODE,
        mod_names, mod_agg_info, args.suppress_progress)

    return


if __name__ == '__main__':
    main()
