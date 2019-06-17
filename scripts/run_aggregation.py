#!/usr/bin/env python3
import os
# set blas library environment variables (without these the cblas calls
# can completely halt processing)
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import argparse
from time import sleep

from megalodon import (
    aggregate, backends, logging, mapping, mods, snps, megalodon_helper as mh)


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--taiyaki-model-filename',
        help='Taiyaki model checkpoint file. Default: Load default model ' +
        '({})'.format(mh.MODEL_PRESET_DESC))
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
        '--output-suffix', default='re_aggregated',
        help='Suffix to apply to aggregated results, to avoid ' +
        'overwriting results. Default: %(default)s')
    parser.add_argument(
        '--processes', type=int, default=1,
        help='Number of parallel processes. Default: %(default)d')
    parser.add_argument(
        '--read-ids-filename',
        help='File containing read ids to process (one per ' +
        'line). Default: All reads')
    parser.add_argument(
        '--reference',
        help='Reference FASTA or minimap2 index file used for mapping ' +
        'called reads. Used to annotate VCF file with contig names.')
    parser.add_argument(
        '--suppress-progress', action='store_true',
        help='Suppress progress bar output.')
    parser.add_argument(
        '--write-vcf-log-probs', action='store_true',
        help='Write alt log prbabilities out in non-standard VCF field.')

    return parser

def main():
    args = get_parser().parse_args()
    log_suffix = ('aggregation' if args.output_suffix is None else
                  'aggregation.' + args.output_suffix)
    logging.init_logger(args.output_directory, out_suffix=log_suffix)
    logger = logging.get_logger()
    tai_model_fn = mh.get_model_fn(args.taiyaki_model_filename)
    model_info = backends.ModelInfo(tai_model_fn)
    mod_names = (model_info.mod_long_names
                 if mh.MOD_NAME in args.outputs else [])
    mod_agg_info = mods.AGG_INFO(
        mods.BIN_THRESH_NAME, args.mod_binary_threshold)
    aligner = mapping.alignerPlus(
        str(args.reference), preset=str('map-ont'), best_n=1)
    aligner.add_ref_lens()
    valid_read_ids = None
    if args.read_ids_filename is not None:
        with open(args.read_ids_filename) as read_ids_fp:
            valid_read_ids = set(line.strip() for line in read_ids_fp)
    aggregate.aggregate_stats(
        args.outputs, args.output_directory, args.processes,
        args.write_vcf_log_probs, args.heterozygous_factors,
        snps.HAPLIOD_MODE if args.haploid else snps.DIPLOID_MODE,
        mod_names, mod_agg_info, args.suppress_progress,
        aligner.ref_names_and_lens, valid_read_ids, args.output_suffix)

    if mh.SNP_NAME in args.outputs:
        logger.info('Sorting output variant file')
        variant_fn = mh.add_fn_suffix(
            mh.get_megalodon_fn(args.output_directory, mh.SNP_NAME),
            args.output_suffix)
        sort_variant_fn = mh.add_fn_suffix(variant_fn, 'sorted')
        snps.sort_variants(variant_fn, sort_variant_fn)
        logger.info('Indexing output variant file')
        index_var_fn = snps.index_variants(sort_variant_fn)

    return


if __name__ == '__main__':
    main()
