#!/usr/bin/env python3
import os
import sys

from megalodon import (
    aggregate, logging, mods, variants, megalodon_helper as mh)
from ._extras_parsers import get_parser_aggregate_run


# set blas library environment variables (without these the cblas calls
# can completely halt processing)
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'

LOGGER = logging.get_logger()


def _main(args):
    log_suffix = ('aggregation' if args.output_suffix is None else
                  'aggregation.' + args.output_suffix)
    logging.init_logger(args.megalodon_directory, out_suffix=log_suffix)
    LOGGER.debug('Command: """' + ' '.join(sys.argv) + '"""')

    if args.mod_aggregate_method == mh.MOD_EM_NAME:
        mod_agg_info = mods.AGG_INFO(mh.MOD_EM_NAME, None)
    elif args.mod_aggregate_method == mh.MOD_BIN_THRESH_NAME:
        mod_agg_info = mods.AGG_INFO(
            mh.MOD_BIN_THRESH_NAME, args.mod_binary_threshold)
    valid_read_ids = mh.parse_read_ids(args.read_ids_filename)
    aggregate.aggregate_stats(
        args.outputs, args.megalodon_directory, args.processes,
        args.write_vcf_log_probs, args.heterozygous_factors,
        variants.HAPLIOD_MODE if args.haploid else variants.DIPLOID_MODE,
        mod_agg_info, args.write_mod_log_probs, args.mod_output_formats,
        args.suppress_progress, valid_read_ids, args.output_suffix)

    if mh.VAR_NAME in args.outputs:
        LOGGER.info('Sorting output variant file')
        variant_fn = mh.add_fn_suffix(
            mh.get_megalodon_fn(args.megalodon_directory, mh.VAR_NAME),
            args.output_suffix)
        sort_variant_fn = mh.add_fn_suffix(variant_fn, 'sorted')
        variants.sort_variants(variant_fn, sort_variant_fn)
        LOGGER.info('Indexing output variant file')
        variants.index_variants(sort_variant_fn)


if __name__ == '__main__':
    _main(get_parser_aggregate_run().parse_args())
