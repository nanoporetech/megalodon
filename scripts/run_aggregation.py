#!/usr/bin/env python3
import os
# set blas library environment variables (without these the cblas calls
# can completely halt processing)
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import argparse

from megalodon import aggregate, backends, megalodon_helper as mh


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
        default=['basecalls',], choices=tuple(mh.OUTPUT_FNS.keys()),
        help='Output type(s) to produce. Default: %(default)s')
    parser.add_argument(
        '--heterozygous-factor', type=float, default=0.5,
        help='Bayesian prior factor for heterozygous calls (compared to 1.0 ' +
        'for hom ref/alt). Default: %(default)f')
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
        args.flappie_model_name, args.taiyaki_model_filename, None)
    mod_names = ([(mod_b, mod_b) for mod_b in model_info.alphabet[4:]]
                 if mh.MOD_NAME in args.outputs else [])
    aggregate.aggregate_stats(
        args.outputs, args.output_directory, args.processes,
        args.write_vcf_llr, args.heterozygous_factor, mod_names,
        args.suppress_progress)

    return


if __name__ == '__main__':
    main()
