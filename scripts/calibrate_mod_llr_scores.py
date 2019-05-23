import os
import sys
import argparse
from collections import defaultdict

import numpy as np

from megalodon import calibration


# TODO add this as a command line option as in validate script to save plots
DO_PLOT = False


def extract_llrs(llr_fn):
    mod_base_llrs = defaultdict(lambda: ([], []))
    with open(llr_fn) as llr_fp:
        for line in llr_fp:
            is_mod, llr, mod_type = line.split()
            llr = float(llr)
            if np.isnan(llr): continue
            if is_mod == 'True':
                mod_base_llrs[mod_type][0].append(llr)
            else:
                mod_base_llrs[mod_type][1].append(llr)

    return mod_base_llrs


def prep_out(out_fn, overwrite):
    if os.path.exists(out_fn):
        if overwrite:
            os.remove(out_fn)
        else:
            raise NotImplementedError(
                'ERROR: --out-filename exists and --overwrite not set.')
    try:
        open(out_fn, 'w').close()
        os.remove(out_fn)
    except:
        sys.stderr.write(
            '*' * 60 + '\nERROR: Attempt to write to --out-filename location ' +
            'failed with the following error.\n' + '*' * 60 + '\n\n')
        raise

    return


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--ground-truth-llrs', default='mod_calibration_statistics.txt',
        help='Ground truth log-likelihood ratio statistics (produced by ' +
        'generate_ground_truth_snp_mod_scores.py). Default: %(default)s')
    parser.add_argument(
        '--max-input-llr', type=int, default=calibration.DEFAULT_SMOOTH_MAX,
        help='Maximum log-likelihood ratio to compute calibration. ' +
        'Default: %(default)d')
    parser.add_argument(
        '--num-calibration-values', type=int,
        default=calibration.DEFAULT_SMOOTH_NVALS,
        help='Number of discrete calibration values to compute. ' +
        'Default: %(default)d')
    parser.add_argument(
        '--smooth-bandwidth', type=float, default=calibration.DEFAULT_SMOOTH_BW,
        help='Smoothing bandwidth. Default: %(default)f')
    parser.add_argument(
        '--min-density', type=float, default=calibration.DEFAULT_MIN_DENSITY,
        help='Minimum density value to compute calibration. This value ' +
        'dynamically adjusts [--max-input-llr] when it is too large. ' +
        'Default: %(default)f')
    parser.add_argument(
        '--out-filename', default='megalodon_mod_calibration.npz',
        help='Filename to output calibration values. Default: %(default)s')
    parser.add_argument(
        '--overwrite', action='store_true',
        help='Overwrite --out-filename if it exists.')

    return parser


def main():
    args = get_parser().parse_args()

    prep_out(args.out_filename, args.overwrite)

    sys.stderr.write('Parsing log-likelihood ratios\n')
    mod_base_llrs = extract_llrs(args.ground_truth_llrs)

    save_kwargs = {}
    for mod_type, (mod_llrs, can_llrs) in mod_base_llrs.items():
        sys.stderr.write(
            'Computing {} modified base calibration.\n'.format(mod_type))
        mod_calib, mod_llr_range = calibration.compute_calibration(
            np.array(can_llrs), np.array(mod_llrs), args.max_input_llr,
            args.num_calibration_values, args.smooth_bandwidth,
            args.min_density, DO_PLOT)
        save_kwargs[mod_type + '_llr_range'] = mod_llr_range
        save_kwargs[mod_type + '_calibration_table'] = mod_calib

    # save calibration table for reading into SNP table
    sys.stderr.write('Saving calibrations to file.\n')
    mod_bases = list(mod_base_llrs.keys())
    np.savez(
        args.out_filename,
        stratify_type='mod_base',
        smooth_nvals=args.num_calibration_values,
        mod_bases=mod_bases,
        **save_kwargs)

    return


if __name__ == '__main__':
    main()
