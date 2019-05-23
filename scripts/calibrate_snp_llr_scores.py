import os
import sys
import argparse

import numpy as np

from megalodon import calibration


# TODO add this as a command line option as in validate script to save plots
DO_PLOT = False


def extract_llrs(llr_fn, max_indel_len=None):
    (snp_ref_llrs, snp_alt_llrs, ins_ref_llrs, ins_alt_llrs,
     del_ref_llrs, del_alt_llrs) = ([], [], [], [], [], [])
    with open(llr_fn) as llr_fp:
        for line in llr_fp:
            is_ref_correct, llr, ref_seq, alt_seq = line.split()
            llr = float(llr)
            if np.isnan(llr): continue
            if (max_indel_len is not None and
                np.abs(len(ref_seq) - len(alt_seq)) > max_indel_len):
                continue
            if len(ref_seq) == 1 and len(alt_seq) == 1:
                if is_ref_correct == 'True':
                    snp_ref_llrs.append(llr)
                else:
                    snp_alt_llrs.append(llr)
            else:
                if len(ref_seq) > len(alt_seq):
                    if is_ref_correct == 'True':
                        del_ref_llrs.append(llr)
                    else:
                        del_alt_llrs.append(llr)
                else:
                    if is_ref_correct == 'True':
                        ins_ref_llrs.append(llr)
                    else:
                        ins_alt_llrs.append(llr)

    return map(np.array, (snp_ref_llrs, snp_alt_llrs,
                          ins_ref_llrs, ins_alt_llrs,
                          del_ref_llrs, del_alt_llrs))


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
        '--ground-truth-llrs', default='snp_calibration_statistics.txt',
        help='Ground truth log-likelihood ratio statistics (produced by ' +
        'generate_ground_truth_snp_llr_scores.py). Default: %(default)s')
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
        '--out-filename', default='megalodon_snp_calibration.npz',
        help='Filename to output calibration values. Default: %(default)s')
    parser.add_argument(
        '--overwrite', action='store_true',
        help='Overwrite --out-filename if it exists.')

    return parser


def main():
    args = get_parser().parse_args()

    prep_out(args.out_filename, args.overwrite)

    sys.stderr.write('Parsing log-likelihood ratios\n')
    (snp_ref_llrs, snp_alt_llrs, ins_ref_llrs, ins_alt_llrs,
     del_ref_llrs, del_alt_llrs) = extract_llrs(args.ground_truth_llrs)

    sys.stderr.write('Computing single-base SNP calibration.\n')
    snp_calib, snp_llr_range = calibration.compute_calibration(
        snp_ref_llrs, snp_alt_llrs, args.max_input_llr,
        args.num_calibration_values, args.smooth_bandwidth, args.min_density,
        DO_PLOT)
    sys.stderr.write('Computing deletion calibration.\n')
    del_calib, del_llr_range = calibration.compute_calibration(
        del_ref_llrs, del_alt_llrs, args.max_input_llr,
        args.num_calibration_values, args.smooth_bandwidth, args.min_density,
        DO_PLOT)
    sys.stderr.write('Computing insertion calibration.\n')
    ins_calib, ins_llr_range = calibration.compute_calibration(
        ins_ref_llrs, ins_alt_llrs, args.max_input_llr,
        args.num_calibration_values, args.smooth_bandwidth, args.min_density,
        DO_PLOT)

    # save calibration table for reading into SNP table
    sys.stderr.write('Saving calibrations to file.\n')
    np.savez(
        args.out_filename,
        stratify_type='snp_ins_del',
        smooth_nvals=args.num_calibration_values,
        snp_llr_range=snp_llr_range,
        snp_calibration_table=snp_calib,
        deletion_llr_range=del_llr_range,
        deletion_calibration_table=del_calib,
        insertion_llr_range=ins_llr_range,
        insertion_calibration_table=ins_calib
    )

    return


if __name__ == '__main__':
    main()
