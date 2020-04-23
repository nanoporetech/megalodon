import os
import sys
import argparse
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from megalodon import calibration, mods


PROB_COLORS = ("#bcbddc", "#807dba", "#6a51a3")


def plot_calib(
        pdf_fp, mod_base, smooth_ls, s_ref, sm_ref, s_alt, sm_alt,
        mono_prob, prob_alt, prob_threshs, add_prob_thresh):
    f, axarr = plt.subplots(3, sharex=True, figsize=(11, 7))
    axarr[0].plot(smooth_ls, s_ref, color='orange')
    axarr[0].plot(smooth_ls, sm_ref, color='red')
    axarr[0].plot(smooth_ls, s_alt, color='grey')
    axarr[0].plot(smooth_ls, sm_alt, color='blue')
    axarr[0].set_ylabel(
        'Probability Density\nred/orange=canonical\nblue/grey=modified')
    axarr[0].set_title(mod_base + ' Calibration')
    axarr[1].plot(smooth_ls, mono_prob, color='orange')
    axarr[1].plot(
        smooth_ls, 1 / (np.exp(smooth_ls) + 1), color='purple')
    axarr[1].set_ylabel(
        'Emperical Modified\nProbability\norange=calibrated\npurple=raw')
    axarr[2].plot(
        smooth_ls, np.log((1 - prob_alt) / prob_alt), color='red')
    axarr[2].plot(
        smooth_ls, np.log((1 - mono_prob) / mono_prob), color='orange')
    axarr[2].set_ylabel('Calibrated LLR\norage=monotonic')
    axarr[2].set_xlabel('Theoretical LLR (NN Score)')
    if add_prob_thresh:
        # indicate the cutoff points for several common cutoff locations
        thresh_f = np.log((1 - mono_prob) / mono_prob)
        for p, col in zip(prob_threshs, PROB_COLORS):
            llr_x = np.log(p / (1 - p))
            thresh_val = np.argmin(np.abs(thresh_f - llr_x))
            nthresh_val = np.argmin(np.abs(thresh_f + llr_x))
            for i in range(2):
                axarr[i].axvline(x=smooth_ls[thresh_val], color=col)
                axarr[i].axvline(x=smooth_ls[nthresh_val], color=col)
            axarr[2].axvline(x=smooth_ls[thresh_val], color=col)
            axarr[2].axvline(x=smooth_ls[nthresh_val], color=col, label=p)
            axarr[2].legend(loc='upper right', bbox_to_anchor=(1, -0.12),
                            ncol=3)

    pdf_fp.savefig(bbox_inches='tight')
    plt.close()
    return


def extract_llrs(llr_fn):
    llrs_data = np.load(llr_fn)
    mod_bases = llrs_data[mods.GT_ALL_MOD_BASE_STR]
    mod_base_llrs = {}
    for mod_base in mod_bases:
        mod_base_llrs[mod_base] = (
            llrs_data[mods.GT_MOD_LLR_STR.format(mod_base)],
            llrs_data[mods.GT_CAN_LLR_STR.format(mod_base)])

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
    except Exception:
        sys.stderr.write(
            '*' * 60 + '\nERROR: Attempt to write to --out-filename ' +
            'location failed with the following error.\n' + '*' * 60 + '\n\n')
        raise

    return


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--ground-truth-llrs', default='mod_calibration_statistics.npz',
        help='Ground truth log-likelihood ratio statistics (produced by ' +
        'generate_ground_truth_mod_llr_scores.py). Default: %(default)s')
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
        '--smooth-bandwidth', type=float,
        default=calibration.DEFAULT_SMOOTH_BW,
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
        '--out-pdf',
        help='Output pdf filename for modified base calibration ' +
        'visualization. Default: Do not produce plot.')
    parser.add_argument(
        '--pdf-prob-thresholds', nargs=3, type=float, default=[0.75, 0.8, 0.9],
        help='Probability thresholds to mark on output pdf.')
    parser.add_argument(
        '--plot-without-prob-thresholds', action='store_true',
        help='Do not include probability thresholds in plot(s).')
    parser.add_argument(
        '--overwrite', action='store_true',
        help='Overwrite --out-filename if it exists.')

    return parser


def main():
    args = get_parser().parse_args()

    prep_out(args.out_filename, args.overwrite)

    sys.stderr.write('Parsing log-likelihood ratios\n')
    mod_base_llrs = extract_llrs(args.ground_truth_llrs)

    pdf_fp = None if args.out_pdf is None else PdfPages(args.out_pdf)
    save_kwargs = {}
    for mod_base, (mod_llrs, can_llrs) in mod_base_llrs.items():
        sys.stderr.write(
            'Computing {} modified base calibration.\n'.format(mod_base))
        mod_calib, mod_llr_range, plot_data = calibration.compute_calibration(
            can_llrs, mod_llrs, args.max_input_llr,
            args.num_calibration_values, args.smooth_bandwidth,
            args.min_density, pdf_fp is not None)
        save_kwargs[mod_base + '_llr_range'] = mod_llr_range
        save_kwargs[mod_base + '_calibration_table'] = mod_calib
        if pdf_fp is not None:
            plot_calib(pdf_fp, mod_base, *plot_data, args.pdf_prob_thresholds,
                       not args.plot_without_prob_thresholds)
    if pdf_fp is not None:
        pdf_fp.close()

    # save calibration table for reading into mod calibration table
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
