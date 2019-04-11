import os
import sys
import argparse
from collections import defaultdict

import numpy as np
from tqdm import tqdm


DO_PLOT = False
DEFAULT_SMOOTH_BW = 0.8
DEFAULT_SMOOTH_MAX = 200
DEFAULT_SMOOTH_NVALS = 1001
DEFAULT_MIN_DENSITY = 1e-5

if DO_PLOT:
    import matplotlib
    if sys.platform == 'darwin':
        matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt
    #plt.ion()


def determine_min_dens_edge(
        sm_ref, sm_alt, num_calib_vals, min_dens_val, smooth_ls):
    """ Compute positions where the density values are too small to produce
    robust calibration estimates and return range with valid density values
    """
    lower_invalid_dens_pos = 0
    ref_before = np.where(sm_ref[:num_calib_vals // 2] < min_dens_val)[0]
    if len(ref_before) > 0:
        lower_invalid_dens_pos = max(ref_before[-1], lower_invalid_dens_pos)
    alt_before = np.where(sm_alt[:num_calib_vals // 2] < min_dens_val)[0]
    if len(alt_before) > 0:
        lower_invalid_dens_pos = max(alt_before[-1], lower_invalid_dens_pos)

    upper_invalid_dens_pos = 1
    ref_after = np.where(sm_ref[num_calib_vals // 2:] < min_dens_val)[0]
    if len(ref_after) > 0:
        upper_invalid_dens_pos = max((num_calib_vals // 2) - ref_after[0] + 1,
                                     upper_invalid_dens_pos)
    alt_after = np.where(sm_alt[num_calib_vals // 2:] < min_dens_val)[0]
    if len(alt_after) > 0:
        upper_invalid_dens_pos = max((num_calib_vals // 2) - alt_after[0] + 1,
                                     upper_invalid_dens_pos)

    return (np.around(smooth_ls[lower_invalid_dens_pos]).astype(int),
            np.around(smooth_ls[-upper_invalid_dens_pos]).astype(int))



def compute_smooth_mono_density(llrs, num_calib_vals, smooth_bw, smooth_ls):
    def guassian(x):
        return (np.exp(-x ** 2 / (2 * smooth_bw ** 2)) /
                (smooth_bw * np.sqrt(2 * np.pi)))


    smooth_vals = np.zeros(num_calib_vals)
    for llr in tqdm(llrs, smoothing=0):
        smooth_vals += guassian(smooth_ls - llr)
    smooth_vals /= llrs.shape[0]

    peak_site = np.argmax(smooth_vals)
    # force monotonic increasing before peak and monotonic decreasing after
    mono_smooth_vals = np.concatenate([
        np.mean(np.stack([
            np.maximum.accumulate(smooth_vals[:peak_site]),
            np.minimum.accumulate(smooth_vals[:peak_site][::-1])[::-1]]),
                axis=0),
        np.mean(np.stack([
            np.minimum.accumulate(smooth_vals[peak_site:]),
            np.maximum.accumulate(smooth_vals[peak_site:][::-1])[::-1]]),
                axis=0)])

    return mono_smooth_vals, smooth_vals


def compute_calibration(
        ref_llrs, alt_llrs, max_input_llr, num_calib_vals, smooth_bw,
        min_dens_val):
    smooth_ls = np.linspace(-max_input_llr, max_input_llr,
                            num_calib_vals, endpoint=True)
    sys.stderr.write('\tComputing reference emperical density.\n')
    sm_ref, s_ref = compute_smooth_mono_density(
        ref_llrs, num_calib_vals, smooth_bw, smooth_ls)
    sys.stderr.write('\tComputing alternative emperical density.\n')
    sm_alt, s_alt = compute_smooth_mono_density(
        alt_llrs, num_calib_vals, smooth_bw, smooth_ls)

    # the ratio of very small density values can cause invalid or inaccurate
    # calibration values, so check that the max_input_llr is valid or
    # find a valid clipping location according to min_dens_val
    # then recompute smooth values
    new_input_llr_range = determine_min_dens_edge(
        sm_ref, sm_alt, num_calib_vals, min_dens_val, smooth_ls)
    if (new_input_llr_range[0] != -max_input_llr or
        new_input_llr_range[1] != max_input_llr):
        sys.stderr.write(
            '\tSetting new input llr range for more robust calibration ' +
            '({}, {})\n'.format(*new_input_llr_range))
        smooth_ls = np.linspace(new_input_llr_range[0], new_input_llr_range[1],
                                num_calib_vals, endpoint=True)
        sys.stderr.write('\tComputing new reference emperical density.\n')
        sm_ref, s_ref = compute_smooth_mono_density(
            ref_llrs, num_calib_vals, smooth_bw, smooth_ls)
        sys.stderr.write('\tComputing new alternative emperical density.\n')
        sm_alt, s_alt = compute_smooth_mono_density(
            alt_llrs, num_calib_vals, smooth_bw, smooth_ls)


    prob_alt = sm_alt / (sm_ref + sm_alt)
    # compute probability mid-point
    prob_mp = np.argmin(np.abs(prob_alt - 0.5))
    # force monotonic decreasing with reverse maximum before p=0.5 and
    # forward minimum after p=0.5
    mono_prob = np.concatenate([
        np.maximum.accumulate(prob_alt[:prob_mp][::-1])[::-1],
        np.minimum.accumulate(prob_alt[prob_mp:])])

    if DO_PLOT:
        sys.stderr.write('\tPlotting.\n')
        f, axarr = plt.subplots(3, sharex=True)
        axarr[0].plot(smooth_ls, s_ref, color='orange')
        axarr[0].plot(smooth_ls, sm_ref, color='red')
        axarr[0].plot(smooth_ls, s_alt, color='grey')
        axarr[0].plot(smooth_ls, sm_alt, color='blue')
        axarr[1].plot(smooth_ls, mono_prob, color='orange')
        axarr[1].plot(smooth_ls, 1 / (np.exp(smooth_ls) + 1), color='purple')
        axarr[2].plot(smooth_ls, np.log((1 - prob_alt) / prob_alt), color='red')
        axarr[2].plot(smooth_ls, np.log((1 - mono_prob) / mono_prob),
                      color='orange')
        plt.show()

    return np.log((1 - mono_prob) / mono_prob), new_input_llr_range


def extract_llrs(llr_fn, max_indel_len=None):
    snp_llrs, ins_ref_llrs, ins_alt_llrs, del_ref_llrs, del_alt_llrs = (
        [], [], [], [], [])
    with open(llr_fn) as llr_fp:
        for line in llr_fp:
            is_ref_correct, llr, ref_seq, alt_seq = line.split()
            llr = float(llr)
            if np.isnan(llr): continue
            if (max_indel_len is not None and
                np.abs(len(ref_seq) - len(alt_seq)) > max_indel_len):
                continue
            if len(ref_seq) == 1 and len(alt_seq) == 1:
                snp_llrs.append(llr)
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

    return map(np.array, (snp_llrs, ins_ref_llrs, ins_alt_llrs,
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
        'generate_ground_truth_indel_scores.py). Default: %(default)s')
    parser.add_argument(
        '--max-input-llr', type=int, default=DEFAULT_SMOOTH_MAX,
        help='Maximum log-likelihood ratio to compute calibration. ' +
        'Default: %(default)d')
    parser.add_argument(
        '--num-calibration-values', type=int, default=DEFAULT_SMOOTH_NVALS,
        help='Number of discrete calibration values to compute. ' +
        'Default: %(default)d')
    parser.add_argument(
        '--smooth-bandwidth', type=float, default=DEFAULT_SMOOTH_BW,
        help='Smoothing bandwidth. Default: %(default)f')
    parser.add_argument(
        '--min-density', type=float, default=DEFAULT_MIN_DENSITY,
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
    (snp_llrs, ins_ref_llrs, ins_alt_llrs,
     del_ref_llrs, del_alt_llrs) = extract_llrs(args.ground_truth_llrs)

    sys.stderr.write('Computing single-base SNP calibration.\n')
    snp_calib, snp_llr_range = compute_calibration(
        snp_llrs, -snp_llrs, args.max_input_llr, args.num_calibration_values,
        args.smooth_bandwidth, args.min_density)
    sys.stderr.write('Computing deletion calibration.\n')
    del_calib, del_llr_range = compute_calibration(
        del_ref_llrs, del_alt_llrs, args.max_input_llr,
        args.num_calibration_values, args.smooth_bandwidth, args.min_density)
    sys.stderr.write('Computing insertion calibration.\n')
    ins_calib, ins_llr_range = compute_calibration(
        ins_ref_llrs, ins_alt_llrs, args.max_input_llr,
        args.num_calibration_values, args.smooth_bandwidth, args.min_density)

    # save valibration table for reading into SNP table
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
