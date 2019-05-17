import sys

import numpy as np
from tqdm import tqdm


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
        min_dens_val, do_plot=False):
    if do_plot:
        import matplotlib
        if sys.platform == 'darwin':
            matplotlib.use("TkAgg")
        import matplotlib.pyplot as plt
        #plt.ion()

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

    if do_plot:
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


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
