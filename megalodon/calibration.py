import sys

import numpy as np
from tqdm import tqdm


DEFAULT_SMOOTH_BW = 0.8
DEFAULT_SMOOTH_MAX = 200
DEFAULT_SMOOTH_NVALS = 1001
DEFAULT_MIN_DENSITY = 1e-5


##################################
##### Calibration Estimation #####
##################################

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
        min_dens_val, return_plot_info=False):
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

    plot_data = None
    if return_plot_info:
        plot_data = (smooth_ls, s_ref, sm_ref, s_alt, sm_alt, mono_prob,
                     prob_alt)

    return np.log((1 - mono_prob) / mono_prob), new_input_llr_range, plot_data


###############################
##### Calibration Readers #####
###############################

class SnpCalibrator(object):
    def _load_calibration(self):
        calib_data = np.load(self.fn)
        self.stratify_type = str(calib_data['stratify_type'])
        assert self.stratify_type == 'snp_ins_del'

        self.num_calib_vals = np.int(calib_data['smooth_nvals'])

        self.snp_llr_range = calib_data['snp_llr_range'].copy()
        self.snp_step = (self.snp_llr_range[1] - self.snp_llr_range[0]) / (
            self.num_calib_vals - 1)
        self.snp_calib_table = calib_data['snp_calibration_table'].copy()

        self.del_llr_range = calib_data['deletion_llr_range'].copy()
        self.del_step = (self.del_llr_range[1] - self.del_llr_range[0]) / (
            self.num_calib_vals - 1)
        self.del_calib_table = calib_data['deletion_calibration_table'].copy()

        self.ins_llr_range = calib_data['insertion_llr_range'].copy()
        self.ins_step = (self.ins_llr_range[1] - self.ins_llr_range[0]) / (
            self.num_calib_vals - 1)
        self.ins_calib_table = calib_data['insertion_calibration_table'].copy()

        return

    def __init__(self, snps_calib_fn):
        self.fn = snps_calib_fn
        if self.fn is not None:
            self._load_calibration()
        self.calib_loaded = self.fn is not None
        return

    def calibrate_llr(self, llr, snp_ref_seq, snp_alt_seq):
        if not self.calib_loaded:
            return llr
        if len(snp_ref_seq) == len(snp_alt_seq):
            return self.snp_calib_table[np.around((
                np.clip(llr, self.snp_llr_range[0], self.snp_llr_range[1]) -
                self.snp_llr_range[0]) / self.snp_step).astype(int)]
        elif len(snp_ref_seq) > len(snp_alt_seq):
            return self.del_calib_table[np.around((
                np.clip(llr, self.del_llr_range[0], self.del_llr_range[1]) -
                self.del_llr_range[0]) / self.del_step).astype(int)]
        return self.ins_calib_table[np.around((
            np.clip(llr, self.ins_llr_range[0], self.ins_llr_range[1]) -
            self.ins_llr_range[0]) / self.ins_step).astype(int)]

class ModCalibrator(object):
    def _load_calibration(self):
        calib_data = np.load(self.fn)
        self.stratify_type = str(calib_data['stratify_type'])
        assert self.stratify_type == 'mod_base'

        self.num_calib_vals = np.int(calib_data['smooth_nvals'])
        self.mod_bases = calib_data['mod_bases']
        self.mod_base_calibs = {}
        for mod_base in self.mod_bases:
            mod_llr_range = calib_data[mod_base + '_llr_range'].copy()
            mod_step = (mod_llr_range[1] - mod_llr_range[0]) / (
                self.num_calib_vals - 1)
            mod_calib_table = calib_data[mod_base + '_calibration_table'].copy()
            self.mod_base_calibs[mod_base] = (
                mod_llr_range, mod_step, mod_calib_table)

        return

    def __init__(self, mods_calib_fn):
        self.fn = mods_calib_fn
        if self.fn is not None:
            self._load_calibration()
        self.calib_loaded = self.fn is not None
        return

    def calibrate_llr(self, llr, mod_base):
        if not self.calib_loaded or mod_base not in self.mod_base_calibs:
            return llr

        llr_range, step, calib_table = self.mod_base_calibs[mod_base]
        return calib_table[np.around((
            np.clip(llr, llr_range[0], llr_range[1]) -
            llr_range[0]) / step).astype(int)]


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
