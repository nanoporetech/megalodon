import sys

import numpy as np
from tqdm import tqdm

from megalodon import megalodon_helper as mh


DEFAULT_SMOOTH_BW = 0.8
DEFAULT_SMOOTH_MAX = 200
DEFAULT_SMOOTH_NVALS = 1001
DEFAULT_MIN_DENSITY = 5e-6

VAR_CALIB_TYPE = 'snp_type_indel_len'
GENERIC_BASE = 'N'
SNP_CALIB_TMPLT = 'snp_{}_{}_calibration'
SNP_LLR_RNG_TMPLT = 'snp_{}_{}_llr_range'
DEL_CALIB_TMPLT = 'del_{}_calibration'
DEL_LLR_RNG_TMPLT = 'del_{}_llr_range'
INS_CALIB_TMPLT = 'ins_{}_calibration'
INS_LLR_RNG_TMPLT = 'ins_{}_llr_range'


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
    for llr in tqdm(llrs, smoothing=0, dynamic_ncols=True):
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


def compute_mirrored_calibration(
        ref_llrs, max_input_llr, num_calib_vals, smooth_bw,
        min_dens_val, return_plot_info=False):
    smooth_ls = np.linspace(-max_input_llr, max_input_llr,
                            num_calib_vals, endpoint=True)

    sys.stderr.write('\tComputing reference emperical density.\n')
    sm_ref, s_ref = compute_smooth_mono_density(
        ref_llrs, num_calib_vals, smooth_bw, smooth_ls)

    # the ratio of very small density values can cause invalid or inaccurate
    # calibration values, so check that the max_input_llr is valid or
    # find a valid clipping location according to min_dens_val
    # then recompute smooth values
    new_input_llr_range = determine_min_dens_edge(
        sm_ref, sm_ref[::-1], num_calib_vals, min_dens_val, smooth_ls)
    assert new_input_llr_range[0] == -new_input_llr_range[1]
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


    prob_alt = sm_ref[::-1] / (sm_ref + sm_ref[::-1])
    # compute probability mid-point (llr=0 for mirrored)
    prob_mp = int(np.around(num_calib_vals / 2))
    # force monotonic decreasing with reverse maximum before p=0.5 and
    # forward minimum after p=0.5
    mono_prob = np.concatenate([
        np.maximum.accumulate(prob_alt[:prob_mp][::-1])[::-1],
        np.minimum.accumulate(prob_alt[prob_mp:])])

    plot_data = None
    if return_plot_info:
        plot_data = (smooth_ls, s_ref, sm_ref, s_ref[::-1], sm_ref[::-1],
                     mono_prob, prob_alt)

    return np.log((1 - mono_prob) / mono_prob), new_input_llr_range, plot_data


#####################
##### LLR Stats #####
#####################

def compute_log_probs(alt_llrs):
    """ Compute log probabilities from a set of log likelihood ratios all
    against the reference allele
    """
    ref_lp = np.log(1) - np.log1p(np.sum(1 / np.exp(alt_llrs)))
    # set maximum log probability to avoid reporting 0 and 1 probabilities
    return ref_lp - alt_llrs


###############################
##### Calibration Readers #####
###############################

class VarCalibrator(object):
    def _load_calibration(self):
        calib_data = np.load(self.fn)
        self.stratify_type = str(calib_data['stratify_type'])
        assert self.stratify_type == VAR_CALIB_TYPE

        self.num_calib_vals = np.int(calib_data['smooth_nvals'])
        self.max_indel_len = np.int(calib_data['max_indel_len'])

        (self.snp_llr_ranges, self.snp_steps, self.snp_calib_tables,
         self.del_llr_ranges, self.del_steps, self.del_calib_tables,
         self.ins_llr_ranges, self.ins_steps, self.ins_calib_tables) = (
             {} for _ in range(9))
        # load generic snp
        ref_base, alt_base = GENERIC_BASE, GENERIC_BASE
        snp_type_llr_range = calib_data[
            SNP_LLR_RNG_TMPLT.format(ref_base, alt_base)].copy()
        self.snp_llr_ranges[(ref_base, alt_base)] = snp_type_llr_range
        self.snp_steps[(ref_base, alt_base)] = (
            snp_type_llr_range[1] - snp_type_llr_range[0]) / (
                self.num_calib_vals - 1)
        self.snp_calib_tables[(ref_base, alt_base)] = calib_data[
            SNP_CALIB_TMPLT.format(ref_base, alt_base)].copy()
        # load other base combinations
        for ref_base in mh.ALPHABET:
            for alt_base in set(mh.ALPHABET).difference(ref_base):
                snp_type_llr_range = calib_data[
                    SNP_LLR_RNG_TMPLT.format(ref_base, alt_base)].copy()
                self.snp_llr_ranges[(ref_base, alt_base)] = snp_type_llr_range
                self.snp_steps[(ref_base, alt_base)] = (
                    snp_type_llr_range[1] - snp_type_llr_range[0]) / (
                        self.num_calib_vals - 1)
                self.snp_calib_tables[(ref_base, alt_base)] = calib_data[
                    SNP_CALIB_TMPLT.format(ref_base, alt_base)].copy()
        for indel_len in range(1, self.max_indel_len + 1):
            del_type_llr_range = calib_data[
                DEL_LLR_RNG_TMPLT.format(indel_len)].copy()
            self.del_llr_ranges[indel_len] = del_type_llr_range
            self.del_steps[indel_len] = (
                del_type_llr_range[1] - del_type_llr_range[0]) / (
                    self.num_calib_vals - 1)
            self.del_calib_tables[indel_len] = calib_data[
                DEL_CALIB_TMPLT.format(indel_len)].copy()
            ins_type_llr_range = calib_data[
                INS_LLR_RNG_TMPLT.format(indel_len)].copy()
            self.ins_llr_ranges[indel_len] = ins_type_llr_range
            self.ins_steps[indel_len] = (
                ins_type_llr_range[1] - ins_type_llr_range[0]) / (
                    self.num_calib_vals - 1)
            self.ins_calib_tables[indel_len] = calib_data[
                INS_CALIB_TMPLT.format(indel_len)].copy()

        return

    def __init__(self, vars_calib_fn):
        self.fn = vars_calib_fn
        if self.fn is not None:
            self._load_calibration()
        self.calib_loaded = self.fn is not None
        return

    def calibrate_llr(self, llr, read_ref_seq, read_alt_seq):
        def simplify_var_seq(ref_seq, alt_seq):
            while (len(ref_seq) > 0 and len(alt_seq) > 0 and
                   ref_seq[0] == alt_seq[0]):
                ref_seq = ref_seq[1:]
                alt_seq = alt_seq[1:]
            while (len(ref_seq) > 0 and len(alt_seq) > 0 and
                   ref_seq[-1] == alt_seq[-1]):
                ref_seq = ref_seq[:-1]
                alt_seq = alt_seq[:-1]
            return ref_seq, alt_seq

        if not self.calib_loaded:
            return llr
        if len(read_ref_seq) == len(read_alt_seq):
            ref_seq, alt_seq = simplify_var_seq(read_ref_seq, read_alt_seq)
            # default to a "generic" SNP type that is the total of all SNP types
            snp_type = ((ref_seq, alt_seq) if (ref_seq, alt_seq)
                        in self.snp_calib_tables else
                        (GENERIC_BASE, GENERIC_BASE))
            calib_table = self.snp_calib_tables[snp_type]
            step = self.snp_steps[snp_type]
            llr_range = self.snp_llr_ranges[snp_type]
        elif len(read_ref_seq) > len(read_alt_seq):
            del_len = min(
                len(read_ref_seq) - len(read_alt_seq), self.max_indel_len)
            calib_table = self.del_calib_tables[del_len]
            step = self.del_steps[del_len]
            llr_range = self.del_llr_ranges[del_len]
        else:
            ins_len = min(
                len(read_alt_seq) - len(read_ref_seq), self.max_indel_len)
            calib_table = self.ins_calib_tables[ins_len]
            step = self.ins_steps[ins_len]
            llr_range = self.ins_llr_ranges[ins_len]

        if llr < llr_range[0]:
            llr = llr_range[0]
        elif llr > llr_range[1]:
            llr = llr_range[1]
        return calib_table[np.around((llr - llr_range[0]) / step).astype(int)]


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
        if llr < llr_range[0]:
            llr = llr_range[0]
        elif llr > llr_range[1]:
            llr = llr_range[1]
        return calib_table[np.around((llr - llr_range[0]) / step).astype(int)]


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
