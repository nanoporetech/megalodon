import sys
import queue
from time import sleep
import multiprocessing as mp

import numpy as np
from tqdm import tqdm

from megalodon import logging, megalodon_helper as mh


LOGGER = logging.get_logger()

VAR_CALIB_TYPE = 'snp_type_indel_len'
GENERIC_BASE = 'N'
SNP_CALIB_TMPLT = 'snp_{}_{}_calibration'
SNP_LLR_RNG_TMPLT = 'snp_{}_{}_llr_range'
DEL_CALIB_TMPLT = 'del_{}_calibration'
DEL_LLR_RNG_TMPLT = 'del_{}_llr_range'
INS_CALIB_TMPLT = 'ins_{}_calibration'
INS_LLR_RNG_TMPLT = 'ins_{}_llr_range'

# modified base fixed text strings
MOD_STRAT_TYPE_TXT = 'stratify_type'
MOD_BASE_STRAT_TYPE = 'mod_base'
MOD_VALID_STRAT_TYPES = set((MOD_BASE_STRAT_TYPE, ))
SMOOTH_NVALS_TXT = 'smooth_nvals'
MOD_BASES_TXT = 'mod_bases'
LLR_RANGE_SUFFIX = '_llr_range'
CALIB_TABLE_SUFFIX = '_calibration_table'

DEFAULT_SMOOTH_MP_BATCH_SIZE = 10000


##########################
# Calibration Estimation #
##########################

def determine_llr_plateau_edge(
        sm_ref, sm_alt, num_calib_vals, diff_eps, llr_buffer, smooth_ls):
    """ Compute new edges of calibration computation based on sites where
    log likelihood ratios plateau.
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        prob_alt = sm_alt / (sm_ref + sm_alt)
    # compute probability mid-point (llr=0 for mirrored)
    prob_mp = int(np.around(num_calib_vals / 2))
    # force monotonic decreasing with reverse maximum before p=0.5 and
    # forward minimum after p=0.5
    mono_prob = np.concatenate([
        np.maximum.accumulate(prob_alt[:prob_mp][::-1])[::-1],
        np.minimum.accumulate(prob_alt[prob_mp:])])
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        llr = np.log((1 - mono_prob) / mono_prob)
        llr_moved_sites = np.where(np.diff(llr) > diff_eps)[0]
    if len(llr_moved_sites) == 0:
        return 0, 0
    return (np.around(smooth_ls[llr_moved_sites[0]]).astype(int) - llr_buffer,
            np.around(smooth_ls[llr_moved_sites[-1]]).astype(int) + llr_buffer)


def determine_min_dens_edge(sm_ref, sm_alt, min_dens_val, smooth_ls):
    """ Compute positions where the density values are too small to produce
    robust calibration estimates and return range with valid density values.

    Assumes input densities are smooth and monotonic decreasing from a
    single peak.
    """
    ref_peak = np.argmax(sm_ref)
    alt_peak = np.argmax(sm_alt)
    lower_invalid_dens_pos = 0
    ref_before = np.where(sm_ref[:ref_peak] < min_dens_val)[0]
    if len(ref_before) > 0:
        lower_invalid_dens_pos = max(ref_before[-1], lower_invalid_dens_pos)
    alt_before = np.where(sm_alt[:alt_peak] < min_dens_val)[0]
    if len(alt_before) > 0:
        lower_invalid_dens_pos = max(alt_before[-1], lower_invalid_dens_pos)

    upper_invalid_dens_pos = 1
    ref_after = np.where(sm_ref[ref_peak:] < min_dens_val)[0]
    if len(ref_after) > 0:
        upper_invalid_dens_pos = max(
            sm_ref.shape[0] - ref_peak - ref_after[0] + 1,
            upper_invalid_dens_pos)
    alt_after = np.where(sm_alt[alt_peak:] < min_dens_val)[0]
    if len(alt_after) > 0:
        upper_invalid_dens_pos = max(
            sm_alt.shape[0] - alt_peak - alt_after[0] + 1,
            upper_invalid_dens_pos)

    return (np.around(smooth_ls[lower_invalid_dens_pos]).astype(int),
            np.around(smooth_ls[-upper_invalid_dens_pos]).astype(int))


def _compute_smooth_density_worker(llr_q, smooth_llr_q, smooth_bw, smooth_ls):
    def guassian(x):
        return (np.exp(-x ** 2 / (2 * smooth_bw ** 2)) /
                (smooth_bw * np.sqrt(2 * np.pi)))

    while True:
        try:
            batch_llrs = llr_q.get(block=True, timeout=0.01)
        except queue.Empty:
            sleep(0.001)
            continue
        # None stored in queue to indicate exhausted queue
        if batch_llrs is None:
            break
        smooth_vals = np.zeros(smooth_ls.shape[0])
        for llr in batch_llrs:
            smooth_vals += guassian(smooth_ls - llr)
        smooth_llr_q.put((smooth_vals, batch_llrs.shape[0]))


def compute_smooth_density_mp(
        llrs, smooth_bw, smooth_ls, num_proc,
        smooth_mp_batch_size=DEFAULT_SMOOTH_MP_BATCH_SIZE):
    # fill queue with all llrs
    llr_q = mp.Queue()
    for b_start in range(0, llrs.shape[0], smooth_mp_batch_size):
        llr_q.put(llrs[b_start:b_start + smooth_mp_batch_size])
    for _ in range(num_proc):
        llr_q.put(None)

    # start smooth bw processes
    smooth_llr_q = mp.Queue()
    smooth_ps = []
    for _ in range(num_proc):
        p = mp.Process(
            target=_compute_smooth_density_worker,
            args=(llr_q, smooth_llr_q, smooth_bw, smooth_ls),
            daemon=True)
        p.start()
        smooth_ps.append(p)

    bar = tqdm(total=llrs.shape[0], smoothing=0, dynamic_ncols=True)
    # retrieve smooth arrays
    smooth_vals = np.zeros(smooth_ls.shape[0])
    total_nvals = 0
    while any(p.is_alive() for p in smooth_ps):
        try:
            batch_smooth_vals, nvals = smooth_llr_q.get(
                block=True, timeout=0.01)
        except queue.Empty:
            sleep(0.001)
            continue
        smooth_vals += batch_smooth_vals
        total_nvals += nvals
        bar.update(nvals)
    while not smooth_llr_q.empty():
        batch_smooth_vals, nvals = smooth_llr_q.get(block=False)
        smooth_vals += batch_smooth_vals
        total_nvals += nvals
        bar.update(nvals)
    bar.close()

    return smooth_vals / total_nvals


def compute_smooth_mono_density(
        llrs, num_calib_vals, smooth_bw, smooth_ls, num_proc=1):
    def guassian(x):
        return (np.exp(-x ** 2 / (2 * smooth_bw ** 2)) /
                (smooth_bw * np.sqrt(2 * np.pi)))

    if num_proc == 1:
        smooth_vals = np.zeros(num_calib_vals)
        for llr in tqdm(llrs, smoothing=0, dynamic_ncols=True):
            smooth_vals += guassian(smooth_ls - llr)
        smooth_vals /= llrs.shape[0]
    else:
        smooth_vals = compute_smooth_density_mp(
            llrs, smooth_bw, smooth_ls, num_proc)

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
        min_dens_val, diff_eps, llr_buffer, return_plot_info=False,
        num_proc=1):
    smooth_ls = np.linspace(-max_input_llr, max_input_llr,
                            num_calib_vals, endpoint=True)
    LOGGER.info('\tComputing reference emperical density.')
    sm_ref, s_ref = compute_smooth_mono_density(
        ref_llrs, num_calib_vals, smooth_bw, smooth_ls, num_proc)
    LOGGER.info('\tComputing alternative emperical density.')
    sm_alt, s_alt = compute_smooth_mono_density(
        alt_llrs, num_calib_vals, smooth_bw, smooth_ls, num_proc)

    plateau_llr_range = determine_llr_plateau_edge(
        sm_ref, sm_alt, num_calib_vals, diff_eps, llr_buffer, smooth_ls)
    min_dens_llr_range = determine_min_dens_edge(
        sm_ref, sm_alt, min_dens_val, smooth_ls)
    new_input_llr_range = (max(plateau_llr_range[0], min_dens_llr_range[0]),
                           min(plateau_llr_range[1], min_dens_llr_range[1]))
    if new_input_llr_range[1] - new_input_llr_range[0] <= 0:
        raise mh.MegaError('Ground truth smoothed monotonic densities do ' +
                           'not overlap. Consider lowering min_dens_val.')
    if new_input_llr_range[0] != -max_input_llr or \
       new_input_llr_range[1] != max_input_llr:
        LOGGER.info(
            '\tSetting new input llr range for more robust calibration ' +
            '({}, {})'.format(*new_input_llr_range))
        smooth_ls = np.linspace(new_input_llr_range[0], new_input_llr_range[1],
                                num_calib_vals, endpoint=True)
        LOGGER.info('\tComputing new reference emperical density.')
        sm_ref, s_ref = compute_smooth_mono_density(
            ref_llrs, num_calib_vals, smooth_bw, smooth_ls, num_proc)
        LOGGER.info('\tComputing new alternative emperical density.')
        sm_alt, s_alt = compute_smooth_mono_density(
            alt_llrs, num_calib_vals, smooth_bw, smooth_ls, num_proc)

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
        min_dens_val, diff_eps, llr_buffer, return_plot_info=False,
        num_proc=1):
    smooth_ls = np.linspace(-max_input_llr, max_input_llr,
                            num_calib_vals, endpoint=True)

    LOGGER.info('\tComputing reference emperical density.')
    sm_ref, s_ref = compute_smooth_mono_density(
        ref_llrs, num_calib_vals, smooth_bw, smooth_ls, num_proc)

    plateau_llr_range = determine_llr_plateau_edge(
        sm_ref, sm_ref[::-1], num_calib_vals, diff_eps, llr_buffer, smooth_ls)
    min_dens_llr_range = determine_min_dens_edge(
        sm_ref, sm_ref[::-1], min_dens_val, smooth_ls)
    new_input_llr_range = (max(plateau_llr_range[0], min_dens_llr_range[0]),
                           min(plateau_llr_range[1], min_dens_llr_range[1]))
    if new_input_llr_range[0] >= new_input_llr_range[1]:
        raise mh.MegaError('Invalid densities for calibration.')
    if new_input_llr_range[0] != -new_input_llr_range[1]:
        LOGGER.warning(
            'Unexpected new llr range for mirrored calibration: {}'.format(
                str(new_input_llr_range)))
    if new_input_llr_range[0] != -max_input_llr or \
       new_input_llr_range[1] != max_input_llr:
        LOGGER.info(
            '\tSetting new input llr range for more robust calibration ' +
            '({}, {})'.format(*new_input_llr_range))
        smooth_ls = np.linspace(new_input_llr_range[0], new_input_llr_range[1],
                                num_calib_vals, endpoint=True)
        LOGGER.info('\tComputing new reference emperical density.')
        sm_ref, s_ref = compute_smooth_mono_density(
            ref_llrs, num_calib_vals, smooth_bw, smooth_ls, num_proc)

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


#############
# LLR Stats #
#############

def compute_log_probs(alt_llrs):
    """ Compute log probabilities from a set of log likelihood ratios all
    against the reference allele

    Note this function can raise divide by zero and overflow numpy warnings.
    This function should be wrapped at some level by
    with np.errstate(divide='ignore', over='ignore'):
    """
    ref_lp = np.log(1) - np.log1p(np.sum(1 / np.exp(alt_llrs)))
    return ref_lp - alt_llrs


#######################
# Calibration Readers #
#######################

class VarCalibrator(object):
    def _load_calibration(self):
        calib_data = np.load(self.fn)
        self.stratify_type = str(calib_data['stratify_type'])
        assert self.stratify_type == VAR_CALIB_TYPE

        self.num_calib_vals = np.int(calib_data['smooth_nvals'])
        self.max_indel_len = np.int(calib_data['max_indel_len'])

        (self.snp_input_values, self.snp_calib_tables,
         self.del_input_values, self.del_calib_tables,
         self.ins_input_values, self.ins_calib_tables) = (
             {} for _ in range(6))
        # load generic snp
        ref_base, alt_base = GENERIC_BASE, GENERIC_BASE
        snp_type_llr_range = calib_data[
            SNP_LLR_RNG_TMPLT.format(ref_base, alt_base)].copy()
        self.snp_input_values[(ref_base, alt_base)] = np.linspace(
            snp_type_llr_range[0], snp_type_llr_range[1],
            self.num_calib_vals, endpoint=True)
        self.snp_calib_tables[(ref_base, alt_base)] = calib_data[
            SNP_CALIB_TMPLT.format(ref_base, alt_base)].copy()
        # load other base combinations
        for ref_base in mh.ALPHABET:
            for alt_base in set(mh.ALPHABET).difference(ref_base):
                snp_type_llr_range = calib_data[
                    SNP_LLR_RNG_TMPLT.format(ref_base, alt_base)].copy()
                self.snp_input_values[(ref_base, alt_base)] = np.linspace(
                    snp_type_llr_range[0], snp_type_llr_range[1],
                    self.num_calib_vals, endpoint=True)
                self.snp_calib_tables[(ref_base, alt_base)] = calib_data[
                    SNP_CALIB_TMPLT.format(ref_base, alt_base)].copy()
        for indel_len in range(1, self.max_indel_len + 1):
            del_type_llr_range = calib_data[
                DEL_LLR_RNG_TMPLT.format(indel_len)].copy()
            self.del_input_values[indel_len] = np.linspace(
                del_type_llr_range[0], del_type_llr_range[1],
                self.num_calib_vals, endpoint=True)
            self.del_calib_tables[indel_len] = calib_data[
                DEL_CALIB_TMPLT.format(indel_len)].copy()
            ins_type_llr_range = calib_data[
                INS_LLR_RNG_TMPLT.format(indel_len)].copy()
            self.ins_input_values[indel_len] = np.linspace(
                ins_type_llr_range[0], ins_type_llr_range[1],
                self.num_calib_vals, endpoint=True)
            self.ins_calib_tables[indel_len] = calib_data[
                INS_CALIB_TMPLT.format(indel_len)].copy()

    def __init__(self, vars_calib_fn):
        self.fn = vars_calib_fn
        if self.fn is not None:
            self._load_calibration()
        self.calib_loaded = self.fn is not None

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
        seq_len_diff = len(read_ref_seq) - len(read_alt_seq)
        if seq_len_diff == 0:
            ref_seq, alt_seq = simplify_var_seq(read_ref_seq, read_alt_seq)
            # default to a "generic" SNP type (total of all SNP types)
            try:
                calib_table = self.snp_calib_tables[(ref_seq, alt_seq)]
                input_vals = self.snp_input_values[(ref_seq, alt_seq)]
            except KeyError:
                calib_table = self.snp_calib_tables[
                    (GENERIC_BASE, GENERIC_BASE)]
                input_vals = self.snp_input_values[
                    (GENERIC_BASE, GENERIC_BASE)]
        elif seq_len_diff > 0:
            del_len = min(seq_len_diff, self.max_indel_len)
            calib_table = self.del_calib_tables[del_len]
            input_vals = self.del_input_values[del_len]
        else:
            ins_len = min(-seq_len_diff, self.max_indel_len)
            calib_table = self.ins_calib_tables[ins_len]
            input_vals = self.ins_input_values[ins_len]

        idx = np.searchsorted(input_vals, llr, side='left')
        # full closest search would be:
        # if idx > 0 and (idx == self.num_calib_vals or
        #                np.abs(llr - input_vals[idx - 1]) <
        #                np.abs(llr - input_vals[idx])):
        # but for performance just adjust last index
        if idx == self.num_calib_vals:
            idx -= 1
        return calib_table[idx]


class ModCalibrator(object):
    def _load_calibration(self):
        calib_data = np.load(self.fn)
        self.stratify_type = str(calib_data[MOD_STRAT_TYPE_TXT])
        assert self.stratify_type in MOD_VALID_STRAT_TYPES

        self.num_calib_vals = np.int(calib_data[SMOOTH_NVALS_TXT])
        self.mod_bases = calib_data[MOD_BASES_TXT]
        self.mod_base_calibs = {}
        for mod_base in self.mod_bases:
            mod_llr_range = calib_data[mod_base + LLR_RANGE_SUFFIX].copy()
            input_vals = np.linspace(
                mod_llr_range[0], mod_llr_range[1],
                self.num_calib_vals, endpoint=True)
            mod_calib_table = calib_data[
                mod_base + CALIB_TABLE_SUFFIX].copy()
            self.mod_base_calibs[mod_base] = (input_vals, mod_calib_table)

    def __init__(self, mods_calib_fn):
        self.fn = mods_calib_fn
        if self.fn is not None:
            self._load_calibration()
        self.calib_loaded = self.fn is not None

    def calibrate_llr(self, llr, mod_base):
        if not self.calib_loaded or mod_base not in self.mod_base_calibs:
            return llr

        input_vals, calib_table = self.mod_base_calibs[mod_base]
        idx = np.searchsorted(input_vals, llr, side='left')
        # full closest search would be:
        # if idx > 0 and (idx == self.num_calib_vals or
        #                np.abs(llr - input_vals[idx - 1]) <
        #                np.abs(llr - input_vals[idx])):
        # but for performance just adjust last index
        if idx == self.num_calib_vals:
            idx -= 1
        return calib_table[idx]


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
