import queue
import sqlite3
from time import sleep
from collections import namedtuple

import numpy as np

from megalodon import decode, megalodon_helper as mh


FIELD_NAMES = ('read_id', 'chrm', 'strand', 'pos', 'score',
               'mod_base', 'motif')
MOD_DATA = namedtuple('MOD_DATA', FIELD_NAMES)
CREATE_MODS_TBLS = """
CREATE TABLE mods (
    {} TEXT,
    {} TEXT,
    {} INTEGER,
    {} INTEGER,
    {} FLOAT,
    {} TEXT,
    {} TEXT
)""".format(*FIELD_NAMES)

SET_NO_ROLLBACK_MODE='PRAGMA journal_mode = OFF'
SET_ASYNC_MODE='PRAGMA synchronous = OFF'

ADDMANY_MODS = "INSERT INTO mods VALUES (?,?,?,?,?,?,?)"
CREATE_MODS_IDX = "CREATE INDEX mod_pos ON mods (chrm, strand, pos)"

COUNT_UNIQ_MODS = """
SELECT COUNT(*) FROM (
SELECT DISTINCT chrm, strand, pos, mod_base FROM mods)"""
SEL_UNIQ_MODS = 'SELECT DISTINCT chrm, strand, pos, mod_base FROM mods'
SEL_MOD_STATS = '''
SELECT * FROM mods WHERE chrm=? AND strand=? AND pos=? and mod_base=?'''

BIN_THRESH_NAME = 'binary_threshold'
EM_NAME = 'em'
PROP_METHOD_NAMES = set((BIN_THRESH_NAME, EM_NAME))


def score_mod_seq(
        tpost, seq, mod_cats, can_mods_offsets,
        tpost_start=0, tpost_end=None, all_paths=False):
    """Score a section of log transition posteriors against a proposed sequence
    using a global mapping.
    :param tpost: `ndarray` containing log transition posteriors to be scored
    :param seq: `ndarray` containing integers encoding proposed sequence
    :param mod_cats: `ndarray` containing integers encoding proposed modified base labels
    :param can_mods_offsets: `ndarray` containing integers encoding proposed modified base labels
    :param tpost_start: start position within post (Default: 0)
    :param tpost_end: end position within post (Default: full posterior)
    :param all_paths: boolean to produce the forwards all paths score (default Viterbi best path)
    """
    seq = seq.astype(np.uintp)
    if tpost_end is None:
        tpost_end = tpost.shape[0]

    return decode.score_mod_seq(
        tpost, seq, mod_cats, can_mods_offsets, tpost_start, tpost_end,
        all_paths)

def call_read_mods(
        r_ref_pos, edge_buffer, context_bases, r_ref_seq, np_ref_seq, rl_cumsum,
        r_to_q_poss, r_post, post_mapped_start, alphabet_info):
    def iter_motif_sites(r_ref_seq):
        max_pos = len(r_ref_seq) - edge_buffer
        for motif, rel_pos, mod_base, raw_motif in alphabet_info.all_mod_motifs:
            for m_pos in [
                    m.start() + rel_pos for m in motif.finditer(r_ref_seq)]:
                if m_pos < edge_buffer: continue
                if m_pos > max_pos: break
                yield m_pos, mod_base, raw_motif
        return


    # call all mods overlapping this read
    r_mod_calls = []
    for pos, mod_base, raw_motif in iter_motif_sites(r_ref_seq):
        pos_bb, pos_ab = min(context_bases, pos), min(
            context_bases, np_ref_seq.shape[0] - pos - 1)
        pos_ref_seq = np_ref_seq[pos - pos_bb:pos + pos_ab + 1]
        pos_ref_mods = np.zeros_like(pos_ref_seq)
        pos_alt_mods = pos_ref_mods.copy()
        pos_alt_mods[pos_bb] = alphabet_info.str_to_int_mod_labels[mod_base]

        blk_start, blk_end = (rl_cumsum[r_to_q_poss[pos - pos_bb]],
                              rl_cumsum[r_to_q_poss[pos + pos_ab]])
        if blk_end - blk_start < (context_bases * 2) + 1:
            # no valid mapping over large inserted query bases
            # i.e. need as many "events/strides" as bases for valid mapping
            continue

        loc_ref_score = score_mod_seq(
            r_post, pos_ref_seq, pos_ref_mods, alphabet_info.can_mods_offsets,
            post_mapped_start + blk_start, post_mapped_start + blk_end,
            alphabet_info.mod_all_paths)
        loc_alt_score = score_mod_seq(
            r_post, pos_ref_seq, pos_alt_mods, alphabet_info.can_mods_offsets,
            post_mapped_start + blk_start, post_mapped_start + blk_end,
            alphabet_info.mod_all_paths)
        if loc_ref_score is None or loc_alt_score is None:
            raise mh.MegaError('Score computation error (memory error)')

        m_ref_pos = (pos + r_ref_pos.start if r_ref_pos.strand == 1 else
                     r_ref_pos.end - pos - 1)
        r_mod_calls.append((
            m_ref_pos, loc_ref_score - loc_alt_score, raw_motif, mod_base))

    return r_mod_calls

def _get_mods_queue(mods_q, mods_conn, mods_db_fn, mods_txt_fn, db_safety):
    mods_db = sqlite3.connect(mods_db_fn)
    mods_db_c = mods_db.cursor()
    if db_safety < 2:
        mods_db_c.execute(SET_ASYNC_MODE)
    if db_safety < 1:
        mods_db_c.execute(SET_NO_ROLLBACK_MODE)
    mods_db_c.execute(CREATE_MODS_TBLS)
    mods_txt_fp = None if mods_txt_fn is None else open(mods_txt_fn, 'w')

    while True:
        try:
            # note strand is +1 for fwd or -1 for rev
            r_mod_calls, (read_id, chrm, strand) = mods_q.get(block=False)
            mods_db_c.executemany(ADDMANY_MODS, [
                (read_id, chrm, strand, pos, score, mod_base, raw_motif)
                for pos, score, raw_motif, mod_base in r_mod_calls])
            if mods_txt_fp is not None:
                # would involve batching and creating several conversion tables
                # for var strings (read_if and chrms).
                mods_txt_fp.write('\n'.join((
                    '{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                        read_id, chrm, strand, pos, score, raw_motif, mod_base)
                    for pos, score, raw_motif, mod_base in r_mod_calls)) + '\n')
                mods_txt_fp.flush()
        except queue.Empty:
            if mods_conn.poll():
                break
            sleep(0.1)
            continue

    while not mods_q.empty():
        r_mod_calls, (read_id, chrm, strand) = mods_q.get(block=False)
        mods_db_c.execute(ADDMANY_MODS, [
            (read_id, chrm, strand, pos, score, mod_base, raw_motif)
            for pos, score, raw_motif, mod_base in r_mod_calls])
        if mods_txt_fp is not None:
            mods_txt_fp.write('\n'.join((
                '{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                    read_id, chrm, strand, pos, score, raw_motif, mod_base)
                for pos, score, raw_motif, mod_base in r_mod_calls)) + '\n')
            mods_txt_fp.flush()
    if mods_txt_fp is not None: mods_txt_fp.close()
    mods_db_c.execute(CREATE_MODS_IDX)
    mods_db.commit()
    mods_db.close()

    return




##################################
##### Mods Aggregation Class #####
##################################

class AggMods(mh.AbstractAggregationClass):
    """ Class to assist in database queries for per-site aggregation of
    modified base calls over reads.
    """
    def _load_calibration(self, mods_calib_fn):
        mod_calib_data = np.load(mods_calib_fn)
        self.stratify_type = str(mod_calib_data['stratify_type'])
        self.max_input_llhr = np.int(mod_calib_data['smooth_max'])
        self.num_calib_vals = np.int(mod_calib_data['smooth_nvals'])
        self.discrete_step = 2 * self.max_input_llhr / (self.num_calib_vals - 1)
        return mod_calib_data['global_calibration_table'].copy()

    def __init__(
            self, mods_db_fn, mods_calib_fn=None,
            prop_method=BIN_THRESH_NAME, binary_thresh=0):
        # open as read only database
        self.mods_db = sqlite3.connect(mods_db_fn, uri=True)
        self.n_uniq_mods = None
        self.calib_table = (None if mods_calib_fn is None else
                            self._load_calibration(mods_calib_fn))
        assert prop_method in PROP_METHOD_NAMES
        self.prop_method = prop_method
        self.binary_thresh = binary_thresh
        if type(self.binary_thresh) in (float, int):
            self.binary_thresh = [binary_thresh, binary_thresh]
        return

    def num_uniq(self):
        if self.n_uniq_mods is None:
            self.n_uniq_mods = self.mods_db.execute(
                COUNT_UNIQ_MODS).fetchone()[0]
        return self.n_uniq_mods

    def iter_uniq(self):
        for q_val in self.mods_db.execute(SEL_UNIQ_MODS):
            yield q_val
        return

    def get_per_read_mod_stats(self, mod_loc):
        return [MOD_DATA(*snp_stats) for mod_stats in self.mods_db.execute(
            SEL_MOD_STATS, mod_loc)]

    def calibrate_llhrs(self, llhrs):
        return self.calib_table[np.around((
            np.clip(llhrs, -self.max_input_llhr, self.max_input_llhr) +
            self.max_input_llhr) / self.discrete_step).astype(int)]

    def est_binary_thresh(pos_scores):
        pos_mod = np.less(pos_scores, self.binary_thresh[0])
        mod_cov = pos_mod.sum()
        valid_cov = np.logical_or(
            pos_mod, np.greater(pos_scores, binary_thresh[1])).sum()
        if valid_cov == 0:
            return 0, valid_cov
        return mod_cov / float(valid_cov), valid_cov

    def est_em_prop(
            pos_scores, max_iters=5, conv_tol=0.005,
            init_thresh=0, min_prop=0.01, max_prop=0.99):
        curr_mix_prop = np.clip(np.mean(pos_scores < init_thresh),
                                min_prop, max_prop)
        for _ in range(max_iters):
            prev_mix_prop = curr_mix_prop
            if prev_mix_prop < min_prop:
                return 0.0, pos_scores.shape[0]
            if prev_mix_prop > max_prop:
                return 1.0, pos_scores.shape[0]
            curr_mix_prop = np.mean(1 / (1 + (
                np.exp(pos_scores) * (1 - prev_mix_prop) / prev_mix_prop)))
            if np.abs(curr_mix_prop - prev_mix_prop) < conv_tol:
                break

        return curr_mix_prop, pos_scores.shape[0]

    def emp_em(pos_scores, max_iters):
        """ Emperical EM estimation

        This method does not actually work yet. Need to compute
        emperical proportion array
        """
        # pre-computed array from ground truth
        # rows are mixture rate, cols are per-read scores
        mod_prob_tbl = np.array()

        curr_mix_prop = np.mean(pos_scores < 0)
        for _ in range(max_iters):
            curr_mix_prop = np.mean(mod_prob_tbl[curr_mix_prop, pos_scores])

        return curr_mix_prop, pos_scores.shape[0]

    def compute_mod_stats(self, mod_loc, prop_method=None):
        pr_mod_stats = self.get_per_read_mod_stats(mod_loc)
        loc_llhrs = np.array([r_stats.score for r_stats in pr_snp_stats])
        if self.calib_table is not None:
            llhrs = self.calibrate_llhrs(llhrs)
        if prop_method is None:
            prop_method = self.prop_method
        if prop_method == BIN_THRESH_NAME:
            prop_est = self.est_binary_thresh(llhrs)
        elif prop_method == EM_NAME:
            prop_est = self.est_em_prop(llhrs)
        else:
            raise NotImplementedError(
                'No modified base proportion estimation method: {}'.format(
                    prop_method))
        return prop_est

    def close(self):
        self.mods_db.close()
        return


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
