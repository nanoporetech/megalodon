import re
import sys
import queue
import sqlite3
import datetime
from time import sleep
from collections import defaultdict, namedtuple, OrderedDict

import numpy as np

from megalodon import (
    calibration, decode, logging, mapping, megalodon_helper as mh)
from megalodon._version import MEGALODON_VERSION


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
SELECT DISTINCT chrm, strand, pos FROM mods)"""
SEL_UNIQ_MODS = 'SELECT DISTINCT chrm, strand, pos FROM mods'
SEL_MOD_STATS = '''
SELECT * FROM mods WHERE chrm=? AND strand=? AND pos=?'''

BIN_THRESH_NAME = 'binary_threshold'
EM_NAME = 'em'
AGG_METHOD_NAMES = set((BIN_THRESH_NAME, EM_NAME))
AGG_INFO = namedtuple('AGG_INFO', ('method', 'binary_threshold'))
DEFAULT_AGG_INFO = AGG_INFO(BIN_THRESH_NAME, [-1, 1])

FIXED_VCF_MI = [
    'INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
    'INFO=<ID=SN,Number=1,Type=String,Description="Strand">',
]
MOD_MI_TMPLTS = [
    'FORMAT=<ID={0}DP,Number=1,Type=Integer,Description=' +
    '"Valid Read Depth for {1}">',
    'FORMAT=<ID={0},Number=1,Type=Float,Description='+
    '"{1} Modified Base Proportion">']


################################
##### Per-read Mod Scoring #####
################################

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
        r_ref_pos, edge_buffer, r_ref_seq, np_ref_seq, rl_cumsum,
        r_to_q_poss, r_post, post_mapped_start, mods_info):
    def iter_motif_sites(r_ref_seq):
        max_pos = len(r_ref_seq) - edge_buffer
        for motif, rel_pos, mod_base, raw_motif in mods_info.all_mod_motifs:
            for m_pos in [
                    m.start() + rel_pos for m in motif.finditer(r_ref_seq)]:
                if m_pos < edge_buffer: continue
                if m_pos > max_pos: break
                yield m_pos, mod_base, raw_motif
        return


    # call all mods overlapping this read
    r_mod_scores = []
    for pos, mod_base, raw_motif in iter_motif_sites(r_ref_seq):
        pos_bb, pos_ab = min(mods_info.mod_context_bases, pos), min(
            mods_info.mod_context_bases, np_ref_seq.shape[0] - pos - 1)
        pos_ref_seq = np_ref_seq[pos - pos_bb:pos + pos_ab + 1]
        pos_ref_mods = np.zeros_like(pos_ref_seq)
        pos_alt_mods = pos_ref_mods.copy()
        pos_alt_mods[pos_bb] = mods_info.str_to_int_mod_labels[mod_base]

        blk_start, blk_end = (rl_cumsum[r_to_q_poss[pos - pos_bb]],
                              rl_cumsum[r_to_q_poss[pos + pos_ab]])
        if blk_end - blk_start < (mods_info.mod_context_bases * 2) + 1:
            # no valid mapping over large inserted query bases
            # i.e. need as many "events/strides" as bases for valid mapping
            continue

        loc_ref_score = score_mod_seq(
            r_post, pos_ref_seq, pos_ref_mods, mods_info.can_mods_offsets,
            post_mapped_start + blk_start, post_mapped_start + blk_end,
            mods_info.mod_all_paths)
        loc_alt_score = score_mod_seq(
            r_post, pos_ref_seq, pos_alt_mods, mods_info.can_mods_offsets,
            post_mapped_start + blk_start, post_mapped_start + blk_end,
            mods_info.mod_all_paths)
        if loc_ref_score is None or loc_alt_score is None:
            raise mh.MegaError('Score computation error (memory error)')

        m_ref_pos = (pos + r_ref_pos.start if r_ref_pos.strand == 1 else
                     r_ref_pos.end - pos - 1)
        # calibrate llr scores
        calib_llr = mods_info.calibrate_llr(
            loc_ref_score - loc_alt_score, mod_base)
        r_mod_scores.append((m_ref_pos, calib_llr, raw_motif, mod_base))

    return r_mod_scores


###############################
##### Per-read Mod Output #####
###############################

def annotate_mods(r_start, ref_seq, r_mod_scores, strand):
    """ Annotate reference sequence with called modified bases.

    Note: Reference sequence is in read orientation and mod calls are in
    genome coordiates.
    """
    mod_seqs = []
    prev_pos = 0
    if strand == -1:
        ref_seq = ref_seq[::-1]
    for mod_pos, score, _, mod_base in sorted(r_mod_scores):
        # called canonical
        # TODO: handle models with more than one mod per canonical base
        if score >= 0: continue
        mod_seqs.append(ref_seq[prev_pos:mod_pos - r_start] + mod_base)
        prev_pos = mod_pos - r_start + 1
    mod_seqs.append(ref_seq[prev_pos:])
    mod_seq = ''.join(mod_seqs)
    if strand == -1:
        mod_seq = mod_seq[::-1]

    return mod_seq

def _get_mods_queue(
        mods_q, mods_conn, mods_db_fn, mods_txt_fn, db_safety,
        pr_refs_fn, pr_ref_filts):
    def get_mod_call():
        # note strand is +1 for fwd or -1 for rev
        r_mod_scores, (
            read_id, chrm, strand, r_start, ref_seq, read_len, q_st, q_en,
            cigar) = mods_q.get(block=False)
        mods_db.executemany(ADDMANY_MODS, [
            (read_id, chrm, strand, pos, score, mod_base, raw_motif)
            for pos, score, raw_motif, mod_base in r_mod_scores])
        if mods_txt_fp is not None:
            # would involve batching and creating several conversion tables
            # for var strings (read_if and chrms).
            mods_txt_fp.write('\n'.join((
                '{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                    read_id, chrm, strand, pos, score, raw_motif, mod_base)
                for pos, score, raw_motif, mod_base in r_mod_scores)) + '\n')
            mods_txt_fp.flush()
        if pr_refs_fn is not None:
            if not mapping.read_passes_filters(
                    pr_ref_filts, read_len, q_st, q_en, cigar):
                return

            pr_refs_fp.write('>{}\n{}\n'.format(read_id, annotate_mods(
                r_start, ref_seq, r_mod_scores, strand)))
            pr_refs_fp.flush()

        return


    mods_db = sqlite3.connect(mods_db_fn)
    if db_safety < 2:
        mods_db.execute(SET_ASYNC_MODE)
    if db_safety < 1:
        mods_db.execute(SET_NO_ROLLBACK_MODE)
    mods_db.execute(CREATE_MODS_TBLS)
    mods_txt_fp = None if mods_txt_fn is None else open(mods_txt_fn, 'w')

    if pr_refs_fn is not None:
        pr_refs_fp = open(pr_refs_fn, 'w')

    while True:
        try:
            get_mod_call()
        except queue.Empty:
            if mods_conn.poll():
                break
            sleep(0.1)
            continue

    while not mods_q.empty():
        get_mod_call()
    if mods_txt_fp is not None: mods_txt_fp.close()
    if pr_refs_fn is not None: pr_refs_fp.close()
    mods_db.execute(CREATE_MODS_IDX)
    mods_db.commit()
    mods_db.close()

    return


####################
##### Mod Info #####
####################

class ModInfo(object):
    single_letter_code = {
        'A':'A', 'C':'C', 'G':'G', 'T':'T', 'B':'[CGT]',
        'D':'[AGT]', 'H':'[ACT]', 'K':'[GT]', 'M':'[AC]',
        'N':'[ACGT]', 'R':'[AG]', 'S':'[CG]', 'V':'[ACG]',
        'W':'[AT]', 'Y':'[CT]'}

    def _parse_mod_motifs(self, all_mod_motifs_raw):
        # note only works for mod_refactor models currently
        self.all_mod_motifs = []
        if all_mod_motifs_raw is None or len(all_mod_motifs_raw) == 0:
            for can_base, mod_bases in self.can_base_mods.items():
                for mod_base in mod_bases:
                    self.all_mod_motifs.append((
                        re.compile(can_base), 0, mod_base, can_base))
        else:
            # parse detection motifs
            for mod_base, raw_motif, pos in all_mod_motifs_raw:
                assert len(mod_base) == 1, (
                    'Modfied base must be a single character. Got {}'.format(
                        mod_base))
                assert mod_base in self.str_to_int_mod_labels, (
                    'Modified base label ({}) not found in model ' +
                    'alphabet ({}).').format(
                        mod_base, list(self.str_to_int_mod_labels.keys()))
                pos = int(pos)
                can_base = next(
                    can_base for can_base, can_mods in
                    self.can_base_mods.items() if mod_base in can_mods)
                assert (can_base == raw_motif[pos]), (
                    'Invalid modified base motif. Raw motif modified ' +
                    'position ({}) base ({}) does not match ' +
                    'collapsed alphabet value ({}).').format(
                        pos, raw_motif[pos], can_base)
                motif = re.compile(''.join(
                    self.single_letter_code[letter] for letter in raw_motif))
                self.all_mod_motifs.append((motif, pos, mod_base, raw_motif))

        return

    def __init__(
            self, model_info, all_mod_motifs_raw=None, mod_all_paths=False,
            write_mods_txt=None, mod_context_bases=None,
            do_output_mods=False, do_pr_ref_mods=False, mods_calib_fn=None):
        logger = logging.get_logger()
        # this is pretty hacky, but these attributes are stored here as
        # they are generally needed alongside other alphabet info
        # don't want to pass all of these parameters around individually though
        # as this would make function signatures too complicated
        self.mod_all_paths = mod_all_paths
        self.write_mods_txt = write_mods_txt
        self.mod_context_bases = mod_context_bases
        self.do_output_mods = do_output_mods
        self.do_pr_ref_mods = do_pr_ref_mods
        self.mod_long_names = model_info.mod_long_names
        self.calib_table = calibration.ModCalibrator(mods_calib_fn)

        self.alphabet = model_info.can_alphabet
        self.ncan_base = len(self.alphabet)
        try:
            self.alphabet = self.alphabet.decode()
        except:
            pass
        if model_info.is_cat_mod:
            # TODO also output "(alt to C)" for each mod
            logger.info(
                'Using canoncical alphabet {} and modified bases {}.'.format(
                    self.alphabet, ' '.join(
                        '{}={}'.format(*mod_b)
                        for mod_b in model_info.mod_long_names)))
        else:
            logger.info(
                'Using canoncical alphabet {}.'.format(self.alphabet))

        self.nbase = len(self.alphabet)
        self.n_can_state = (self.ncan_base + self.ncan_base) * (
            self.ncan_base + 1)
        if model_info.is_cat_mod:
            self.nmod_base = model_info.n_mods
            self.can_base_mods = model_info.can_base_mods
            self.can_mods_offsets = model_info.can_indices
            self.str_to_int_mod_labels = model_info.str_to_int_mod_labels
            assert (
                model_info.output_size - self.n_can_state ==
                self.nmod_base + 1), (
                    'Alphabet ({}) and model number of modified bases ({}) ' +
                    'do not agree.').format(
                        self.alphabet,
                        model_info.output_size - self.n_can_state - 1)
        else:
            self.nmod_base = 0
            self.can_base_mods = {}
            self.can_mods_offsets = None
            self.str_to_int_mod_labels = None

        # parse mod motifs or use "swap" base if no motif provided
        self._parse_mod_motifs(all_mod_motifs_raw)

        return

    def calibrate_llr(self, score, mod_base):
        return self.calib_table.calibrate_llr(score, mod_base)


#########################
##### modVCF Writer #####
#########################

class ModSite(object):
    """ Modified base site for entry into ModVcfWriter.
    Currently only handles a single sample.
    """
    def __init__(
            self, chrom, pos, strand, ref, mods, id='.', qual='.', filter='.',
            info=None, sample_dict=None):
        self.chrom = chrom
        self.pos = int(pos)
        self.ref = ref.upper()
        self.alt = mods
        self.id = str(id)
        self.qual = qual
        self.filter = str(filter)
        self.strand = strand
        if info is None:
            info = {}
        if 'STRD' not in info:
            info['STRD'] = strand
        self.info_dict = info
        if sample_dict is None:
            sample_dict = OrderedDict()
        self.sample_dict = sample_dict
        return

    @property
    def _sorted_format_keys(self):
        sorted_keys = sorted(self.sample_dict.keys())
        return sorted_keys
    @property
    def format(self):
        return ':'.join(map(str, self._sorted_format_keys))
    @property
    def sample(self):
        return ':'.join((str(self.sample_dict[k])
                         for k in self._sorted_format_keys))
    @property
    def info(self):
        str_tags = []
        for key, value in self.info_dict.items():
            # If key is of type 'Flag', print only key, else 'key=value'
            if value is True:
                str_tags.append(key)
            else:
                if isinstance(value, (tuple, list)):
                    value = ','.join(map(str, value))
                str_tags.append('{}={}'.format(key, value))
        return ':'.join(str_tags)

    def add_tag(self, tag, value=None):
        self.info_dict[tag] = value

    def add_sample_field(self, tag, value=None):
        self.sample_dict[tag] = value

    def add_mod_props(self, mod_props):
        with np.errstate(divide='ignore'):
            can_pl = -10 * np.log10(1 - sum((x[1] for x in mod_props)))
        self.qual = '{:.0f}'.format(
            np.abs(np.around(np.minimum(can_pl, mh.MAX_PL_VALUE))))
        for mod_name, mod_prop, _ in mod_props:
            self.add_sample_field(mod_name, '{:.4f}'.format(mod_prop))
        return

    def __eq__(self, mod2):
        return (self.chrm, self.pos, self.strand) == (
            mod2.chrm, mod2.pos, mod2.strand)
    def __ne__(self, var2):
        return (self.chrm, self.pos, self.strand) != (
            mod2.chrm, mod2.pos, mod2.strand)
    def __lt__(self, var2):
        return (self.chrm, self.pos, self.strand) < (
            mod2.chrm, mod2.pos, mod2.strand)
    def __le__(self, var2):
        return (self.chrm, self.pos, self.strand) <= (
            mod2.chrm, mod2.pos, mod2.strand)
    def __gt__(self, var2):
        return (self.chrm, self.pos, self.strand) > (
            mod2.chrm, mod2.pos, mod2.strand)
    def __ge__(self, var2):
        return (self.chrm, self.pos, self.strand) >= (
            mod2.chrm, mod2.pos, mod2.strand)


class ModVcfWriter(object):
    """ modVCF writer class
    """
    version_options = set(['4.2',])
    def __init__(
            self, filename, mods, mode='w',
            header=('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
                    'INFO', 'FORMAT', 'SAMPLE'),
            extra_meta_info=FIXED_VCF_MI, version='4.2', ref_fn=None,
            ref_names_and_lens=None):
        self.filename = filename
        self.mods = mods
        self.mode = mode
        self.header = header
        if version not in self.version_options:
            raise ValueError('version must be one of {}'.format(
                self.version_options))
        self.version = version
        contig_mis = [] if ref_names_and_lens is None else [
            mh.CONTIG_MI.format(ref_name, ref_len)
            for ref_name, ref_len in zip(*ref_names_and_lens)]
        self.meta = [
            mh.VCF_VERSION_MI.format(self.version),
            mh.FILE_DATE_MI.format(
                datetime.date.today().strftime("%Y%m%d")),
            mh.SOURCE_MI.format(MEGALODON_VERSION),
            mh.REF_MI.format(ref_fn)] + contig_mis + extra_meta_info + [
                mod_tmplt.format(*mod_name) for mod_name in self.mods
                for mod_tmplt in MOD_MI_TMPLTS]

        self.handle = open(self.filename, self.mode, encoding='utf-8')
        self.handle.write('\n'.join('##' + line for line in self.meta) + '\n')
        self.handle.write('#' + '\t'.join(self.header) + '\n')
        return

    def close(self):
        self.handle.close()
        return

    def write_mod_site(self, mod_site):
        elements = [getattr(mod_site, field.lower()) for field in self.header]
        elements = ['.' if e == '' else e for e in elements]
        # VCF POS field is 1-based
        elements[self.header.index('POS')] += 1
        self.handle.write('{}\n'.format('\t'.join(map(str, elements))))
        self.handle.flush()

        return


##################################
##### Mods Aggregation Class #####
##################################

class AggMods(mh.AbstractAggregationClass):
    """ Class to assist in database queries for per-site aggregation of
    modified base calls over reads.
    """
    def __init__(self, mods_db_fn, agg_info=DEFAULT_AGG_INFO):
        # open as read only database
        self.mods_db = sqlite3.connect(mods_db_fn, uri=True)
        self.n_uniq_mods = None
        assert agg_info.method in AGG_METHOD_NAMES
        self.agg_method = agg_info.method
        self.binary_thresh = agg_info.binary_threshold
        if type(self.binary_thresh) in (float, int):
            self.binary_thresh = [self.binary_thresh, self.binary_thresh]
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
        return [MOD_DATA(*mod_stats) for mod_stats in self.mods_db.execute(
            SEL_MOD_STATS, mod_loc)]

    def est_binary_thresh(self, pos_scores):
        pos_mod = np.less(pos_scores, self.binary_thresh[0])
        mod_cov = pos_mod.sum()
        valid_cov = np.logical_or(
            pos_mod, np.greater(pos_scores, self.binary_thresh[1])).sum()
        if valid_cov == 0:
            return 0, valid_cov
        return mod_cov / float(valid_cov), valid_cov

    def est_em_prop(
            self, pos_scores, max_iters=5, conv_tol=0.005,
            init_thresh=0, min_prop=0.01, max_prop=0.99):
        """ Estimate proportion of modified bases at a position via EM
        computation
        """
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

    def emp_em(self, pos_scores, max_iters):
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

    def compute_mod_stats(self, mod_loc, agg_method=None):
        if agg_method is None:
            agg_method = self.agg_method
        if agg_method not in AGG_METHOD_NAMES:
            raise NotImplementedError(
                'No modified base proportion estimation method: {}'.format(
                    agg_method))

        pr_mod_stats = self.get_per_read_mod_stats(mod_loc)
        mod_type_stats = defaultdict(list)
        for r_stats in pr_mod_stats:
            mod_type_stats[r_stats.mod_base].append(r_stats)
        mt_stats = []
        for mod_base, mt_reads in mod_type_stats.items():
            mt_llrs = np.array([r_stats.score for r_stats in mt_reads])
            if agg_method == BIN_THRESH_NAME:
                prop_est, valid_cov = self.est_binary_thresh(mt_llrs)
            else:
                prop_est, valid_cov = self.est_em_prop(mt_llrs)
            mt_stats.append((mod_base, prop_est, valid_cov))
        r0_stats = pr_mod_stats[0]
        strand = '+' if r0_stats.strand == 1 else '-'
        site_motifs = ','.join(sorted(set(rs.motif for rs in pr_mod_stats)))
        mod_site = ModSite(
            chrom=r0_stats.chrm, pos=r0_stats.pos, strand=strand,
            ref=site_motifs, mods=','.join(mod_type_stats.keys()),
            id='_'.join(map(str, (r0_stats.chrm, r0_stats.pos, strand))))
        mod_site.add_mod_props(mt_stats)
        mod_site.add_tag('DP', '{}'.format(len(pr_mod_stats)))
        for mod_base, _, valid_cov in mt_stats:
            mod_site.add_sample_field(
                '{}DP'.format(mod_base), '{}'.format(int(valid_cov)))
        return mod_site

    def close(self):
        self.mods_db.close()
        return


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
