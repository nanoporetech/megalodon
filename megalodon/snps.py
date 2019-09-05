import os
import sys
import queue
import sqlite3
import datetime
from time import sleep
from array import array
import multiprocessing as mp
from operator import itemgetter
from itertools import product, combinations, groupby
from collections import defaultdict, namedtuple, OrderedDict

import pysam
import numpy as np
from scipy import stats

from megalodon import (calibration, decode, logging, mapping,
                       megalodon_helper as mh)
from megalodon._version import MEGALODON_VERSION


_DEBUG_PER_READ = False
_RAISE_VARIANT_PROCESSING_ERRORS = False

VARIANT_DATA = namedtuple('VARIANT_DATA', (
    'np_ref', 'np_alts', 'id', 'chrom', 'start', 'stop',
    'ref', 'alts', 'ref_start', ))
# set default value of None for ref, alts and ref_start
VARIANT_DATA.__new__.__defaults__ = (None, None, None)

DIPLOID_MODE = 'diploid'
HAPLIOD_MODE = 'haploid'

FIELD_NAMES = ('read_id', 'chrm', 'strand', 'pos', 'score',
               'ref_seq', 'alt_seq', 'snp_id', 'test_start', 'test_end')
SNP_DATA = namedtuple('SNP_DATA', FIELD_NAMES)
CREATE_SNPS_TBLS = """
CREATE TABLE snps (
    {} TEXT,
    {} TEXT,
    {} INTEGER,
    {} INTEGER,
    {} FLOAT,
    {} TEXT,
    {} TEXT,
    {} TEXT,
    {} INTEGER,
    {} INTEGER
)""".format(*FIELD_NAMES)

SET_NO_ROLLBACK_MODE='PRAGMA journal_mode = OFF'
SET_ASYNC_MODE='PRAGMA synchronous = OFF'

ADDMANY_SNPS = "INSERT INTO snps VALUES (?,?,?,?,?,?,?,?,?,?)"
CREATE_SNPS_IDX = '''
CREATE INDEX snp_pos ON snps (chrm, test_start, test_end)'''

COUNT_UNIQ_SNPS = """
SELECT COUNT(*) FROM (
SELECT DISTINCT chrm, test_start, test_end FROM snps)"""
SEL_UNIQ_SNP_ID = '''
SELECT DISTINCT chrm, test_start, test_end FROM snps'''
SEL_SNP_STATS = '''
SELECT * FROM snps WHERE chrm IS ? AND test_start IS ? AND test_end IS ?'''

SAMPLE_NAME = 'SAMPLE'
# specified by sam format spec
WHATSHAP_MAX_QUAL = 40
WHATSHAP_RG_ID = '1'
FIXED_VCF_MI = [
    'phasing=none',
    'INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
    'FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    'FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
    'FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
    ('FORMAT=<ID=GL,Number=G,Type=Float,' +
     'Description="Log10 likelihoods for genotypes">'),
    ('FORMAT=<ID=PL,Number=G,Type=Integer,' +
     'Description="Normalized, Phred-scaled likelihoods for genotypes">')
]
FORMAT_LOG_PROB_MI = (
    'FORMAT=<ID=LOG_PROBS,Number=A,Type=String,' +
    'Description="Per-read log10 likelihoods for alternative ' +
    'alleles (semi-colon separated)">')


############################
##### Helper Functions #####
############################

def logsumexp(x):
    x_max = x.max()
    return np.log(np.sum(np.exp(x - x_max))) + x_max


################################
##### Per-read SNP Scoring #####
################################

def write_per_read_debug(
        snp_ref_pos, snp_id, read_ref_pos, np_s_snp_ref_seq, np_s_snp_alt_seqs,
        np_s_context_seqs, loc_contexts_ref_lps, loc_contexts_alts_lps,
        w_context, logger):
    ref_seq = mh.int_to_seq(np_s_snp_ref_seq)
    if read_ref_pos.strand == -1:
        ref_seq = mh.revcomp(ref_seq)
    alts_seq = [mh.int_to_seq(np_alt) for np_alt in np_s_snp_alt_seqs]
    if read_ref_pos.strand == -1:
        alts_seq = [mh.revcomp(alt_seq) for alt_seq in alts_seq]
    ','.join(alts_seq)

    context_seqs = []
    for up_context_seq, dn_context_seq in np_s_context_seqs:
        up_seq, dn_seq = (mh.int_to_seq(up_context_seq),
                          mh.int_to_seq(dn_context_seq))
        if read_ref_pos.strand == -1:
            up_seq, dn_seq = mh.revcomp(dn_seq), mh.revcomp(up_seq)
        context_seqs.append((up_seq, dn_seq))
    out_txt = '\n'
    for ref_lp, alt_lps, (up_seq, dn_seq) in zip(
            loc_contexts_ref_lps, zip(*loc_contexts_alts_lps), context_seqs):
        out_txt += ('VARIANT_FULL_DATA: {}\t{}\t{}\t{}\t{}[{}]{}\t{}\t' +
                    '{:.2f}\t{}\t{}\n').format(
                        read_ref_pos.chrm, read_ref_pos.strand, snp_ref_pos,
                        snp_id, up_seq, ref_seq, dn_seq, ','.join(alts_seq),
                        ref_lp, ','.join(('{:.2f}'.format(alt_lp)
                                          for alt_lp in alt_lps)),
                        'WITH_CONTEXT' if w_context else 'NO_CONTEXT')
    logger.debug(out_txt)
    return

def score_seq(tpost, seq, tpost_start=0, tpost_end=None,
              all_paths=False):
    """Score a section of log transition posteriors against a proposed sequence
    using a global mapping.
    :param tpost: `ndarray` containing log transition posteriors to be scored
    :param seq: `ndarray` containing integers encoding proposed sequence
    :param tpost_start: start position within post (Default: 0)
    :param tpost_end: end position within post (Default: full posterior)
    :param all_paths: boolean to produce the forwards all paths score
        (default Viterbi best path)
    """
    seq = seq.astype(np.uintp)
    if tpost_end is None:
        tpost_end = post.shape[0]
    if seq.shape[0] >= tpost_end - tpost_start:
        raise mh.MegaError('Mapped signal too short for proposed sequence.')

    score = decode.score_seq(tpost, seq, tpost_start, tpost_end, all_paths)
    if np.isnan(score):
        raise mh.MegaError('Score computation error (likely  memory error).')

    return score

def call_read_snps(
        snps_data, read_ref_pos, strand_read_np_ref_seq, rl_cumsum, r_to_q_poss,
        r_post, post_mapped_start):
    if read_ref_pos.end - read_ref_pos.start <= 2 * snps_data.edge_buffer:
        raise mh.MegaError('Mapped region too short for variant calling.')

    # convert to forward strand sequence in order to annotate with variants
    read_ref_fwd_seq = (strand_read_np_ref_seq if read_ref_pos.strand == 1 else
                        mh.revcomp_np(strand_read_np_ref_seq))
    # call all snps overlapping this read
    r_snp_calls = []
    logger = logging.get_logger('per_read_snps')
    read_cached_scores = {}
    read_variants = snps_data.fetch_read_variants(
        read_ref_pos, read_ref_fwd_seq)
    filt_read_variants = []
    # first pass over variants assuming the reference ground truth
    # (not including context variants)
    for (np_s_snp_ref_seq, np_s_snp_alt_seqs, np_s_context_seqs,
         s_ref_start, s_ref_end, variant) in snps_data.iter_snps(
             read_variants, read_ref_pos, read_ref_fwd_seq, context_max_dist=0):
        blk_start  = rl_cumsum[r_to_q_poss[s_ref_start]]
        blk_end = rl_cumsum[r_to_q_poss[s_ref_end]]
        if blk_end - blk_start <= max(
                len(up_seq) + len(dn_seq)
                for up_seq, dn_seq in np_s_context_seqs) + max(
                        np_s_snp_ref_seq.shape[0], max(
                            snp_alt_seq.shape[0]
                            for snp_alt_seq in np_s_snp_alt_seqs)):
            # no valid mapping over large inserted query bases
            # i.e. need as many "events/strides" as bases for valid mapping
            continue

        np_ref_seq = np.concatenate([
            np_s_context_seqs[0][0], np_s_snp_ref_seq, np_s_context_seqs[0][1]])
        loc_ref_lp = score_seq(
            r_post, np_ref_seq, post_mapped_start + blk_start,
            post_mapped_start + blk_end, snps_data.all_paths)

        loc_alt_lps = []
        loc_alt_llrs = []
        if _DEBUG_PER_READ:
            loc_contexts_alts_lps = []
        for np_s_snp_alt_seq, var_alt_seq in zip(
                np_s_snp_alt_seqs, variant.alts):
            np_alt_seq = np.concatenate([
                np_s_context_seqs[0][0], np_s_snp_alt_seq,
                np_s_context_seqs[0][1]])
            loc_alt_lp = score_seq(
                r_post, np_alt_seq, post_mapped_start + blk_start,
                post_mapped_start + blk_end, snps_data.all_paths)
            loc_alt_lps.append(loc_alt_lp)
            if _DEBUG_PER_READ:
                loc_contexts_alts_lps.append(np.array([loc_alt_lp,]))
            # calibrate log probs
            loc_alt_llrs.append(snps_data.calibrate_llr(
                loc_ref_lp - loc_alt_lp, variant.ref, var_alt_seq))

        # due to calibration mutli-allelic log likelihoods could result in
        # inferred negative reference likelihood, so re-normalize here
        loc_alt_log_ps = calibration.compute_log_probs(np.array(loc_alt_llrs))

        if _DEBUG_PER_READ:
            write_per_read_debug(
                variant.start, variant.id, read_ref_pos,
                np_s_snp_ref_seq, np_s_snp_alt_seqs, np_s_context_seqs,
                np.array([loc_ref_lp,]), loc_contexts_alts_lps, False, logger)

        if sum(np.exp(loc_alt_log_ps)) >= snps_data.context_min_alt_prob:
            filt_read_variants.append(variant)
            read_cached_scores[(variant.id, variant.start, variant.stop)] = (
                loc_ref_lp, loc_alt_lps)
        else:
            r_snp_calls.append((
                variant.ref_start, loc_alt_log_ps, variant.ref,
                variant.alts, variant.id, variant.start,
                variant.start + variant.np_ref.shape[0]))

    # second round for variants with some evidence for alternative alleles
    # process with other potential variants as context
    for (np_s_snp_ref_seq, np_s_snp_alt_seqs, np_s_context_seqs,
         s_ref_start, s_ref_end, variant) in snps_data.iter_snps(
             filt_read_variants, read_ref_pos, read_ref_fwd_seq):
        ref_cntxt_ref_lp, ref_cntxt_alt_lps = read_cached_scores[(
            variant.id, variant.start, variant.stop)]

        blk_start  = rl_cumsum[r_to_q_poss[s_ref_start]]
        blk_end = rl_cumsum[r_to_q_poss[s_ref_end]]
        if blk_end - blk_start <= max(
                len(up_seq) + len(dn_seq)
                for up_seq, dn_seq in np_s_context_seqs) + max(
                        np_s_snp_ref_seq.shape[0], max(
                            snp_alt_seq.shape[0]
                            for snp_alt_seq in np_s_snp_alt_seqs)):
            # if some context sequences are too long for signal
            # just use cached lps
            # TODO could also filter out invalid context sequences
            r_snp_calls.append((
                variant.start, ref_cntxt_alt_lps, variant.ref,
                variant.alts, variant.id, variant.start,
                variant.start + variant.np_ref.shape[0]))
            continue

        # skip first (reference) context seq as this was cached
        ref_context_seqs = (
            np.concatenate([up_context_seq, np_s_snp_ref_seq, dn_context_seq])
            for up_context_seq, dn_context_seq in np_s_context_seqs[1:])
        loc_contexts_ref_lps = np.array([ref_cntxt_ref_lp] + [score_seq(
            r_post, ref_seq, post_mapped_start + blk_start,
            post_mapped_start + blk_end, snps_data.all_paths)
                                         for ref_seq in ref_context_seqs])
        loc_ref_lp = logsumexp(loc_contexts_ref_lps)

        loc_alt_llrs = []
        if _DEBUG_PER_READ:
            loc_contexts_alts_lps = []
        for np_s_snp_alt_seq, var_alt_seq, ref_cntxt_alt_lp in zip(
                np_s_snp_alt_seqs, variant.alts, ref_cntxt_alt_lps):
            alt_context_seqs = (
                np.concatenate([
                    up_context_seq, np_s_snp_alt_seq, dn_context_seq])
                for up_context_seq, dn_context_seq in np_s_context_seqs[1:])
            loc_contexts_alt_lps = np.array([ref_cntxt_alt_lp,] + [
                score_seq(r_post, alt_seq, post_mapped_start + blk_start,
                          post_mapped_start + blk_end, snps_data.all_paths)
                for alt_seq in alt_context_seqs])
            loc_alt_lp = logsumexp(loc_contexts_alt_lps)
            if _DEBUG_PER_READ:
                loc_contexts_alts_lps.append(loc_contexts_alt_lps)
            # calibrate log probs
            loc_alt_llrs.append(snps_data.calibrate_llr(
                loc_ref_lp - loc_alt_lp, variant.ref, var_alt_seq))

        # due to calibration mutli-allelic log likelihoods could result in
        # inferred negative reference likelihood, so re-normalize here
        loc_alt_log_ps = calibration.compute_log_probs(np.array(loc_alt_llrs))

        if _DEBUG_PER_READ:
            write_per_read_debug(
                variant.start, variant.id, read_ref_pos,
                np_s_snp_ref_seq, np_s_snp_alt_seqs, np_s_context_seqs,
                loc_contexts_ref_lps, loc_contexts_alts_lps, True, logger)

        r_snp_calls.append((
            variant.ref_start, loc_alt_log_ps, variant.ref,
            variant.alts, variant.id, variant.start,
            variant.start + variant.np_ref.shape[0]))

    # re-sort variants after adding context-included computations
    return sorted(r_snp_calls, key=lambda x: x[0])


###############################
##### Per-read SNP Output #####
###############################

def log_prob_to_phred(log_prob):
    with np.errstate(divide='ignore'):
        return -10 * np.log10(1 - np.exp(log_prob))

def simplify_snp_seq(ref_seq, alt_seq):
    trim_before = trim_after = 0
    while (len(ref_seq) > 0 and len(alt_seq) > 0 and
           ref_seq[0] == alt_seq[0]):
        trim_before += 1
        ref_seq = ref_seq[1:]
        alt_seq = alt_seq[1:]
    while (len(ref_seq) > 0 and len(alt_seq) > 0 and
           ref_seq[-1] == alt_seq[-1]):
        trim_after += 1
        ref_seq = ref_seq[:-1]
        alt_seq = alt_seq[:-1]

    return ref_seq, alt_seq, trim_before, trim_after

def iter_non_overlapping_snps(r_snp_calls):
    def get_max_prob_allele_snp(snp_grp):
        """ For overlapping SNPs return the snp with the highest probability
        single allele as this one will be added to the reference sequence.

        More complex chained SNPs could be handled, but are not here.
        For example, a 5 base deletion covering 2 single base swap SNPs could
        validly result in 2 alternative single base swap alleles, but the
        logic here would only allow one of those alternatives since they
        are covered by the same reference deletion. There are certainly many
        more edge cases than this and each one would require specific logic.
        This likely covers the majority of valid cases and limiting to
        50 base indels by default limits the scope of this issue.
        """
        most_prob_snp = None
        for snp_data in snp_grp:
            ref_lp = np.log1p(-np.exp(snp_data[1]).sum())
            snp_max_lp = max(ref_lp, max(snp_data[1]))
            if most_prob_snp is None or snp_max_lp > most_prob_snp[0]:
                most_prob_snp = (snp_max_lp, ref_lp, snp_data)

        _, ref_lp, (snp_pos, alt_lps, snp_ref_seq,
                    snp_alt_seqs, _, _, _) = most_prob_snp
        return snp_pos, alt_lps, snp_ref_seq, snp_alt_seqs, ref_lp


    if len(r_snp_calls) == 0: return
    r_snp_calls_iter = iter(r_snp_calls)
    # initialize snp_grp with first snp
    snp_data = next(r_snp_calls_iter)
    prev_snp_end = snp_data[0] + len(snp_data[2])
    snp_grp = [snp_data]
    for snp_data in sorted(r_snp_calls_iter, key=itemgetter(0)):
        if snp_data[0] < prev_snp_end:
            prev_snp_end = max(snp_data[0] + len(snp_data[2]), prev_snp_end)
            snp_grp.append(snp_data)
        else:
            yield get_max_prob_allele_snp(snp_grp)
            prev_snp_end = snp_data[0] + len(snp_data[2])
            snp_grp = [snp_data]

    # yeild last snp grp data
    yield get_max_prob_allele_snp(snp_grp)
    return

def annotate_snps(r_start, ref_seq, r_snp_calls, strand):
    """ Annotate reference sequence with called snps.

    Note: Reference sequence is in read orientation and snp calls are in
    genome coordiates.
    """
    snp_seqs, snp_quals, snp_cigar = [], [], []
    prev_pos, curr_match = 0, 0
    # ref_seq is read-centric so flop order to process snps in genomic order
    if strand == -1:
        ref_seq = ref_seq[::-1]
    for (snp_pos, alt_lps, snp_ref_seq, snp_alt_seqs,
         ref_lp) in iter_non_overlapping_snps(r_snp_calls):
        prev_len = snp_pos - r_start - prev_pos
        # called canonical
        if ref_lp >= max(alt_lps):
            snp_seqs.append(
                ref_seq[prev_pos:snp_pos - r_start + len(snp_ref_seq)])
            snp_quals.extend(
                ([WHATSHAP_MAX_QUAL] * prev_len) +
                ([min(log_prob_to_phred(ref_lp), WHATSHAP_MAX_QUAL)] *
                 len(snp_ref_seq)))
            curr_match += prev_len + len(snp_ref_seq)
        else:
            alt_seq = snp_alt_seqs[np.argmax(alt_lps)]
            # complement since ref_seq is complement seq
            # (not reversed; see loop init)
            read_alt_seq = alt_seq if strand == 1 else mh.comp(alt_seq)
            snp_seqs.append(ref_seq[prev_pos:snp_pos - r_start] + read_alt_seq)
            snp_quals.extend(
                ([WHATSHAP_MAX_QUAL] * prev_len) +
                ([min(log_prob_to_phred(max(alt_lps)), WHATSHAP_MAX_QUAL)] *
                 len(alt_seq)))

            # add cigar information for snp or indel
            t_ref_seq, t_alt_seq, t_before, t_after = simplify_snp_seq(
                snp_ref_seq, alt_seq)
            curr_match += t_before
            snp_cigar.append((7, curr_match + prev_len))
            if len(t_alt_seq) == len(t_ref_seq):
                snp_cigar.append((8, len(t_alt_seq)))
            elif len(t_alt_seq) > len(t_ref_seq):
                # left justify mismatch bases in complex insertion
                if len(t_ref_seq) != 0:
                    snp_cigar.append((8, len(t_ref_seq)))
                snp_cigar.append((1, len(t_alt_seq) - len(t_ref_seq)))
            else:
                # left justify mismatch bases in complex deletion
                if len(t_alt_seq) != 0:
                    snp_cigar.append((8, len(t_alt_seq)))
                snp_cigar.append((2, len(t_ref_seq) - len(t_alt_seq)))
            curr_match = t_after
        prev_pos = snp_pos - r_start + len(snp_ref_seq)

    snp_seqs.append(ref_seq[prev_pos:])
    snp_seq = ''.join(snp_seqs)
    if strand == -1:
        snp_seq = snp_seq[::-1]
    len_remain = len(ref_seq) - prev_pos
    snp_quals.extend([WHATSHAP_MAX_QUAL] * len_remain)
    if strand == -1:
        snp_quals = snp_quals[::-1]
    snp_quals = list(map(int, snp_quals))
    snp_cigar.append((7, len_remain + curr_match))
    if strand == -1:
        snp_cigar = snp_cigar[::-1]

    return snp_seq, snp_quals, snp_cigar

def _get_snps_queue(
        snps_q, snps_conn, snps_db_fn, snps_txt_fn, db_safety, pr_refs_fn,
        pr_ref_filts, whatshap_map_fn, ref_names_and_lens, ref_fn):
    def write_whatshap_alignment(
            read_id, snp_seq, snp_quals, chrm, strand, r_st, snp_cigar):
        a = pysam.AlignedSegment()
        a.query_name = read_id
        a.flag = 0 if strand == 1 else 16
        a.reference_id = whatshap_map_fp.get_tid(chrm)
        a.reference_start = r_st
        a.template_length = len(snp_seq)
        a.mapping_quality = WHATSHAP_MAX_QUAL
        a.set_tags([('RG', WHATSHAP_RG_ID)])

        # convert to reference based sequence
        if strand == -1:
            snp_seq = mh.revcomp(snp_seq)
            snp_quals = snp_quals[::-1]
            snp_cigar = snp_cigar[::-1]
        a.query_sequence = snp_seq
        a.query_qualities = array('B', snp_quals)
        a.cigartuples = snp_cigar
        whatshap_map_fp.write(a)

        return

    def get_snp_call(
            r_snp_calls, read_id, chrm, strand, r_start, ref_seq, read_len,
            q_st, q_en, cigar):
        # note strand is +1 for fwd or -1 for rev
        snps_db.executemany(ADDMANY_SNPS, [
            (read_id, chrm, strand, pos, alt_lp,
             snp_ref_seq, snp_alt_seq, snp_id, test_start, test_end)
            for pos, alt_lps, snp_ref_seq, snp_alt_seqs, snp_id,
            test_start, test_end in r_snp_calls
            for alt_lp, snp_alt_seq in zip(alt_lps, snp_alt_seqs)])
        if snps_txt_fp is not None and len(r_snp_calls) > 0:
            snps_txt_fp.write('\n'.join((
                ('\t'.join('{}' for _ in field_names)).format(
                    read_id, chrm, strand, pos,
                    np.log1p(-np.exp(alt_lps).sum()), alt_lp,
                    snp_ref_seq, snp_alt_seq, snp_id)
                for pos, alt_lps, snp_ref_seq, snp_alt_seqs, snp_id,
                test_start, test_end in r_snp_calls
                for alt_lp, snp_alt_seq in zip(alt_lps, snp_alt_seqs))) + '\n')
            snps_txt_fp.flush()
        if do_ann_snps:
            if not mapping.read_passes_filters(
                    pr_ref_filts, read_len, q_st, q_en, cigar):
                return
            snp_seq, snp_quals, snp_cigar = annotate_snps(
                r_start, ref_seq, r_snp_calls, strand)
            if pr_refs_fn is not None:
                pr_refs_fp.write('>{}\n{}\n'.format(read_id, snp_seq))
                pr_refs_fp.flush()
            if whatshap_map_fn is not None:
                write_whatshap_alignment(
                    read_id, snp_seq, snp_quals, chrm, strand, r_start,
                    snp_cigar)

        return


    logger = logging.get_logger('snps_getter')
    snps_db = sqlite3.connect(snps_db_fn)
    if db_safety < 2:
        snps_db.execute(SET_ASYNC_MODE)
    if db_safety < 1:
        snps_db.execute(SET_NO_ROLLBACK_MODE)
    snps_db.execute(CREATE_SNPS_TBLS)
    if snps_txt_fn is None:
        snps_txt_fp = None
    else:
        snps_txt_fp = open(snps_txt_fn, 'w')
        field_names = (
            'read_id', 'chrm', 'strand', 'pos', 'ref_log_prob', 'alt_log_prob',
            'ref_seq', 'alt_seq', 'snp_id')
        snps_txt_fp.write('\t'.join(field_names) + '\n')

    if pr_refs_fn is not None:
        pr_refs_fp = open(pr_refs_fn, 'w')

    if whatshap_map_fn is not None:
        _, map_fmt = os.path.splitext(whatshap_map_fn)
        if map_fmt == '.bam': w_mode = 'wb'
        elif map_fmt == '.cram': w_mode = 'wc'
        elif map_fmt == '.sam': w_mode = 'w'
        else:
            raise mh.MegaError('Invalid mapping output format')
        header = {
            'HD': {'VN': '1.4'},
            'SQ': [{'LN': ref_len, 'SN': ref_name}
                   for ref_name, ref_len in sorted(zip(*ref_names_and_lens))],
            'RG': [{'ID':WHATSHAP_RG_ID, 'SM':SAMPLE_NAME},]}
        whatshap_map_fp = pysam.AlignmentFile(
            whatshap_map_fn, w_mode, header=header, reference_filename=ref_fn)

    do_ann_snps = whatshap_map_fn is not None or pr_refs_fn is not None

    while True:
        try:
            r_snp_calls, (read_id, chrm, strand, r_start, ref_seq, read_len,
                          q_st, q_en, cigar) = snps_q.get(block=False)
        except queue.Empty:
            if snps_conn.poll():
                break
            sleep(0.1)
            continue
        try:
            get_snp_call(
                r_snp_calls, read_id, chrm, strand, r_start, ref_seq, read_len,
                q_st, q_en, cigar)
        except Exception as e:
            logger.debug((
                'Error processing variant output for read: {}\nSet' +
                ' _RAISE_VARIANT_PROCESSING_ERRORS in megalodon/snps.py to ' +
                'see full error.\nError type: {}').format(read_id, str(e)))
            if _RAISE_VARIANT_PROCESSING_ERRORS: raise

    while not snps_q.empty():
        r_snp_calls, (read_id, chrm, strand, r_start, ref_seq, read_len,
                      q_st, q_en, cigar) = snps_q.get(block=False)
        try:
            get_snp_call(
                r_snp_calls, read_id, chrm, strand, r_start, ref_seq, read_len,
                q_st, q_en, cigar)
        except Exception as e:
            logger.debug((
                'Error processing variant output for read: {}\nSet' +
                ' _RAISE_VARIANT_PROCESSING_ERRORS in megalodon/snps.py to ' +
                'see full error.\nError type: {}').format(read_id, str(e)))
            if _RAISE_VARIANT_PROCESSING_ERRORS: raise
    if snps_txt_fp is not None: snps_txt_fp.close()
    if pr_refs_fn is not None: pr_refs_fp.close()
    if whatshap_map_fn is not None: whatshap_map_fp.close()
    snps_db.execute(CREATE_SNPS_IDX)
    snps_db.commit()
    snps_db.close()

    return


######################
##### VCF Reader #####
######################

class SnpData(object):
    def check_vars_match_ref(
            self, vars_idx, contigs, aligner, num_contigs=5,
            num_sites_per_contig=50):
        """ Validate that the reference sequences in the variant file matches
        a reference sequence file.
        """
        for contig in contigs[:num_contigs]:
            for var_data in list(vars_idx.fetch(contig))[:num_sites_per_contig]:
                ref_seq = aligner.seq(contig, var_data.start, var_data.stop)
                if ref_seq != var_data.ref:
                    # variant reference sequence does not match reference
                    logger = logging.get_logger()
                    logger.debug((
                        'Reference sequence does not match variant reference ' +
                        'sequence at {} expected "{}" got "{}"').format(
                            snp_ref_pos, var_data.ref, ref_seq))
                    return False

        return True

    def __init__(
            self, variant_fn, max_indel_size, all_paths,
            write_snps_txt, context_bases, snps_calib_fn=None,
            call_mode=DIPLOID_MODE, do_pr_ref_snps=False, aligner=None,
            keep_snp_fp_open=False, do_validate_reference=True,
            edge_buffer=mh.DEFAULT_EDGE_BUFFER,
            context_min_alt_prob=mh.DEFAULT_CONTEXT_MIN_ALT_PROB):
        logger = logging.get_logger('snps')
        self.max_indel_size = max_indel_size
        self.all_paths = all_paths
        self.write_snps_txt = write_snps_txt
        self.snps_calib_fn = snps_calib_fn
        self.calib_table = calibration.SnpCalibrator(self.snps_calib_fn)
        self.context_bases = context_bases
        if len(self.context_bases) != 2:
            raise mh.MegaError(
                'Must provide 2 context bases values (for single base SNPs ' +
                'and indels).')
        self.call_mode = call_mode
        self.do_pr_ref_snps = do_pr_ref_snps
        self.edge_buffer = edge_buffer
        self.context_min_alt_prob = context_min_alt_prob
        self.variant_fn = variant_fn
        self.variants_idx = None
        if self.variant_fn is None:
            return

        logger.info('Loading variants.')
        vars_idx = pysam.VariantFile(self.variant_fn)
        try:
            contigs = list(vars_idx.header.contigs.keys())
            vars_idx.fetch(next(iter(contigs)), 0, 0)
        except ValueError:
            logger.warn(
                'Variants file must be indexed. Performing indexing now.')
            vars_idx.close()
            self.variant_fn = index_variants(self.variant_fn)
            vars_idx = pysam.VariantFile(self.variant_fn)
        if aligner is None:
            raise mh.MegaError(
                'Must provide aligner if SNP filename is provided')
        if len(set(aligner.ref_names_and_lens[0]).intersection(contigs)) == 0:
            raise mh.MegaError((
                'Reference and variant files contain no chromosomes/contigs ' +
                'in common.\n\t\tFirst 3 reference contigs:\t{}\n\t\tFirst 3 ' +
                'variant file contigs:\t{}').format(
                    ', '.join(aligner.ref_names_and_lens[0][:3]),
                    ', '.join(contigs[:3])))
        if do_validate_reference and not self.check_vars_match_ref(
                vars_idx, contigs, aligner):
            raise mh.MegaError(
                'Reference sequence file does not match reference sequence ' +
                'from variants file.')

        if keep_snp_fp_open:
            self.variants_idx = vars_idx
        else:
            vars_idx.close()
            self.variants_idx = None

        return

    @property
    def substitution_context(self):
        return self.context_bases[0]
    @property
    def indel_context(self):
        return self.context_bases[1]

    def reopen_variant_index(self):
        if self.variant_fn is not None:
            self.variants_idx = pysam.VariantFile(self.variant_fn)
        return

    @staticmethod
    def compute_variant_distance(var1, var2):
        # if the variants overlap return None
        if not (var1.start >= var2.stop or
                var2.start >= var1.stop):
            return None
        if var1.start >= var2.stop:
            return var1.start - var2.stop
        return var2.start - var1.stop

    @staticmethod
    def any_variants_overlap(variants):
        for var1, var2 in combinations(variants, 2):
            if not (var1.start >= var2.stop or
                    var2.start >= var1.stop):
                return True
        return False

    @staticmethod
    def annotate_context_seqs(
            context_vars, up_context_seq, dn_context_seq, context_ref_start,
            context_rel_var_start, context_rel_var_end):
        # annotate upstream sequence
        ann_up_seq = []
        prev_end = 0
        for context_var, np_alt_seq in context_vars:
            if context_var.stop - context_ref_start > context_rel_var_start:
                break
            ann_up_seq.extend((up_context_seq[
                prev_end:context_var.start - context_ref_start], np_alt_seq))
            prev_end = context_var.stop - context_ref_start
        ann_up_seq.append(up_context_seq[prev_end:])
        ann_up_seq = np.concatenate(ann_up_seq)

        # annotate downstream sequence
        ann_dn_seq = []
        prev_end = 0
        for context_var, np_alt_seq in context_vars:
            if context_var.start - context_ref_start < context_rel_var_end:
                continue
            ann_dn_seq.extend((dn_context_seq[
                prev_end:context_var.start - context_ref_start -
                context_rel_var_end], np_alt_seq))
            prev_end = (context_var.stop - context_ref_start -
                        context_rel_var_end)
        ann_dn_seq.append(dn_context_seq[prev_end:])
        ann_dn_seq = np.concatenate(ann_dn_seq)

        return ann_up_seq, ann_dn_seq

    @staticmethod
    def iter_variant_combos_by_distance(variant, context_variants):
        """ Yield combinations of variants ordered by inclusion of variants
        closer to the variant of interest first.
        """
        def iter_alt_variant_seqs(variants):
            """ Single variants can have multiple alternative sequences,
            so iterate over those here
            """
            vars_w_alts = [[(var, np_alt) for np_alt in var.np_alts]
                           for var in variants]
            for vars_w_alt in product(*vars_w_alts):
                yield sorted(vars_w_alt, key=lambda x: x[0].start)
            return

        dist_vars = defaultdict(list)
        for context_var in context_variants:
            var_dist = SnpData.compute_variant_distance(variant, context_var)
            if var_dist is not None:
                dist_vars[var_dist].append(context_var)

        used_vars = []
        # sort list by distance to a variant in order to report most
        # proximal variants first
        for var_dist, dist_context_vars in sorted(dist_vars.items()):
            # loop over numbers of current distance variants (1 or more) and
            # previously used variants (0 or more)
            for n_dist_vars, n_used_vars in product(
                    range(1, len(dist_context_vars) + 1),
                    range(len(used_vars) + 1)):
                # loop over actual selection of variants
                for curr_vars, used_vars_i in product(
                        combinations(dist_context_vars, n_dist_vars),
                        combinations(used_vars, n_used_vars)):
                    # including selected alt seq
                    for vars_w_alt in iter_alt_variant_seqs(
                            list(curr_vars) + list(used_vars_i)):
                        yield vars_w_alt
            used_vars.extend(dist_context_vars)

        return

    @staticmethod
    def iter_context_variants(variants, context_max_dist):
        """ Iterate variants as well as variants within context_max_dist
        """
        vars_iter = iter(variants)
        def next_var_or_none():
            try:
                return next(vars_iter)
            except StopIteration:
                return None

        curr_vars = [next(vars_iter),]
        next_var = next(vars_iter)
        if next_var is None:
            if curr_vars[0] is not None:
                yield curr_vars[0], []
            return
        curr_var_idx = 0

        while next_var is not None:
            if curr_var_idx == len(curr_vars):
                curr_vars.append(next_var)
                next_var = next_var_or_none()
            curr_var = curr_vars[curr_var_idx]
            # add relevant variants
            while (next_var is not None and
                   next_var.start - context_max_dist <= curr_var.stop):
                curr_vars.append(next_var)
                next_var = next_var_or_none()
            # remove variants that end before the variant of interest
            n_vars_removed = sum(
                var.stop + context_max_dist < curr_var.start + 1
                for var in curr_vars)
            curr_vars = [
                var for var in curr_vars
                if var.stop + context_max_dist >= curr_var.start + 1]
            curr_var_idx -= n_vars_removed - 1
            # yeild variants in range of current variant
            yield curr_var, curr_vars

        # yield final vars from the read
        while len(curr_vars) < curr_var_idx:
            curr_var = curr_vars[curr_var_idx]
            # remove variants that end before the variant of interest
            n_vars_removed = sum(
                var.stop + context_max_dist < curr_var.start
                for var in curr_vars)
            curr_vars = [
                var for var in curr_vars
                if var.stop + context_max_dist >= curr_var.start]
            curr_var_idx -= n_vars_removed - 1
            # yeild variants in range of current variant
            yield curr_var, curr_vars

        return

    @staticmethod
    def add_indel_context_base(
            np_ref_seq, np_alt_seqs, var_start, read_ref_fwd_seq,
            read_ref_pos):
        if np_ref_seq.shape[0] == 0 or any(
                np_alt.shape[0] == 0 for np_alt in np_alt_seqs):
            upstrm_base = mh.ALPHABET[read_ref_fwd_seq[
                var_start - read_ref_pos.start - 1]]
            var_start -= 1
        else:
            upstrm_base = ''
        ref_seq = upstrm_base + mh.int_to_seq(np_ref_seq)
        alt_seqs = tuple((upstrm_base + mh.int_to_seq(np_alt)
                          for np_alt in np_alt_seqs))
        return ref_seq, alt_seqs, var_start

    def merge_variants(
            self, grouped_read_vars, read_ref_fwd_seq, read_ref_pos):
        """ Merge atomized variants into multi-allelic sites.
        if this is not done, allele probabilities will not be normalized
        correctly.
        """
        variants = []
        for _, site_vars in sorted(grouped_read_vars.items()):
            if len(site_vars) == 1:
                # add upstream seq to simple indels
                (out_ref_seq, out_alt_seqs,
                 out_start) = self.add_indel_context_base(
                     site_vars[0].np_ref, site_vars[0].np_alts,
                     site_vars[0].start, read_ref_fwd_seq, read_ref_pos)
                site_var = site_vars[0]._replace(
                    ref=out_ref_seq, alts=out_alt_seqs, ref_start=out_start)
                variants.append(site_var)
                continue

            # join all valid ids
            # skip None ids ('.' in VCF)
            site_var_ids = set(
                var_id for var in site_vars
                if var.id is not None and var.id != '.'
                for var_id in var.id.split(';'))
            # if all ids are None leave id as None
            site_var_ids = (';'.join(sorted(site_var_ids))
                            if len(site_var_ids) > 0 else None)
            # assuming atomized variants with single alt
            # np arrays aren't hash-able, so test for equality manually here
            np_alts = []
            for var in site_vars:
                if any((len(var.np_alts[0]) == len(added_np_alt) and
                        np.all(var.np_alts[0] == added_np_alt))
                       for added_np_alt in np_alts):
                    continue
                np_alts.append(var.np_alts[0])

            # add upstream seq to simple indels
            out_ref_seq, out_alt_seqs, out_start = self.add_indel_context_base(
                site_vars[0].np_ref, np_alts, site_vars[0].start,
                read_ref_fwd_seq, read_ref_pos)

            variants.append(VARIANT_DATA(
                np_ref=site_vars[0].np_ref, np_alts=np_alts,
                id=site_var_ids, chrom=site_vars[0].chrom,
                start=site_vars[0].start, stop=site_vars[0].stop,
                ref=out_ref_seq, alts=out_alt_seqs, ref_start=out_start))

        return variants

    @staticmethod
    def expand_ambig_variant(
            np_ref_seq, np_alt_seq, var_start, read_ref_fwd_seq, read_ref_pos):
        # don't try to expand complex variants
        if np_ref_seq.shape[0] != 0 and np_alt_seq.shape[0] != 0:
            return np_ref_seq, np_alt_seq, var_start
        elif np_ref_seq.shape[0] == 0:
            # expand ambiguous insertion sequence
            indel_seq = np_alt_seq
        else:
            # expand ambiguous deletion sequence
            indel_seq = np_ref_seq

        expand_dnstrm_seq = []
        expand_start = expand_pos = (
            var_start - read_ref_pos.start + np_ref_seq.shape[0])
        # if shifting the insertion one base downstream would be
        # equivalent add this position
        while (expand_pos < read_ref_fwd_seq.shape[0] and
               indel_seq[(expand_pos - expand_start) %
                          indel_seq.shape[0]] ==
               read_ref_fwd_seq[expand_pos]):
            expand_dnstrm_seq.append(read_ref_fwd_seq[expand_pos])
            expand_pos += 1
        if expand_pos == read_ref_fwd_seq.shape[0]:
            raise mh.MegaError(
                'Variant is ambiguous up to the end of the read.')

        # mirror for upstream ambiguous sequence
        expand_upstrm_seq = []
        expand_start = expand_pos = var_start - read_ref_pos.start - 1
        # if shifting the insertion one base upstream would be
        # equivalent add this position
        while (expand_pos >= 0 and
               indel_seq[indel_seq.shape[0] - 1 - (
                   (expand_start - expand_pos) % indel_seq.shape[0])] ==
               read_ref_fwd_seq[expand_pos]):
            expand_upstrm_seq.insert(0, read_ref_fwd_seq[expand_pos])
            expand_pos -= 1
        if expand_pos == -1:
            raise mh.MegaError(
                'Variant is ambiguous up to the start of the read.')

        np_ref_seq = np.concatenate([
            expand_upstrm_seq, np_ref_seq, expand_dnstrm_seq]).astype(np.uintp)
        np_alt_seq = np.concatenate([
            expand_upstrm_seq, np_alt_seq, expand_dnstrm_seq]).astype(np.uintp)
        var_start -= len(expand_upstrm_seq)
        return np_ref_seq, np_alt_seq, var_start

    def iter_atomized_variants(
            self, var, np_ref_seq, np_alt_seq, read_ref_fwd_seq, read_ref_pos):
        # substitutions
        if np_alt_seq.shape[0] == np_ref_seq.shape[0]:
            # convert all substitutions into single base substitutions
            for sub_offset, (np_ref_base, np_alt_base) in enumerate(zip(
                    np_ref_seq, np_alt_seq)):
                if np_ref_base == np_alt_base: continue
                yield (
                    (var.start + sub_offset, var.start + sub_offset + 1),
                    VARIANT_DATA(
                        np_ref=np.array([np_ref_base], dtype=np.uintp),
                        np_alts=(np.array([np_alt_base], dtype=np.uintp),),
                        id=var.id, chrom=var.chrom,
                        start=var.start + sub_offset,
                        stop=var.start + sub_offset + 1))
        else:
            # skip large indels
            if np.abs(np_ref_seq.shape[0] -
                      np_alt_seq.shape[0]) > self.max_indel_size:
                return

            # trim context bases from seq
            np_ref_seq, np_alt_seq, start_trim, _ = simplify_snp_seq(
                np_ref_seq, np_alt_seq)
            var_start = var.start + start_trim
            try:
                # expand seqs to include ambiguous locations
                np_ref_seq, np_alt_seq, var_start = self.expand_ambig_variant(
                    np_ref_seq, np_alt_seq, var_start, read_ref_fwd_seq,
                    read_ref_pos)
            except mh.MegaError:
                # if variant is ambiguous to the end of the read, then skip it
                return

            yield ((var_start, var_start + np_ref_seq.shape[0]),
                   VARIANT_DATA(
                       np_ref=np_ref_seq, np_alts=(np_alt_seq,), id=var.id,
                       chrom=var.chrom, start=var_start,
                       stop=var_start + np_ref_seq.shape[0]))
        return

    def atomize_variants(self, fetch_res, read_ref_fwd_seq, read_ref_pos):
        grouped_read_vars = defaultdict(list)
        for var in fetch_res:
            # fetch results include any overlap where only inclusive overlap
            # are valid here
            if (var.stop + self.edge_buffer > read_ref_pos.end or
                var.start - self.edge_buffer < read_ref_pos.start):
                continue

            np_ref_seq = mh.seq_to_int(var.ref)
            for alt_seq in var.alts:
                np_alt_seq = mh.seq_to_int(alt_seq)
                for start_stop, atom_var in self.iter_atomized_variants(
                        var, np_ref_seq, np_alt_seq, read_ref_fwd_seq,
                        read_ref_pos):
                    grouped_read_vars[start_stop].append(atom_var)
        return grouped_read_vars

    def fetch_read_variants(self, read_ref_pos, read_ref_fwd_seq):
        try:
            fetch_res = self.variants_idx.fetch(
                read_ref_pos.chrm, read_ref_pos.start + self.edge_buffer,
                read_ref_pos.end - self.edge_buffer)
        except ValueError:
            raise mh.MegaError('Mapped location not valid for variants file.')

        grouped_read_vars = self.atomize_variants(
            fetch_res, read_ref_fwd_seq, read_ref_pos)
        read_variants = self.merge_variants(
            grouped_read_vars, read_ref_fwd_seq, read_ref_pos)
        return read_variants

    def iter_snps(
            self, read_variants, read_ref_pos, read_ref_fwd_seq,
            max_contexts=16, context_max_dist=mh.CONTEXT_MAX_DIST):
        """Iterator over SNPs overlapping the read mapped position.

        Args:
            read_variants: List of variant objects (from fetch_read_variants)
            read_ref_pos: Read reference mapping position
                (`megalodon.mapping.MAP_POS`)
            read_ref_fwd_seq: Mapped reference sequence. Forward strand sequence
                no matter the mapped strand.
            max_contexts: Maximum number of context variant combinations to
                include around each variant.

        Yields:
            snp_ref_seq: Reference variant sequence on read strand
            snp_alt_seqs: Alternative variant sequences on read strand
            context_seqs: Sequences surrounding the variant on read strand
            context_start: Start of variant context in read coordinates
            context_end: End of variant context in read coordinates
            variant_ref: Reference variant sequence on reference
            variant_alts: Alternate variant sequences on reference
            variant_id: string idnentifier for the variant
            pos: variant position (0-based coordinate)

        SNPs within edge buffer of the end of the mapping will be ignored.

        If more than max_contexts snps exist within context_basss then only
        the max_contexts most proximal to the variant in question will be
        returned.
        """
        def revcomp_variant(context_seqs, np_var_ref, var_alts):
            rc_context_seqs = [
                (mh.revcomp_np(dn_context_seq), mh.revcomp_np(up_context_seq))
                for up_context_seq, dn_context_seq in context_seqs]
            return rc_context_seqs, mh.revcomp_np(np_var_ref), [
                mh.revcomp_np(np_alt) for np_alt in np_var_alts]

        def extract_variant_contexts(variant, context_vars):
            # compute various relative coordinates (in reference coordinates
            # not on read strand, to make it easier to work with variants)
            # select single base substitution or indel context width
            var_context_bases = self.substitution_context if all(
                variant.np_ref.shape[0] == np_alt.shape[0]
                for np_alt in variant.np_alts) else self.indel_context
            context_ref_start = variant.start - var_context_bases
            if context_ref_start < read_ref_pos.start:
                context_ref_start = read_ref_pos.start
            context_ref_end = variant.stop + var_context_bases
            if context_ref_end > read_ref_pos.end:
                context_ref_end = read_ref_pos.end
            context_read_start = context_ref_start - read_ref_pos.start
            context_read_end = context_ref_end - read_ref_pos.start

            context_rel_var_start = variant.start - context_ref_start
            context_rel_var_end = (
                context_rel_var_start + variant.np_ref.shape[0])

            context_ref_seq = read_ref_fwd_seq[
                context_read_start:context_read_end]
            # first context is always reference sequence
            up_context_seq = context_ref_seq[:context_rel_var_start]
            dn_context_seq = context_ref_seq[context_rel_var_end:]
            context_seqs = [(up_context_seq, dn_context_seq),]

            if max_contexts == 1 or len(context_vars) == 0:
                return (context_ref_start, context_read_start, context_read_end,
                        context_seqs)

            for context_vars_i in self.iter_variant_combos_by_distance(
                    variant, context_vars):
                if self.any_variants_overlap(list(zip(*context_vars_i))[0]):
                    continue
                context_seqs.append(self.annotate_context_seqs(
                    context_vars_i, up_context_seq, dn_context_seq,
                    context_ref_start, context_rel_var_start,
                    context_rel_var_end))
                if len(context_seqs) >= max_contexts:
                    break

            return (context_ref_start, context_read_start, context_read_end,
                    context_seqs)


        logger = logging.get_logger('snps')
        for variant, context_variants in self.iter_context_variants(
                read_variants, context_max_dist):
            (context_ref_start, context_read_start, context_read_end,
             np_context_seqs) = extract_variant_contexts(
                 variant, context_variants)

            # revcomp seqs for strand and convert to numpy arrays
            np_var_ref, np_var_alts = variant.np_ref, variant.np_alts
            if read_ref_pos.strand == -1:
                np_context_seqs, np_var_ref, np_var_alts = revcomp_variant(
                    np_context_seqs, np_var_ref, np_var_alts)
                context_read_start, context_read_end = (
                    read_ref_fwd_seq.shape[0] - context_read_end,
                    read_ref_fwd_seq.shape[0] - context_read_start)
            if np.concatenate([np_var_ref,] + list(np_var_alts) + [
                    seq for cntxt_seqs in np_context_seqs
                    for seq in cntxt_seqs]).max() > len(mh.ALPHABET):
                # some sequence contained invalid characters
                logger.debug(
                    'Invalid sequence encountered for variant ' +
                    '"{}" at {}:{}'.format(
                        variant.id, variant.chrom, variant.start))
                continue

            yield (
                np_var_ref, np_var_alts, np_context_seqs,
                context_read_start, context_read_end, variant)

        return

    def calibrate_llr(self, llr, var_ref_seq, var_alt_seq):
        return self.calib_table.calibrate_llr(
            llr, var_ref_seq, var_alt_seq)


######################
##### VCF Writer #####
######################

class Variant(object):
    """ Variant for entry into VcfWriter.
    Currently only handles a single sample.
    """
    def __init__(
            self, chrom, pos, ref, alts, id='.', qual='.', filter='.',
            info=None, sample_dict=None):
        self.chrom = chrom
        self.pos = int(pos)
        self.ref = ref.upper()
        self.alts = [alt.upper() for alt in alts]
        self.alt = ','.join(self.alts)
        self.id = '' if id is None else str(id)
        self.qual = qual
        self.filter = str(filter)
        if info is None:
            info = {}
        self.info_dict = info
        if sample_dict is None:
            sample_dict = OrderedDict()
        self.sample_dict = sample_dict
        return

    @property
    def _sorted_format_keys(self):
        sorted_keys = sorted(self.sample_dict.keys())
        if 'GT' in sorted_keys:
            sorted_keys.insert(0, sorted_keys.pop(sorted_keys.index('GT')))
        if 'LOG_PROBS' in sorted_keys:
            # move log probs to end of format field for easier human readability
            sorted_keys.append(sorted_keys.pop(sorted_keys.index('LOG_PROBS')))
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

    def add_haploid_probs(self, probs, gts):
        # phred scaled likelihoods
        with np.errstate(divide='ignore'):
            gl = np.maximum(mh.MIN_GL_VALUE, np.log10(probs))
        raw_pl = -10 * gl
        # "normalized" PL values stored as decsribed by VCF format
        # abs to remove negative 0 from file
        pl = np.abs(np.minimum(raw_pl - raw_pl.min(), mh.MAX_PL_VALUE))
        s_pl = np.sort(pl)

        # add sample tags
        self.add_sample_field('GT', gts[np.argmax(probs)])
        try:
            qual = int(np.around(np.minimum(raw_pl[0], mh.MAX_PL_VALUE)))
        except ValueError:
            logger = logging.get_logger()
            logger.debug(
                'NAN quality value encountered. gts:{}, probs:{}'.format(
                    str(gts), str(probs)))
            qual = mh.MAX_PL_VALUE
        self.qual = '{:.0f}'.format(qual) if qual > 0 else '.'
        self.add_sample_field('GQ', '{:.0f}'.format(np.around(s_pl[1])))
        self.add_sample_field(
            'GL', ','.join('{:.2f}' for _ in range(probs.shape[0])).format(*gl))
        self.add_sample_field(
            'PL', ','.join('{:.0f}' for _ in range(probs.shape[0])).format(
                *np.around(pl)))
        return

    def add_diploid_probs(self, probs, gts):
        # phred scaled likelihoods
        with np.errstate(divide='ignore'):
            gl = np.maximum(mh.MIN_GL_VALUE, np.log10(probs))
        raw_pl = -10 * gl
        # "normalized" PL values stored as decsribed by VCF format
        pl = np.minimum(raw_pl - raw_pl.min(), mh.MAX_PL_VALUE)
        s_pl = np.sort(pl)

        # add sample tags
        self.add_sample_field('GT', gts[np.argmax(probs)])
        try:
            qual = int(np.minimum(np.around(raw_pl[0]), mh.MAX_PL_VALUE))
        except ValueError:
            logger = logging.get_logger()
            logger.debug(
                'NAN quality value encountered. gts:{}, probs:{}'.format(
                    str(gts), str(probs)))
            qual = mh.MAX_PL_VALUE
        self.qual = '{:.0f}'.format(qual) if qual > 0 else '.'
        self.add_sample_field('GQ', '{:.0f}'.format(np.around(s_pl[1])))
        self.add_sample_field(
            'GL', ','.join('{:.2f}' for _ in range(probs.shape[0])).format(*gl))
        self.add_sample_field(
            'PL', ','.join('{:.0f}' for _ in range(probs.shape[0])).format(
                *np.around(pl)))
        return

    def __eq__(self, var2):
        return (self.chrm, self.pos, self.id) == (var2.chrm, var2.pos, var2.id)
    def __ne__(self, var2):
        return (self.chrm, self.pos, self.id) != (var2.chrm, var2.pos, var2.id)
    def __lt__(self, var2):
        return (self.chrm, self.pos, self.id) < (var2.chrm, var2.pos, var2.id)
    def __le__(self, var2):
        return (self.chrm, self.pos, self.id) <= (var2.chrm, var2.pos, var2.id)
    def __gt__(self, var2):
        return (self.chrm, self.pos, self.id) > (var2.chrm, var2.pos, var2.id)
    def __ge__(self, var2):
        return (self.chrm, self.pos, self.id) >= (var2.chrm, var2.pos, var2.id)


class VcfWriter(object):
    """ VCF writer class
    """
    version_options = set(['4.1',])
    def __init__(
            self, filename, mode='w',
            header=('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
                    'INFO', 'FORMAT', SAMPLE_NAME),
            extra_meta_info=FIXED_VCF_MI, version='4.1', ref_fn=None,
            ref_names_and_lens=None, write_vcf_lp=False):
        self.filename = filename
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
            mh.FILE_DATE_MI.format(datetime.date.today().strftime("%Y%m%d")),
            mh.SOURCE_MI.format(MEGALODON_VERSION),
            mh.REF_MI.format(ref_fn)] + contig_mis + extra_meta_info
        if write_vcf_lp:
            self.meta.append(FORMAT_LOG_PROB_MI)

        self.handle = open(self.filename, self.mode, encoding='utf-8')
        self.handle.write('\n'.join('##' + line for line in self.meta) + '\n')
        self.handle.write('#' + '\t'.join(self.header) + '\n')
        return

    def close(self):
        self.handle.close()
        return

    def write_variant(self, variant):
        elements = [getattr(variant, field.lower()) for field in self.header]
        elements = ['.' if e == '' else e for e in elements]
        # VCF POS field is 1-based
        elements[self.header.index('POS')] += 1
        self.handle.write('{}\n'.format('\t'.join(map(str, elements))))
        self.handle.flush()

        return


#################################
##### SNP Aggregation Class #####
#################################

class AggSnps(mh.AbstractAggregationClass):
    """ Class to assist in database queries for per-site aggregation of
    SNP calls over reads.
    """
    def __init__(self, snps_db_fn, write_vcf_log_probs=False):
        # open as read only database
        self.snps_db = sqlite3.connect(snps_db_fn, uri=True)
        self.n_uniq_snps = None
        self.write_vcf_log_probs = write_vcf_log_probs
        return

    def num_uniq(self):
        if self.n_uniq_snps is None:
            self.n_uniq_snps = self.snps_db.execute(
                COUNT_UNIQ_SNPS).fetchone()[0]
        return self.n_uniq_snps

    def iter_uniq(self):
        for q_val in self.snps_db.execute(SEL_UNIQ_SNP_ID):
            yield q_val
        return

    def get_per_read_snp_stats(self, snp_loc):
        return [SNP_DATA(*snp_stats) for snp_stats in self.snps_db.execute(
            SEL_SNP_STATS, snp_loc)]

    def compute_diploid_probs(self, ref_lps, alts_lps, het_factor=1.0):
        def compute_het_lp(a1, a2):
            # order by the inverse log likelihood ratio
            llr_ord = np.argsort(all_lps[a2] - all_lps[a1])
            s_a1_lps = all_lps[a1][llr_ord]
            s_a2_lps = all_lps[a2][llr_ord]
            with np.errstate(divide='ignore'):
                # compute log probability of heterozygous genotype by binomial
                # weighted sum of maximum likelihoods
                het_lp = logsumexp(np.array([
                    np.log(stats.binom.pmf(i, all_lps.shape[1], 0.5)) +
                    np.sum(s_a1_lps[:i]) + np.sum(s_a2_lps[i:])
                    for i in range(all_lps.shape[1] + 1)]))
            return het_lp


        all_lps = np.concatenate([ref_lps.reshape(1, -1), alts_lps], axis=0)
        genotype_lps, het_gts, gts = [], [], []
        # order genotypes as described here under the GL genotype fields section
        # http://samtools.github.io/hts-specs/VCFv4.1.pdf
        for a2 in range(all_lps.shape[0]):
            for a1 in range(a2 + 1):
                gts.append('{}/{}'.format(a1, a2))
                if a1 == a2:
                    genotype_lps.append(np.sum(all_lps[a1]))
                    het_gts.append(False)
                else:
                    genotype_lps.append(compute_het_lp(a1, a2))
                    het_gts.append(True)

        log_prior_weights = np.array([
            0.0 if het_gt else all_lps.shape[1] * np.log(het_factor)
            for het_gt in het_gts])
        log_prior_weights = log_prior_weights - logsumexp(log_prior_weights)
        snp_lps = np.array(genotype_lps) + log_prior_weights
        post_snp_lps = snp_lps - logsumexp(snp_lps)
        return np.exp(post_snp_lps), gts

    def compute_haploid_probs(self, ref_lps, alts_lps):
        snp_lps = np.concatenate([[ref_lps.sum()], alts_lps.sum(axis=1)])
        post_snp_lps = snp_lps - logsumexp(snp_lps)
        return np.exp(post_snp_lps), list(map(str, range(snp_lps.shape[0])))

    def compute_snp_stats(
            self, snp_loc, het_factors, call_mode=DIPLOID_MODE,
            valid_read_ids=None):
        assert call_mode in (HAPLIOD_MODE, DIPLOID_MODE), (
            'Invalid SNP aggregation ploidy call mode: {}.'.format(call_mode))

        pr_snp_stats = self.get_per_read_snp_stats(snp_loc)
        alt_seqs = sorted(set(r_stats.alt_seq for r_stats in pr_snp_stats))
        pr_alt_lps = defaultdict(dict)
        for r_stats in pr_snp_stats:
            if (valid_read_ids is not None and
                r_stats.read_id not in valid_read_ids):
                continue
            pr_alt_lps[r_stats.read_id][r_stats.alt_seq] = r_stats.score
        if len(pr_alt_lps) == 0:
            raise mh.MegaError('No valid reads cover SNP')

        alt_seq_lps = [[] for _ in range(len(alt_seqs))]
        for read_lps in pr_alt_lps.values():
            for i, alt_seq in enumerate(alt_seqs):
                try:
                    alt_seq_lps[i].append(read_lps[alt_seq])
                except KeyError:
                    raise mh.MegaError(
                        'Alternative SNP seqence must exist for all reads.')
        alts_lps = np.stack(alt_seq_lps, axis=0)
        with np.errstate(all='ignore'):
            ref_lps = np.log1p(-np.exp(alts_lps).sum(axis=0))

        r0_stats = pr_snp_stats[0]
        snp_var = Variant(
            chrom=r0_stats.chrm, pos=r0_stats.pos, ref=r0_stats.ref_seq,
            alts=alt_seqs, id=r0_stats.snp_id)
        snp_var.add_tag('DP', '{}'.format(ref_lps.shape[0]))
        snp_var.add_sample_field('DP', '{}'.format(ref_lps.shape[0]))

        if self.write_vcf_log_probs:
            snp_var.add_sample_field('LOG_PROBS', ','.join(
                ';'.join('{:.2f}'.format(lp) for lp in alt_i_lps)
                for alt_i_lps in alts_lps))

        if call_mode == DIPLOID_MODE:
            het_factor = (
                het_factors[0] if len(r0_stats.ref_seq) == 1 and
                len(r0_stats.alt_seq) == 1 else
                het_factors[1])
            diploid_probs, gts = self.compute_diploid_probs(
                ref_lps, alts_lps, het_factor)
            snp_var.add_diploid_probs(diploid_probs, gts)
        elif call_mode == HAPLIOD_MODE:
            haploid_probs, gts = self.compute_haploid_probs(ref_lps, alts_lps)
            snp_var.add_haploid_probs(haploid_probs, gts)

        return snp_var

    def close(self):
        self.snps_db.close()
        return


#########################################
##### Whatshap Mapping Post-process #####
#########################################

def get_whatshap_command(index_variant_fn, whatshap_sort_fn, phase_fn):
    return ('Run following command to obtain phased variants:\n\t\t' +
            'whatshap phase --indels --distrust-genotypes -o {} {} {}').format(
                phase_fn, index_variant_fn, whatshap_sort_fn)

def sort_variants(in_vcf_fn, out_vcf_fn):
    in_vcf_fp = pysam.VariantFile(in_vcf_fn)
    with pysam.VariantFile(
            out_vcf_fn, 'w', header=in_vcf_fp.header) as out_vcf_fp:
        for rec in sorted(in_vcf_fp.fetch(), key=lambda r: (r.chrom, r.start)):
            out_vcf_fp.write(rec)
    return

def index_variants(variant_fn):
    try:
        return pysam.tabix_index(
            variant_fn, force=True, preset='vcf', keep_original=True)
    except OSError:
        # file likely not sorted
        sort_variant_fn = mh.add_fn_suffix(variant_fn, 'sorted')
        sort_variants(variant_fn, sort_variant_fn)
        return pysam.tabix_index(
            sort_variant_fn, force=True, preset='vcf', keep_original=True)


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
