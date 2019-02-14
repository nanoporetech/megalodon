#!/usr/bin/env python3
import os
# set blas library environment variables (without these the cblas calls
# can completely halt processing)
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'


import re
import sys
import queue
import shutil
import argparse
import warnings
import threading
import traceback
from time import sleep
import multiprocessing as mp
from collections import namedtuple, defaultdict

import h5py
import mappy
import pysam
import numpy as np
from tqdm import tqdm
from tqdm._utils import _term_move_up

from megalodon import decode, backends, megalodon_helper as mh


COMP_BASES = dict(zip(map(ord, 'ACGT'), map(ord, 'TGCA')))
DEFAULT_CONTEXT_BASES = 10
DEFAULT_EDGE_BUFFER = 5

MAP_POS = namedtuple('MAP_POS', (
    'chrm', 'strand', 'start', 'end', 'q_trim_start', 'q_trim_end'))

SINGLE_LETTER_CODE = {
    'A':'A', 'C':'C', 'G':'G', 'T':'T', 'B':'[CGT]',
    'D':'[AGT]', 'H':'[ACT]', 'K':'[GT]', 'M':'[AC]',
    'N':'[ACGT]', 'R':'[AG]', 'S':'[CG]', 'V':'[ACG]',
    'W':'[AT]', 'Y':'[CT]'}

BC_NAME = 'basecalls'
BC_OUT_FMTS = ('fasta',)
MAP_NAME = 'mappings'
MAP_OUT_FMTS = ('bam', 'cram', 'sam')
PR_SNP_NAME = 'per_read_snps'
PR_MOD_NAME = 'per_read_mods'
ALIGN_OUTPUTS = set((MAP_NAME, PR_SNP_NAME, PR_MOD_NAME))
OUTPUT_FNS = {
    BC_NAME:'basecalls',
    MAP_NAME:['mappings', 'mappings.summary.txt'],
    PR_SNP_NAME:'per_read_snp_calls.txt',
    PR_MOD_NAME:'per_read_modified_base_calls.txt',
}

_DO_PROFILE = False
_UNEXPECTED_ERROR_CODE = 'Unexpected error'
_UNEXPECTED_ERROR_FN = 'unexpected_snp_calling_errors.{}.err'
_MAX_NUM_UNEXP_ERRORS = 50
_MAX_QUEUE_SIZE = 1000

# mappy aligner with extra attrs
class alignerPlus(mappy.Aligner):
    ref_names_and_lens = None
    out_fmt = None
    ref_fn = None
    pass


#############################
##### Signal Extraction #####
#############################

def med_mad(data, factor=None, axis=None, keepdims=False):
    """Compute the Median Absolute Deviation, i.e., the median
    of the absolute deviations from the median, and the median

    :param data: A :class:`ndarray` object
    :param factor: Factor to scale MAD by. Default (None) is to be consistent
    with the standard deviation of a normal distribution
    (i.e. mad( N(0,\sigma^2) ) = \sigma).
    :param axis: For multidimensional arrays, which axis to calculate over
    :param keepdims: If True, axis is kept as dimension of length 1

    :returns: a tuple containing the median and MAD of the data
    """
    if factor is None:
        factor = 1.4826
    dmed = np.median(data, axis=axis, keepdims=True)
    dmad = factor * np.median(abs(data - dmed), axis=axis, keepdims=True)
    if axis is None:
        dmed = dmed.flatten()[0]
        dmad = dmad.flatten()[0]
    elif not keepdims:
        dmed = dmed.squeeze(axis)
        dmad = dmad.squeeze(axis)
    return dmed, dmad

def extract_read_data(fast5_fn, scale=True):
    """Extract raw signal and read id from single fast5 file.
    :returns: tuple containing numpy array with raw signal  and read unique identifier
    """
    # TODO support mutli-fast5
    try:
        fast5_data = h5py.File(fast5_fn, 'r')
    except:
        raise mh.MegaError('Error opening read file')
    try:
        raw_slot = next(iter(fast5_data['/Raw/Reads'].values()))
    except:
        fast5_data.close()
        raise mh.MegaError('Raw signal not found in /Raw/Reads slot')
    try:
        raw_sig = raw_slot['Signal'][:].astype(np.float32)
    except:
        fast5_data.close()
        raise mh.MegaError('Raw signal not found in Signal dataset')
    read_id = None
    try:
        read_id = raw_slot.attrs.get('read_id')
        read_id = read_id.decode()
    except:
        pass
    fast5_data.close()

    if scale:
        med, mad = med_mad(raw_sig)
        raw_sig = (raw_sig - med) / mad

    return raw_sig, read_id


####################################
##### Read mapping and parsing #####
####################################

def nstate_to_nbase(nstate):
    return int(np.sqrt(0.25 + (0.5 * nstate)) - 0.5)

def revcomp(seq):
    return seq.translate(COMP_BASES)[::-1]

def _map_read_worker(aligner, map_conn, mo_q, add_chr_ref):
    # get mappy aligner thread buffer
    map_thr_buf = mappy.ThreadBuffer()

    while True:
        try:
            r_seq, read_id = map_conn.recv()
        except:
            # exit gracefully
            return
        if r_seq is None:
            break
        try:
            r_algn = next(aligner.map(str(r_seq), buf=map_thr_buf))
        except StopIteration:
            # alignment not produced
            map_conn.send([None, None])
            continue

        ref_seq = aligner.seq(r_algn.ctg, r_algn.r_st, r_algn.r_en)
        if r_algn.strand == -1:
            ref_seq = revcomp(ref_seq)
        chrm = 'chr' + r_algn.ctg if add_chr_ref else r_algn.ctg
        r_algn_data = [
            chrm, r_algn.strand, r_algn.r_st, r_algn.r_en,
            r_algn.q_st, r_algn.q_en, r_algn.cigar]
        map_conn.send([ref_seq, r_algn_data])

        if mo_q is not None:
            mo_q.put((
                read_id, r_seq, r_algn.ctg, r_algn.strand, r_algn.r_st,
                r_algn.q_st, r_algn.q_en, r_algn.cigar))

        del r_seq, r_algn, ref_seq, r_algn_data

    return

def parse_cigar(r_cigar, strand):
    # get each base calls genomic position
    r_to_q_poss = []
    # process cigar ops in read direction
    curr_r_pos, curr_q_pos = 0, 0
    cigar_ops = r_cigar if strand == 1 else r_cigar[::-1]
    for op_len, op in cigar_ops:
        if op == 1:
            # inserted bases into ref
            curr_q_pos += op_len
        elif op in (2, 3):
            # deleted ref bases
            r_to_q_poss.extend(((r_pos, curr_q_pos) for r_pos in range(
                curr_r_pos, curr_r_pos + op_len)))
            curr_r_pos += op_len
        elif op in (0, 7, 8):
            # aligned bases
            r_to_q_poss.extend(zip(
                range(curr_r_pos, curr_r_pos + op_len),
                range(curr_q_pos, curr_q_pos + op_len)))
            curr_q_pos += op_len
            curr_r_pos += op_len
        elif op == 6:
            # padding (shouldn't happen in mappy)
            pass

    return dict(r_to_q_poss)

def map_read(r_seq, read_id, caller_conn):
    """Map read (query) sequence and return:
    1) reference sequence (endcoded as int labels)
    2) mapping from reference to read positions (after trimming)
    3) reference mapping position (including read trimming positions)
    """
    # send seq to _map_read_worker and receive mapped seq and pos
    caller_conn.send((r_seq, read_id))
    r_ref_seq, r_algn = caller_conn.recv()
    if r_ref_seq is None:
        raise mh.MegaError('No alignment')
    chrm, strand, r_st, r_en, q_st, q_en, r_cigar = r_algn

    r_to_q_poss = parse_cigar(r_cigar, strand)
    r_pos = MAP_POS(
        chrm=chrm, strand=strand, start=r_st, end=r_en,
        q_trim_start=q_st, q_trim_end=q_en)

    return r_ref_seq, r_to_q_poss, r_pos


###########################
##### Read processing #####
###########################

def encode_snp_seq(seq):
    """ Return sequence encoded as base 5 converted to integer

    Note base 5 is used so the sequence length and value are encoded
    """
    if seq == '':
        return 0
    return sum((mh.ALPHABET.find(b) + 1) * (5 ** i)
               for i, b in enumerate(seq[::-1]))

def decode_snp_seq(val):
    """ Decode value (encoded via encode_snp_seq) into sequence
    """
    if val == 0:
        return ''
    seq = ''
    for bi in np.base_repr(val, 5):
        if bi == '0':
            raise mh.MegaError('Invalid SNP seq encoding')
        seq += mh.ALPHABET[int(bi) - 1]
    return seq

def simplify_and_encode_snp(snp_ref_seq, snp_alt_seq, ref_pos, max_snp_size):
    """ Simplify SNP when extra bases are included (for indels)
    """
    snp_ref_seq = snp_ref_seq.upper()
    snp_alt_seq = snp_alt_seq.upper()
    # handle cases containing non-canonical base values (e.g. dash for deletion;
    # assume this means full ref or alt deletion)
    if not all(rb in mh.ALPHABET for rb in snp_ref_seq):
        if not all(ab in mh.ALPHABET for ab in snp_alt_seq):
            raise mh.MegaError('Invalid SNP')
        if len(snp_alt_seq) > max_snp_size:
            raise mh.MegaError('SNP too long')
        return 0, encode_snp_seq(snp_alt_seq), ref_pos
    elif not all(ab in mh.ALPHABET for ab in snp_alt_seq):
        if len(snp_ref_seq) > max_snp_size:
            raise mh.MegaError('SNP too long')
        return encode_snp_seq(snp_ref_seq), 0, ref_pos

    # trim base positions that are equal
    while (len(snp_ref_seq) > 0 and len(snp_alt_seq) > 0 and
           snp_ref_seq[0] == snp_alt_seq[0]):
        snp_ref_seq = snp_ref_seq[1:]
        snp_alt_seq = snp_alt_seq[1:]
        ref_pos += 1
    while (len(snp_ref_seq) > 0 and len(snp_alt_seq) > 0 and
           snp_ref_seq[-1] == snp_alt_seq[-1]):
        snp_ref_seq = snp_ref_seq[:-1]
        snp_alt_seq = snp_alt_seq[:-1]
    if len(snp_ref_seq) == 0 and len(snp_alt_seq) == 0:
        raise mh.MegaError('Invalid SNP')

    if np.abs(len(snp_ref_seq) - len(snp_alt_seq)) > max_snp_size:
        raise mh.MegaError('SNP too long')
    return encode_snp_seq(snp_ref_seq), encode_snp_seq(snp_alt_seq), ref_pos

def get_overlapping_snps(r_ref_pos, snps_to_test, edge_buffer):
    """Return SNPs overlapping the read mapped position.

    SNPs within edge buffer of the end of the mapping will be ignored.
    """
    ovlp_snps = []
    try:
        chrm_poss, chrm_ref_es, chrm_alt_es = snps_to_test[r_ref_pos.chrm]
    except KeyError:
        raise mh.MegaError(
            'No SNPs on mapped chromosome/record (see --prepend-chr-*)')
    start_idx, end_idx = np.searchsorted(
        chrm_poss, (r_ref_pos.start + edge_buffer, r_ref_pos.end - edge_buffer))
    if start_idx >= end_idx:
        raise mh.MegaError('No overlapping SNPs')

    for pos, ref_es, alt_es, snp_id in zip(chrm_poss[start_idx:end_idx],
                                           chrm_ref_es[start_idx:end_idx],
                                           chrm_alt_es[start_idx:end_idx],
                                           range(start_idx, end_idx)):
        snp_ref_seq = decode_snp_seq(ref_es)
        snp_alt_seq = decode_snp_seq(alt_es)
        if r_ref_pos.strand == 1:
            read_pos = pos - r_ref_pos.start
        else:
            read_pos = r_ref_pos.end - pos - 1
            snp_ref_seq = revcomp(snp_ref_seq)
            snp_alt_seq = revcomp(snp_alt_seq)
        ovlp_snps.append((read_pos, snp_ref_seq, snp_alt_seq, snp_id))

    return ovlp_snps

############################
##### Sequence Scoring #####
############################

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

def score_seq(tpost, seq, tpost_start=0, tpost_end=None,
              all_paths=False):
    """Score a section of log transition posteriors against a proposed sequence
    using a global mapping.
    :param tpost: `ndarray` containing log transition posteriors to be scored
    :param seq: `ndarray` containing integers encoding proposed sequence
    :param tpost_start: start position within post (Default: 0)
    :param tpost_end: end position within post (Default: full posterior)
    :param all_paths: boolean to produce the forwards all paths score (default Viterbi best path)
    """
    seq = seq.astype(np.uintp)
    if tpost_end is None:
        tpost_end = post.shape[0]

    return decode.score_seq(tpost, seq, tpost_start, tpost_end, all_paths)

def call_read_snps(
        r_ref_pos, snps_to_test, edge_buffer, context_bases, r_ref_seq,
        rl_cumsum, r_to_q_poss, r_post, post_mapped_start, all_paths):
    # call all snps overlapping this read
    r_snp_calls = []
    for r_snp_pos, snp_ref_seq, snp_alt_seq, snp_id in get_overlapping_snps(
            r_ref_pos, snps_to_test, edge_buffer):
        pos_bb, pos_ab = min(context_bases, r_snp_pos), min(
            context_bases, r_ref_seq.shape[0] - r_snp_pos - len(snp_ref_seq))
        pos_ref_seq = r_ref_seq[r_snp_pos - pos_bb:
                                r_snp_pos + pos_ab + len(snp_ref_seq)]
        if any(pos_ref_seq[pos_bb:pos_bb + len(snp_ref_seq)] !=
               np.array([mh.ALPHABET.find(b) for b in snp_ref_seq])):
            raise mh.MegaError(
                'Reference SNP sequence does not match reference FASTA.')
        pos_alt_seq = np.concatenate([
            pos_ref_seq[:pos_bb],
            np.array([mh.ALPHABET.find(b) for b in snp_alt_seq], dtype=np.uintp),
            pos_ref_seq[pos_bb + len(snp_ref_seq):]])
        blk_start, blk_end = (rl_cumsum[r_to_q_poss[r_snp_pos - pos_bb]],
                              rl_cumsum[r_to_q_poss[r_snp_pos + pos_ab]])
        if blk_end - blk_start < (context_bases * 2) + 1:
            # no valid mapping over large inserted query bases
            # i.e. need as many "events/strides" as bases for valid mapping
            continue
        loc_ref_score = score_seq(
            r_post, pos_ref_seq, post_mapped_start + blk_start,
            post_mapped_start + blk_end, all_paths)
        loc_alt_score = score_seq(
            r_post, pos_alt_seq, post_mapped_start + blk_start,
            post_mapped_start + blk_end, all_paths)
        if loc_ref_score is None or loc_alt_score is None:
            raise mh.MegaError('Score computation error (memory error)')

        snp_ref_pos = (r_snp_pos + r_ref_pos.start if r_ref_pos.strand == 1 else
                       r_ref_pos.end - r_snp_pos - len(snp_ref_seq))
        fwd_strand_ref_seq = (snp_ref_seq if r_ref_pos.strand == 1 else
                              revcomp(snp_ref_seq))
        fwd_strand_alt_seq = (snp_alt_seq if r_ref_pos.strand == 1 else
                              revcomp(snp_alt_seq))
        if len(fwd_strand_ref_seq) == 0 or len(fwd_strand_alt_seq) == 0:
            fwd_ref_base = (
                mh.ALPHABET[r_ref_seq[pos_bb - 1]] if r_ref_pos.strand == 1 else
                revcomp(mh.ALPHABET[r_ref_seq[pos_bb + len(snp_ref_seq)]]))
            fwd_strand_ref_seq = fwd_ref_base + fwd_strand_ref_seq
            fwd_strand_alt_seq = fwd_ref_base + fwd_strand_alt_seq
        r_snp_calls.append((
            snp_ref_pos, loc_ref_score - loc_alt_score, fwd_strand_ref_seq,
            fwd_strand_alt_seq, snp_id))

    return r_snp_calls

def rle(x, tol=0):
    """  Run length encoding of array x

    Note: where matching is done with some tolerance, the first element
    of the run is chosen as representative.

    :param x: array
    :param tol: tolerance of match (for continuous arrays)

    :returns: tuple of array containing elements of x and array containing
    length of run
    """

    delta_x = np.ediff1d(x, to_begin=1)
    starts = np.where(np.absolute(delta_x) > tol)[0]
    last_runlength = len(x) - starts[-1]
    runlength = np.ediff1d(starts, to_end=last_runlength)

    return x[starts], runlength

def decode_post(r_post, collapse_alphabet=mh.ALPHABET):
    """Decode a posterior using Viterbi algorithm for transducer.
    :param r_post: numpy array containing transducer posteriors.
    :param collapse_alphabet: alphabet corresponding to flip-flop labels.
    :returns: tuple containing (base calls, score and raw block positions).
    """
    nblock, nstate = r_post.shape[:2]
    nbase = len(set(collapse_alphabet))
    if nbase != nstate_to_nbase(nstate):
        raise mh.MegaError(
            'Incompatible decoding alphabet and posterior states.')

    path = np.zeros(nblock + 1, dtype=np.uintp)
    qpath = np.zeros(nblock + 1, dtype=np.float32)

    score = decode.crf_flipflop_viterbi(r_post, path, qpath)

    runval, runlen = rle(path)
    basecall = ''.join(collapse_alphabet[int(b) % nbase] for b in runval)

    return basecall, score, runlen

def process_read(
        raw_sig, read_id, model_info, bc_q, caller_conn, snps_to_test,
        snp_all_paths, snps_q, mods_q, alphabet_info,
        context_bases=DEFAULT_CONTEXT_BASES, edge_buffer=DEFAULT_EDGE_BUFFER):
    if model_info.is_cat_mod:
        bc_weights, mod_weights = model_info.run_model(
            raw_sig, alphabet_info.n_can_state)
    else:
        bc_weights = model_info.run_model(raw_sig)

    r_post = decode.crf_flipflop_trans_post(bc_weights, log=True)
    r_seq, score, runlen = decode_post(r_post, alphabet_info.collapse_alphabet)
    if bc_q is not None:
        bc_q.put((read_id, r_seq))

    # if no mapping connection return after basecalls are passed out
    if caller_conn is None: return

    # map read and record mapping from reference to query positions
    r_ref_seq, r_to_q_poss, r_ref_pos = map_read(r_seq, read_id, caller_conn)
    np_ref_seq = np.array([
        mh.ALPHABET.find(b) for b in r_ref_seq], dtype=np.uintp)

    # get mapped start in post and run len to mapped bit of output
    post_mapped_start = sum(runlen[:r_ref_pos.q_trim_start])
    rl_cumsum = np.cumsum(np.concatenate([
        [0], runlen[r_ref_pos.q_trim_start:r_ref_pos.q_trim_end]]))

    if snps_q is not None:
        try:
            snps_q.put((
                call_read_snps(
                    r_ref_pos, snps_to_test, edge_buffer, context_bases,
                    np_ref_seq, rl_cumsum, r_to_q_poss, r_post,
                    post_mapped_start, snp_all_paths),
                (read_id, r_ref_pos.chrm, r_ref_pos.strand)))
        except mh.MegaError:
            pass
    if mods_q is not None:
        r_post_w_mods = np.concatenate([r_post, mod_weights], axis=1)
        try:
            mods_q.put((
                call_read_mods(
                    r_ref_pos, edge_buffer, context_bases, r_ref_seq,
                    np_ref_seq, rl_cumsum, r_to_q_poss, r_post_w_mods,
                    post_mapped_start, alphabet_info),
                (read_id, r_ref_pos.chrm, r_ref_pos.strand)))
        except mh.MegaError:
            pass

    return


############################
##### Multi-processing #####
############################

def _get_mods_queue(mods_q, mods_conn, mods_fn):
    mods_fp = open(mods_fn, 'w')

    while True:
        try:
            # note strand is +1 for fwd or -1 for rev
            r_mod_calls, (read_id, chrm, strand) = mods_q.get(block=False)
            # TODO record stats for VFC output
            # TODO write sqlite output
            # would involve batching and creating several conversion tables
            # for var strings (read_if and chrms).
            mods_fp.write('\n'.join((
                '{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                    read_id, chrm, strand, pos, score, raw_motif, mod_base)
                for pos, score, raw_motif, mod_base in r_mod_calls)) + '\n')
            mods_fp.flush()
        except queue.Empty:
            if mods_conn.poll():
                break
            sleep(0.1)
            continue

    while not mods_q.empty():
        r_mod_calls, (read_id, chrm, strand) = mods_q.get(block=False)
        mods_fp.write('\n'.join((
            '{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                read_id, chrm, strand, pos, score, raw_motif, mod_base)
            for pos, score, raw_motif, mod_base in r_mod_calls)) + '\n')
        mods_fp.flush()
    mods_fp.close()

    return

def _get_snps_queue(snps_q, snps_conn, snp_id_tbl, snps_fn):
    snps_fp = open(snps_fn, 'w')

    while True:
        try:
            # note strand is +1 for fwd or -1 for rev
            r_snp_calls, (read_id, chrm, strand) = snps_q.get(block=False)
            # TODO record stats for VFC output
            # TODO write sqlite output
            # would involve batching and creating several conversion tables
            # for var strings (read_if and chrms).
            snps_fp.write('\n'.join((
                '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                    read_id, chrm, strand, pos, score, ref_seq, alt_seq,
                    snp_id_tbl[chrm][snp_i])
                for pos, score, ref_seq, alt_seq, snp_i in r_snp_calls)) + '\n')
            snps_fp.flush()
        except queue.Empty:
            if snps_conn.poll():
                break
            sleep(0.1)
            continue

    while not snps_q.empty():
        r_snp_calls, (read_id, chrm, strand) = snps_q.get(block=False)
        snps_fp.write('\n'.join((
            '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                read_id, chrm, strand, pos, score, ref_seq, alt_seq,
                snp_id_tbl[chrm][snp_i])
            for pos, score, ref_seq, alt_seq, snp_i in r_snp_calls)) + '\n')
        snps_fp.flush()
    snps_fp.close()

    return

def _get_map_queue(
        mo_q, map_conn, out_dir, ref_names_and_lens, map_fmt, ref_fn):
    def write_alignment(
            read_id, r_seq, chrm, strand, r_st, q_st, q_en, r_cigar):
        r_seq = r_seq[q_st:q_en]
        if strand == -1:
            r_cigar = r_cigar[::-1]

        a = pysam.AlignedSegment()
        a.query_name = read_id
        a.query_sequence = r_seq if strand == 1 else revcomp(r_seq)
        a.flag = 0 if strand == 1 else 16
        a.reference_id = map_fp.get_tid(chrm)
        a.reference_start = r_st
        a.cigartuples = [(op, op_l) for op_l, op in r_cigar]
        a.template_length = q_en - q_st
        map_fp.write(a)

        nalign, nmatch, ndel, nins = [0,] * 4
        for op_len, op in r_cigar:
            if op not in (4, 5): nalign += op_len
            if op in (0, 7): nmatch += op_len
            elif op in (2, 3): ndel += op_len
            elif op == 1: nins += op_len
        # compute alignment stats
        summ_fp.write('{}\t{:.2f}\t{}\t{}\t{}\t{}\n'.format(
            read_id, 100 * nmatch / float(nalign), nalign, nmatch, ndel, nins))
        summ_fp.flush()

        return


    map_bn, summ_fn = OUTPUT_FNS[MAP_NAME]

    summ_fp = open(os.path.join(out_dir, summ_fn), 'w')
    summ_fp.write('read_id\tpct_identity\tnum_align\tnum_match\t' +
                  'num_del\tnum_ins\n')

    map_fn = os.path.join(out_dir, map_bn + '.' + map_fmt)
    if map_fmt == 'bam': w_mode = 'wb'
    elif map_fmt == 'cram': w_mode = 'wc'
    elif map_fmt == 'sam': w_mode = 'w'
    else:
        raise mh.MegaError('Invalid mapping output format\n')
    map_fp = pysam.AlignmentFile(
        map_fn, w_mode, reference_names=ref_names_and_lens[0],
        reference_lengths=ref_names_and_lens[1], reference_filename=ref_fn)

    try:
        while True:
            try:
                write_alignment(*mo_q.get(block=False))
            except queue.Empty:
                if map_conn.poll():
                    break
                sleep(0.1)
                continue

        while not mo_q.empty():
            write_alignment(*mo_q.get(block=False))
    finally:
        map_fp.close()
        summ_fp.close()

    return

def _get_bc_queue(bc_q, bc_conn, out_dir, bc_fmt):
    bc_fp = open(os.path.join(out_dir, OUTPUT_FNS[BC_NAME] + '.' + bc_fmt), 'w')

    while True:
        try:
            # TODO add quality output to add fastq option
            read_id, r_seq = bc_q.get(block=False)
            bc_fp.write('>{}\n{}\n'.format(read_id, r_seq))
            bc_fp.flush()
        except queue.Empty:
            if bc_conn.poll():
                break
            sleep(0.1)
            continue

    while not bc_q.empty():
        read_id, r_seq = bc_q.get(block=False)
        bc_fp.write('>{}\n{}\n'.format(read_id, r_seq))
        bc_fp.flush()
    bc_fp.close()

    return

def get_read_files(fast5s_dir):
    """Get all fast5 files recursively below this directory
    """
    fast5_fns = []
    # walk through directory structure searching for fast5 files
    for root, _, fns in os.walk(fast5s_dir):
        for fn in fns:
            if not fn.endswith('.fast5'): continue
            fast5_fns.append(os.path.join(root, fn))

    return fast5_fns

def _fill_files_queue(fast5_q, fast5_fns, num_ts):
    for fast5_fn in fast5_fns:
        fast5_q.put(fast5_fn)
    for _ in range(num_ts):
        fast5_q.put(None)

    return

def _process_reads_worker(
        fast5_q, bc_q, mo_q, snps_q, failed_reads_q, model_info, caller_conn,
        snps_to_test, snp_all_paths, alphabet_info, mods_q):
    model_info.prep_model_worker()

    while True:
        try:
            fast5_fn = fast5_q.get(block=False)
        except queue.Empty:
            sleep(0.1)
            continue

        if fast5_fn is None:
            if caller_conn is not None:
                caller_conn.send(None)
            break

        try:
            raw_sig, read_id = extract_read_data(fast5_fn)
            process_read(
                raw_sig, read_id, model_info, bc_q, caller_conn, snps_to_test,
                snp_all_paths, snps_q, mods_q, alphabet_info)
            failed_reads_q.put((False, None, None, None))
        except KeyboardInterrupt:
            failed_reads_q.put((True, 'Keyboard interrupt', fast5_fn, None))
            return
        except mh.MegaError as e:
            failed_reads_q.put((True, str(e), fast5_fn, None))
        except:
            failed_reads_q.put((
                True, _UNEXPECTED_ERROR_CODE, fast5_fn, traceback.format_exc()))

    return

if _DO_PROFILE:
    _process_reads_wrapper = _process_reads_worker
    def _process_reads_worker(*args):
        import cProfile
        cProfile.runctx('_process_reads_wrapper(*args)', globals(), locals(),
                        filename='read_processing.prof')
        return


##################################
##### Dynamic error updating #####
##################################

def format_fail_summ(header, fail_summ=[], num_reads=0, num_errs=None):
    summ_errs = sorted(fail_summ)[::-1]
    if num_errs is not None:
        summ_errs = summ_errs[:num_errs]
        if len(summ_errs) < num_errs:
            summ_errs.extend([(None, '') for _ in range(
                num_errs - len(summ_errs))])
    errs_str = '\n'.join("{:8.1f}% ({:>7} reads)".format(
        100 * n_fns / float(num_reads), n_fns) + " : " + '{:<80}'.format(err)
        if (n_fns is not None and num_reads > 0) else
        '     -----' for n_fns, err in summ_errs)
    return '\n'.join((header, errs_str))

def prep_errors_bar(num_update_errors, tot_reads, suppress_progress):
    if num_update_errors > 0:
        # add lines for dynamic error messages
        sys.stderr.write(
            '\n'.join(['' for _ in range(num_update_errors + 2)]))
    bar, prog_prefix, bar_header = None, None, None
    if suppress_progress:
        num_update_errors = 0
    else:
        bar = tqdm(total=tot_reads, smoothing=0)
    if num_update_errors > 0:
        prog_prefix = ''.join(
            [_term_move_up(),] * (num_update_errors + 1)) + '\r'
        bar_header = (
            str(num_update_errors) + ' most common unsuccessful read types:')
        # write failed read update header
        bar.write(prog_prefix + format_fail_summ(
            bar_header, num_errs=num_update_errors), file=sys.stderr)

    return bar, prog_prefix, bar_header

def _get_fail_queue(
        failed_reads_q, f_conn, tot_reads, num_update_errors,
        suppress_progress):
    reads_called = 0
    failed_reads = defaultdict(list)
    bar, prog_prefix, bar_header = prep_errors_bar(
        num_update_errors, tot_reads, suppress_progress)
    while True:
        try:
            is_err, err_type, fast5_fn, err_tb = failed_reads_q.get(block=False)
            if is_err:
                failed_reads[err_type].append(fast5_fn)
                if err_type == _UNEXPECTED_ERROR_CODE:
                    if len(failed_reads[err_type]) == 1:
                        unexp_err_fp = open(_UNEXPECTED_ERROR_FN.format(
                            np.random.randint(10000)), 'w')
                    if len(failed_reads[err_type]) >= _MAX_NUM_UNEXP_ERRORS:
                        unexp_err_fp.close()
                    else:
                        unexp_err_fp.write(
                            fast5_fn + '\n:::\n' + err_tb + '\n\n\n')
                        unexp_err_fp.flush()

            if not suppress_progress: bar.update(1)
            reads_called += 1
            if num_update_errors > 0:
                bar.write(prog_prefix + format_fail_summ(
                    bar_header,
                    [(len(fns), err) for err, fns in failed_reads.items()],
                    reads_called, num_update_errors), file=sys.stderr)
        except queue.Empty:
            # check if all reads are done signal was sent from main thread
            if f_conn.poll():
                break
            try:
                sleep(0.1)
            except KeyboardInterrupt:
                # exit gracefully on keyboard inturrupt
                return
            continue

    if not suppress_progress: bar.close()
    if len(failed_reads[_UNEXPECTED_ERROR_CODE]) >= 1:
        sys.stderr.write((
            '******* WARNING *******\n\tUnexpected errors occured. See full ' +
            'error stack traces for first (up to) {0:d} errors in ' +
            '"{1}"\n').format(_MAX_NUM_UNEXP_ERRORS, unexp_err_fp.name))

    return


###############################
##### All read processing #####
###############################

def process_all_reads(
        fast5s_dir, num_reads, model_info, outputs, out_dir, bc_fmt, aligner,
        add_chr_ref, snps_data, num_ts, num_update_errors, suppress_progress,
        alphabet_info):
    def create_getter_q(getter_func, args):
        q = mp.Queue(maxsize=_MAX_QUEUE_SIZE)
        main_conn, conn = mp.Pipe()
        p = mp.Process(target=getter_func, daemon=True, args=(q, conn, *args))
        p.start()
        return q, p, main_conn


    sys.stderr.write('Searching for reads\n')
    fast5_fns = get_read_files(fast5s_dir)
    if num_reads is not None:
        fast5_fns = fast5_fns[:num_reads]

    sys.stderr.write('Preparing workers and calling reads\n')
    # read filename queue filler
    fast5_q = mp.Queue(maxsize=_MAX_QUEUE_SIZE)
    files_p = mp.Process(
        target=_fill_files_queue, args=(fast5_q, fast5_fns, num_ts),
        daemon=True)
    files_p.start()
    # progress and failed reads getter
    failed_reads_q, f_p, main_f_conn = create_getter_q(
            _get_fail_queue, (len(fast5_fns), num_update_errors,
                              suppress_progress))

    # start output type getters/writers
    (bc_q, bc_p, main_bc_conn, mo_q, mo_p, main_mo_conn, snps_q, snps_p,
     main_snps_conn, mods_q, mods_p, main_mods_conn) = [None,] * 12
    if BC_NAME in outputs:
        bc_q, bc_p, main_bc_conn = create_getter_q(
            _get_bc_queue, (out_dir, bc_fmt))
    if MAP_NAME in outputs:
        mo_q, mo_p, main_mo_conn = create_getter_q(
            _get_map_queue, (out_dir, aligner.ref_names_and_lens,
                             aligner.out_fmt, aligner.ref_fn))
    if PR_SNP_NAME in outputs:
        snps_q, snps_p, main_snps_conn = create_getter_q(
            _get_snps_queue, (snps_data.snp_id_tbl,
                              os.path.join(out_dir, OUTPUT_FNS[PR_SNP_NAME])))
    if PR_MOD_NAME in outputs:
        mods_q, mods_p, main_mods_conn = create_getter_q(
            _get_mods_queue, (os.path.join(out_dir, OUTPUT_FNS[PR_MOD_NAME]),))

    call_snp_ps, map_conns = [], []
    for _ in range(num_ts):
        if aligner is None:
            map_conn, caller_conn = None, None
        else:
            map_conn, caller_conn = mp.Pipe()
        map_conns.append(map_conn)
        p = mp.Process(
            target=_process_reads_worker, args=(
                fast5_q, bc_q, mo_q, snps_q, failed_reads_q, model_info,
                caller_conn, snps_data.snps_to_test, snps_data.all_paths,
                alphabet_info, mods_q))
        p.daemon = True
        p.start()
        call_snp_ps.append(p)
    sleep(0.1)

    # perform mapping in threads for mappy shared memory interface
    if aligner is None:
        map_read_ts = None
    else:
        map_read_ts = []
        for map_conn in map_conns:
            t = threading.Thread(
                target=_map_read_worker,
                args=(aligner, map_conn, mo_q, add_chr_ref))
            t.daemon = True
            t.start()
            map_read_ts.append(t)

    try:
        files_p.join()
    except KeyboardInterrupt:
        sys.stderr.write('Exiting due to keyboard interrupt.\n')
        sys.exit(1)
    for call_snp_p in call_snp_ps:
        call_snp_p.join()
    if map_read_ts is not None:
        for map_t in map_read_ts:
            map_t.join()
    # comm to getter processes to return
    if f_p.is_alive():
        main_f_conn.send(True)
        f_p.join()
    for on, p, main_conn in (
            (BC_NAME, bc_p, main_bc_conn),
            (MAP_NAME, mo_p, main_mo_conn),
            (PR_SNP_NAME, snps_p, main_snps_conn),
            (PR_MOD_NAME, mods_p, main_mods_conn)):
        if on in outputs and p.is_alive():
            main_conn.send(True)
            p.join()

    return


################################
########## Parse SNPs ##########
################################

class Snps(object):
    def __init__(self, snp_fn, do_prepend_chr_vcf, max_snp_size, all_paths):
        self.all_paths = all_paths
        if snp_fn is None:
            self.snps_to_test = None
            self.snp_id_tbl = None
            return

        sys.stderr.write('Loading SNPs.\n')
        raw_snps_to_test = defaultdict(lambda: defaultdict(list))
        warned_invalid_line = False
        n_skipped_snps = 0
        with open(snp_fn) as fp:
            for line in fp:
                if line.startswith('#'): continue
                try:
                    chrm, pos, snp_id, ref_seq, alt_seq = line.split()[:5]
                except:
                    if not warned_invalid_line:
                        sys.stderr.write(
                            'WARNING: Encountered invalid VCF line. Silently ' +
                            'ignoring any further invalid lines.\n\t' + line)
                    warned_invalid_line = True

                try:
                    ref_es, alt_es, pos = simplify_and_encode_snp(
                        ref_seq, alt_seq, int(pos), max_snp_size)
                except mh.MegaError:
                    n_skipped_snps += 1
                    continue
                raw_snps_to_test[chrm][(pos, ref_es, alt_es)].append(snp_id)

        # re-organize parsed data
        self.snps_to_test = {}
        self.snp_id_tbl = {}
        for chrm, chrm_snps in raw_snps_to_test.items():
            if do_prepend_chr_vcf:
                chrm = 'chr' + chrm
            # note conversion to 0-based coordinates
            s_poss, s_ref_es, s_alt_es, s_snp_ids = zip(*sorted(
                (pos - 1, ref_es, alt_es, ';'.join(snp_ids))
                for (pos, ref_es, alt_es), snp_ids in chrm_snps.items()))
            # conversion table to use in stats process
            self.snp_id_tbl[chrm] = s_snp_ids
            # note lock=False makes the non-safe, but they are read-only
            s_poss = mp.Array('i', s_poss, lock=False)
            s_ref_es = mp.Array('I', s_ref_es, lock=False)
            s_alt_es = mp.Array('I', s_alt_es, lock=False)
            # numpy array to search sorted more efficiently
            self.snps_to_test[chrm] = (s_poss, s_ref_es, s_alt_es)

        n_uniq_snps = sum(len(cs_snps) for cs_snps in raw_snps_to_test.values())
        sys.stderr.write((
            ('Loaded {} SNPs. (Skipped {} entries due to incompatible ' +
             'SNP type)\n')).format(n_uniq_snps, n_skipped_snps))

        return

class AlphabetInfo(object):
    def _parse_cat_mods(self):
        self.n_can_state = (self.ncan_base + self.ncan_base) * (
            self.ncan_base + 1)
        self.nmod_base = len(set(self.alphabet)) - self.ncan_base
        # record the canonical base group indices for modified base
        # categorical output layer
        self.can_mods_offsets = np.cumsum([0] + [
            self.collapse_alphabet.count(can_b)
            for can_b in self.alphabet[:self.ncan_base]]).astype(np.uintp)

        # create table of string mod labels to ints taken by
        # categorical modifications model
        self.str_to_int_mod_labels = {}
        can_grouped_mods = defaultdict(int)
        for mod_lab, can_lab in zip(
                self.alphabet[self.ncan_base:],
                self.collapse_alphabet[self.ncan_base:]):
            can_grouped_mods[can_lab] += 1
            self.str_to_int_mod_labels[mod_lab] = can_grouped_mods[can_lab]

        return

    def _parse_mod_motifs(self, all_mod_motifs_raw):
        # note only works for mod_refactor models currently
        self.all_mod_motifs = []
        if all_mod_motifs_raw is None:
            for can_base, mod_base in zip(
                    self.collapse_alphabet[self.ncan_base:],
                    self.alphabet[self.ncan_base:]):
                self.all_mod_motifs.append((
                    re.compile(can_base), 0, mod_base, can_base))
        else:
            # parse detection motifs
            for mod_motifs_raw in all_mod_motifs_raw:
                mod_base, motifs_raw = mod_motifs_raw.split(':')
                assert mod_base in self.alphabet[self.ncan_base:], (
                    'Modified base label ({}) not found in model ' +
                    'alphabet ({}).').format(mod_base, self.alphabet)
                for mod_motifs_raw in motifs_raw.split(','):
                    raw_motif, pos = mod_motifs_raw.split('-')
                    motif = re.compile(''.join(
                        SINGLE_LETTER_CODE[letter] for letter in raw_motif))
                    self.all_mod_motifs.append(
                        (motif, int(pos), mod_base, raw_motif))

        return

    def __init__(
            self, model_info, all_mod_motifs_raw, mod_all_paths,
            override_alphabet):
        self.alphabet = model_info.alphabet
        self.collapse_alphabet = model_info.collapse_alphabet
        self.ncan_base = len(set(self.collapse_alphabet))
        if override_alphabet is not None:
            self.alphabet, self.collapse_alphabet = override_alphabet
        try:
            self.alphabet = self.alphabet.decode()
            self.collapse_alphabet = self.collapse_alphabet.decode()
        except:
            pass
        sys.stderr.write('Using alphabet {} collapsed to {}\n'.format(
            self.alphabet, self.collapse_alphabet))

        self.nbase = len(self.alphabet)
        if model_info.is_cat_mod:
            self._parse_cat_mods()
            assert (
                model_info.output_size - self.n_can_state ==
                self.nmod_base + 1), (
                    'Alphabet ({}) and model number of modified bases ({}) ' +
                    'do not agree. See --override-alphabet.'
                ).format(self.alphabet,
                         model_info.output_size - self.n_can_state - 1)

        # parse mod motifs or use "swap" base if no motif provided
        self._parse_mod_motifs(all_mod_motifs_raw)
        self.mod_all_paths = mod_all_paths

        return

def mkdir(out_dir, overwrite):
    if os.path.exists(out_dir):
        if not overwrite:
            sys.stderr.write(
                '*' * 100 + '\nERROR: --output-directory exists and ' +
                '--overwrite is not set.\n' + '*' * 100 + '\n')
            sys.exit(1)
        if os.path.isfile(out_dir) or os.path.islink(out_dir):
            os.remove(out_dir)
        else:
            shutil.rmtree(out_dir)
    os.mkdir(out_dir)

    return


##########################
########## Main ##########
##########################

##############################
###### Argument Parsing ######
##############################

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'fast5s_dir',
        help='Directory containing raw fast5 (will be searched recursively).')

    mdl_grp = parser.add_argument_group('Model Arguments')
    mdl_grp.add_argument(
        '--taiyaki-model-filename',
        help='Taiyaki model checkpoint file.')
    mdl_grp.add_argument(
        '--flappie-model-name',
        help='Flappie model name.')

    out_grp = parser.add_argument_group('Output Arguments')
    out_grp.add_argument(
        '--outputs', nargs='+',
        default=['basecalls',], choices=tuple(OUTPUT_FNS.keys()),
        help='Output type(s) to produce. Default: %(default)s')
    out_grp.add_argument(
        '--output-directory',
        default='megalodon_results',
        help='Directory to store output results. Default: %(default)s')
    out_grp.add_argument(
        '--overwrite', action='store_true',
        help='Overwrite output directory if it exists.')
    out_grp.add_argument(
        '--num-reads', type=int,
        help='Number of reads to process. Default: All reads')

    bc_grp = parser.add_argument_group('Basecall Arguments')
    bc_grp.add_argument(
        '--basecalls-format', choices=BC_OUT_FMTS, default=BC_OUT_FMTS[0],
        help='Basecalls output format. Choices: {}'.format(
            ', '.join(BC_OUT_FMTS)))

    map_grp = parser.add_argument_group('Mapping Arguments')
    map_grp.add_argument(
        '--reference',
        help='Reference FASTA file used for mapping called reads.')
    map_grp.add_argument(
        '--mappings-format', choices=MAP_OUT_FMTS, default=MAP_OUT_FMTS[0],
        help='Mappings output format. Choices: {}'.format(
            ', '.join(MAP_OUT_FMTS)))

    snp_grp = parser.add_argument_group('SNP Arguments')
    snp_grp.add_argument(
        '--snp-filename',
        help='SNPs to call for each read in VCF format (required for output).')
    snp_grp.add_argument(
        '--prepend-chr-vcf', action='store_true',
        help='Prepend "chr" to chromosome names from VCF to match ' +
        'reference names.')
    snp_grp.add_argument(
        '--max-snp-size', type=int, default=5,
        help='Maximum difference in number of reference and alternate bases. ' +
        'Default: %(default)d')
    snp_grp.add_argument(
        '--snp-all-paths', action='store_true',
        help='Compute forwards algorithm all paths score. (Default: Viterbi ' +
        'best-path score)')

    mod_grp = parser.add_argument_group('Modified Base Arguments')
    mod_grp.add_argument(
        '--mod-motifs', nargs='+',
        help='Restrict modified base calls to specified motifs. Format as ' +
        '"[mod_base]:[motif]-[relative_pos]". For CpG, dcm and dam calling ' +
        'use "Z:CG-0 Z:CCWGG-1 Y:GATC-1". Default call all valid positions.')
    mod_grp.add_argument(
        '--mod-all-paths', action='store_true',
        help='Compute forwards algorithm all paths score for modified base ' +
        'calls. (Default: Viterbi best-path score)')

    misc_grp = parser.add_argument_group('Miscellaneous Arguments')
    misc_grp.add_argument(
        '--override-alphabet', nargs=2,
        help='Override alphabet and collapse alphabet from model.')
    misc_grp.add_argument(
        '--prepend-chr-ref', action='store_true',
        help='Prepend "chr" to chromosome names from reference to match ' +
        'VCF names.')
    misc_grp.add_argument(
        '--suppress-progress', action='store_true',
        help='Suppress progress bar output.')
    misc_grp.add_argument(
        '--verbose-read-progress', type=int, default=0,
        help='Output verbose output on read progress. Outputs N most ' +
        'common points where reads could not be processed further. ' +
        'Default: %(default)d')
    misc_grp.add_argument(
        '--processes', type=int, default=1,
        help='Number of parallel processes. Default: %(default)d')
    misc_grp.add_argument(
        '--device', default=None, type=int,
        help='CUDA device to use (only valid for taiyaki), or None to use CPU')

    return parser


def _main():
    args = get_parser().parse_args()
    if _DO_PROFILE:
        if args.processes > 1:
            msg = ('Running profiling with multiple processes is ' +
                   'not allowed. Setting to single process.')
            args.processes = 1
        else:
            msg = 'Running profiling. This may slow processing.'
        sys.stderr.write(
            '*' * 100 + '\nWARNING: ' + msg + '\n' + '*' * 100 + '\n')

    mkdir(args.output_directory, args.overwrite)
    model_info = backends.ModelInfo(
        args.flappie_model_name, args.taiyaki_model_filename, args.device)

    alphabet_info = AlphabetInfo(
        model_info, args.mod_motifs, args.mod_all_paths, args.override_alphabet)
    if model_info.is_cat_mod and PR_MOD_NAME not in args.outputs:
        sys.stderr.write(
            '*' * 100 + '\nWARNING: Categorical modifications model ' +
            'provided, but {} not requested '.format(PR_MOD_NAME) +
            '(via --outputs). Modified base output will not be produced.\n' +
            '*' * 100 + '\n')
    if args.mod_motifs is not None and PR_MOD_NAME not in args.outputs:
        sys.stderr.write(
            '*' * 100 + '\nWARNING: --mod-motifs provided, but ' +
            '{} not requested '.format(PR_MOD_NAME) +
            '(via --outputs). Argument will be ignored.\n' + '*' * 100 + '\n')

    if PR_SNP_NAME in args.outputs and args.snp_filename is None:
        sys.stderr.write(
            '*' * 100 + '\nERROR: {} output requested, '.format(PR_SNP_NAME) +
            'but --snp-filename provided.\n' + '*' * 100 + '\n')
        sys.exit(1)
    if PR_SNP_NAME in args.outputs and not (
            model_info.is_cat_mod or
            nstate_to_nbase(model_info.output_size) == 4):
        sys.stderr.write(
            '*' * 100 + '\nERROR: SNP calling from standard modified base ' +
            'flip-flop model is not supported.\n' + '*' * 100 + '\n')
        sys.exit(1)
    # snps data object loads with None snp_fn for easier handling downstream
    snps_data = Snps(
        args.snp_filename, args.prepend_chr_vcf, args.max_snp_size,
        args.snp_all_paths)
    if args.snp_filename is not None and PR_SNP_NAME not in args.outputs:
        sys.stderr.write(
            '*' * 100 + '\nWARNING: --snps-filename provided, but ' +
            'per_read_snp not requested (via --outputs). Argument will be ' +
            'ignored.\n' + '*' * 100 + '\n')

    do_align = len(ALIGN_OUTPUTS.intersection(args.outputs)) > 0
    if do_align:
        if args.reference is None:
            sys.stderr.write(
                '*' * 100 + '\nERROR: Output(s) requiring reference ' +
                'alignment requested, but --reference not provided.\n' +
                '*' * 100 + '\n')
            sys.exit(1)
        sys.stderr.write('Loading reference.\n')
        aligner = alignerPlus(
            str(args.reference), preset=str('map-ont'), best_n=1)
        setattr(aligner, 'out_fmt', args.mappings_format)
        setattr(aligner, 'ref_fn', args.reference)
        if MAP_NAME in args.outputs:
            # extract reference names and lengths
            with pysam.FastaFile(args.reference) as ref:
                setattr(aligner, 'ref_names_and_lens',
                        (ref.references, ref.lengths))
    else:
        aligner = None
    if args.reference is not None and not do_align:
        sys.stderr.write(
            '*' * 100 + '\nWARNING: --reference provided, but no --outputs ' +
            'requiring alignment was requested. Argument will be ignored.\n' +
            '*' * 100 + '\n')

    process_all_reads(
        args.fast5s_dir, args.num_reads, model_info, args.outputs,
        args.output_directory, args.basecalls_format, aligner,
        args.prepend_chr_ref, snps_data, args.processes,
        args.verbose_read_progress, args.suppress_progress,
        alphabet_info)

    return

if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
