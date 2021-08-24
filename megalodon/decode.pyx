# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False

# # cython: profile=True

import cython
from libc.stdint cimport int64_t, uintptr_t
from libc.stdlib cimport calloc, free
from libc.math cimport sqrt, HUGE_VALF

import numpy as np
from megalodon.megalodon_helper import MegaError
from megalodon import logging

LOGGER = logging.get_logger()

cdef extern from "math.h":
    float expf(float x)
    float log1pf(float x)
    float fabsf(float x)
    float fmaxf(float x, float y)


ALPHABET = 'ACGT'


####################
# Helper functions #
####################

"""
The overhead from (even inline) helper function calls mean they are
not used, but left here for reference.
"""

cdef inline size_t nstate_to_nbase(size_t nstate):
    """Convert number of flip-flop states to number of bases and check for
    valid states number.
    """
    cdef float nbase_f = sqrt(0.25 + (0.5 * nstate)) - 0.5
    cdef size_t nbase = <size_t>nbase_f
    return nbase

cdef inline size_t trans_lookup(size_t from_b, size_t to_b, size_t nbase):
    if to_b < nbase:
        return to_b * (nbase + nbase) + from_b
    return nbase * (nbase + nbase) + from_b

cdef inline float logaddexp(float x, float y):
    return fmaxf(x, y) + log1pf(expf(-fabsf(x - y)))

cdef inline void log_row_normalise_inplace(
        float * mat, size_t ncol, size_t nrow):
    cdef size_t col, row, offset
    cdef float row_logsum
    for col in range(ncol):
        offset = col * nrow
        row_logsum = mat[offset]
        for row in range(1, nrow):
            row_logsum = logaddexp(row_logsum, mat[offset + row])

        for row in range(0, nrow):
            mat[offset + row] -= row_logsum

cdef inline void get_stay_step(
        size_t * seq, size_t nseq, size_t nbase,
        size_t * stay_indices, size_t * step_indices):
    cdef size_t * flop_mask_states = <size_t *> calloc(nseq, sizeof(size_t))
    stay_indices[0] = trans_lookup(seq[0], seq[0], nbase)
    flop_mask_states[0] = seq[0]
    cdef size_t seq_pos, prev_fm_state, curr_fm_state
    for seq_pos in range(1, nseq):
        prev_fm_state = flop_mask_states[seq_pos - 1]
        if seq[seq_pos] == prev_fm_state:
            curr_fm_state = seq[seq_pos] + nbase
        else:
            curr_fm_state = seq[seq_pos]
        flop_mask_states[seq_pos] = curr_fm_state
        stay_indices[seq_pos] = trans_lookup(
            curr_fm_state, curr_fm_state, nbase)
        step_indices[seq_pos - 1] = trans_lookup(
            prev_fm_state, curr_fm_state, nbase)
    free(flop_mask_states)


######################
# Decoding functions #
######################

cdef inline void decode_forward_step(
    float * curr_logprob,
    size_t nbase,
    float * prev_fwd,
    float * curr_fwd,
    size_t * tb
):
    cdef size_t nff_state = nbase + nbase

    # to flip base states
    cdef size_t offset_to_state, to_state, from_state
    cdef float score
    for to_state in range(nbase):
        offset_to_state = to_state * nff_state
        # from 0 base state
        curr_fwd[to_state] = curr_logprob[offset_to_state + 0] + prev_fwd[0]
        tb[to_state] = 0
        for from_state in range(1, nff_state):
            # rest of from base states (either flip or flop)
            score = (curr_logprob[offset_to_state + from_state] +
                     prev_fwd[from_state])
            if score > curr_fwd[to_state]:
                curr_fwd[to_state] = score
                tb[to_state] = from_state

    cdef size_t offset_to_flop = nff_state * nbase
    # to flop base state
    for to_state in range(nbase, nff_state):
        # stay state
        curr_fwd[to_state] = prev_fwd[to_state] + curr_logprob[
            offset_to_flop + to_state]
        tb[to_state] = to_state
        # move from flip to flop state
        from_state = to_state - nbase
        score = (curr_logprob[offset_to_flop + from_state] +
                 prev_fwd[from_state])
        if score > curr_fwd[to_state]:
            curr_fwd[to_state] = score
            tb[to_state] = from_state


cdef float flipflop_decode_trans(
    float[:, ::1] tpost,
    size_t nblk,
    size_t nbase,
    uintptr_t[::1] path,
    float[:] qpath
):
    cdef size_t nff_state = nbase + nbase
    cdef size_t ntrans_state = nff_state * (nbase + 1)

    cdef uintptr_t[:, ::1] tb = np.empty((nblk, nff_state), dtype=np.uintp)

    cdef float * curr = <float *> calloc(nff_state, sizeof(float))
    cdef float * prev = <float *> calloc(nff_state, sizeof(float))
    cdef float * tmp
    cdef size_t idx
    for idx in range(nff_state):
        curr[idx] = 0

    # Forwards Viterbi pass
    cdef size_t blk
    for blk in range(nblk):
        # swap curr and prev
        tmp = prev
        prev = curr
        curr = tmp
        decode_forward_step(&tpost[blk,0], nbase, prev, curr, &tb[blk,0])

    cdef size_t max_idx = 0
    cdef float score = curr[0]
    for idx in range(1, nff_state):
        if curr[idx] > score:
            score = curr[idx]
            max_idx = idx

    free(prev)
    free(curr)

    # Traceback
    cdef size_t trans_idx
    path[nblk] = max_idx
    for blk in range(nblk, 0, -1):
        path[blk - 1] = tb[blk - 1, path[blk]]
        trans_idx = (path[blk] * (nbase + nbase) + path[blk - 1]
                     if path[blk] < nbase else
                     nbase * (nbase + nbase) + path[blk - 1])
        qpath[blk] = tpost[blk - 1, trans_idx]
    qpath[0] = np.nan

    return score


cdef inline void decode_forward_step_fb(
        float * curr_logprob, size_t nbase,
        float * prev_fwd, float * curr_fwd):
    cdef size_t nff_state = nbase + nbase

    # to flip base states
    cdef size_t offset_to_state, to_state, from_state
    cdef float score
    for to_state in range(nbase):
        offset_to_state = to_state * nff_state
        # from 0 base state
        curr_fwd[to_state] = curr_logprob[offset_to_state + 0] + prev_fwd[0]
        for from_state in range(1, nff_state):
            # rest of from base states (either flip or flop)
            score = (curr_logprob[offset_to_state + from_state] +
                     prev_fwd[from_state])
            curr_fwd[to_state] = fmaxf(curr_fwd[to_state], score) + log1pf(
                expf(-fabsf(curr_fwd[to_state] - score)))

    cdef size_t offset_to_flop = nff_state * nbase
    # to flop base state
    for to_state in range(nbase, nff_state):
        # stay state
        curr_fwd[to_state] = prev_fwd[to_state] + curr_logprob[
            offset_to_flop + to_state]
        # move from flip to flop state
        from_state = to_state - nbase
        score = (curr_logprob[offset_to_flop + from_state] +
                 prev_fwd[from_state])
        curr_fwd[to_state] = fmaxf(curr_fwd[to_state], score) + log1pf(
            expf(-fabsf(curr_fwd[to_state] - score)))

    cdef float row_logsum, tpost_val
    row_logsum = curr_fwd[0]
    for from_state in range(1, nff_state):
        tpost_val = curr_fwd[from_state]
        row_logsum = fmaxf(row_logsum, tpost_val) + log1pf(
            expf(-fabsf(row_logsum - tpost_val)))

    for to_state in range(nff_state):
        curr_fwd[to_state] -= row_logsum


cdef void decode_forward(
        float[:, ::1] logprob, size_t nbase, size_t nblk,
        float[:, ::1] fwd):
    cdef size_t nff_state = nbase + nbase
    cdef size_t ntrans_state = nff_state * (nbase + 1)

    cdef size_t idx
    for idx in range(nff_state):
        fwd[0, idx] = -np.log(nff_state)

    cdef size_t blk
    for blk in range(nblk):
        decode_forward_step_fb(&logprob[blk,0], nbase, &fwd[blk,0],
                               &fwd[blk + 1,0])


cdef void flipflop_trans_post(
        float[:, ::1] logprob, size_t nbase, size_t nblk, float[:, ::1] tpost):
    cdef size_t nff_state = nbase + nbase
    cdef size_t ntrans_state = nff_state * (nbase + 1)

    cdef float[:, ::1] fwd = np.empty((nblk + 1, nff_state), dtype=np.float32)

    decode_forward(logprob, nbase, nblk, fwd)

    cdef float * curr = <float *> calloc(nff_state, sizeof(float))
    cdef float * prev = <float *> calloc(nff_state, sizeof(float))
    cdef float * tmp
    cdef size_t idx
    for idx in range(nff_state):
        prev[idx] = -HUGE_VALF

    # Backwards pass
    cdef size_t blk, st, to_state, from_state
    cdef float score, row_logsum, tpost_val
    for blk in range(nblk, 0, -1):
        # swap curr and prev
        tmp = prev
        prev = curr
        curr = tmp

        # to flip state
        for to_state in range(nbase):
            for st in range(nff_state):
                tpost[blk - 1, to_state * nff_state + st] = (
                    fwd[blk - 1, st] + prev[to_state] +
                    logprob[blk - 1, to_state * nff_state + st])
        # to flop state
        for to_state in range(nbase, nff_state):
            # from flip
            tpost[blk - 1, nff_state * nbase + to_state - nbase] = (
                fwd[blk - 1, to_state - nbase] + prev[to_state] +
                logprob[blk - 1, nff_state * nbase + to_state - nbase])
            # stay in flop
            tpost[blk - 1, nff_state * nbase + to_state] = (
                fwd[blk - 1, to_state] + prev[to_state] +
                logprob[blk - 1, nff_state * nbase + to_state])

        # Update backwards vector
        # to flop state
        for to_state in range(nbase, nff_state):
            # Move from flip to flop state
            curr[to_state - nbase] = logprob[
                blk - 1, nff_state * nbase + to_state - nbase] + prev[to_state]
            # Stay in flop state
            curr[to_state] = logprob[
                blk - 1, nff_state * nbase + to_state] + prev[to_state]
        # to flip state
        for to_state in range(nbase):
            for from_state in range(nff_state):
                score = (
                    logprob[blk - 1, to_state * nff_state + from_state] +
                    prev[to_state])
                curr[from_state] = fmaxf(curr[from_state], score) + log1pf(
                    expf(-fabsf(curr[from_state] - score)))

        # normalize backwards scores
        row_logsum = curr[0]
        for from_state in range(1, nff_state):
            tpost_val = curr[from_state]
            row_logsum = fmaxf(row_logsum, tpost_val) + log1pf(
                expf(-fabsf(row_logsum - tpost_val)))

        for to_state in range(nff_state):
            curr[to_state] -= row_logsum

    free(curr)
    free(prev)


################################
# Standard flip-flop functions #
################################

def crf_flipflop_trans_post(float[:, ::1] logprob, log=True):
    """ Get posteriors from transition weights
    """
    cdef size_t nblock, nparam, nbase
    nblock, nparam = logprob.shape[0], logprob.shape[1]
    nbase = nstate_to_nbase(nparam)

    tpost_np = np.zeros_like(logprob)
    cdef float[:, ::1] tpost = tpost_np
    flipflop_trans_post(logprob, nbase, nblock, tpost)
    if not log:
        # perform exp inplace
        np.exp(tpost, out=tpost)

    return tpost_np


def crf_flipflop_viterbi(
        float[:, ::1] tpost, uintptr_t[::1] path, float[::1] qpath):
    """ Fast flip-flop Viterbi decoding calling C implementation
    """
    cdef size_t nblock, nparam, nbase
    nblock, nparam = tpost.shape[0], tpost.shape[1]
    nbase = nstate_to_nbase(nparam)

    score = flipflop_decode_trans(tpost, nblock, nbase, path, qpath)

    return score


##################################
# Constrained decoding functions #
##################################

cdef flipflop_constrain_traceback(
    int[::1] int_seq,
    const int[:, ::1] traceback,
    const int[:, ::1] seq_band,
    const int[::1] base_offsets,
    const float[::1] final_fwd_scores
):
    """Perform traceback to extract sequence for most likely path

    Args:
        int_seq (int array): Sequence array to be populated by this function.
        traceback (2D int array): Elements contain flip-flop base to go to at
            the next block back in time. Value equal to index indicates a stay
            in this base.
        seq_band (int array): 2D Int array containing band start and end
            positions for each base. Positions are in num_blocks dimension of
            trans_logprob. Shape should be 2 by seq_len
        base_offsets (int array): Offset to the start of each base in
            traceback/trans_logprob ragged arrays.
        final_fwd_scores (float array): Flip-flop base states at final position.
            Maximum value determines final base in sequence.
    """
    cdef int seq_len = seq_band.shape[1]
    cdef int nblocks = seq_band[1, seq_len - 1]
    cdef int curr_seq_pos = seq_len - 1

    # identify final flip-flop base from final_fwd_scores
    int_seq[curr_seq_pos] = np.argmax(final_fwd_scores)
    cdef int band_pos, curr_block, tb_pos, tb_base
    tb_base = traceback[traceback.shape[0] - 1, int_seq[curr_seq_pos]]
    if tb_base == -1:
        LOGGER.debug(
            "Invalid traceback position at start:\n"
            f"forward score: {np.array(final_fwd_scores)}\n"
            f"traceback: {np.array(traceback[traceback.shape[0] - 1])}"
        )
        raise MegaError(f"Invalid traceback position at start")
    if tb_base != int_seq[curr_seq_pos]:
        curr_seq_pos -= 1
        int_seq[curr_seq_pos] = tb_base
    for curr_block in range(nblocks - 2, 0, -1):
        band_pos = curr_block - seq_band[0, curr_seq_pos]
        if band_pos < 0:
            LOGGER.debug(
                "Constrained basecall traceback outside band (stay)  "
                f"block:{curr_block}/{nblocks-1}  "
                f"seq:{curr_seq_pos}/{seq_len-1}  "
                "band(len="
                f"{seq_band[1, curr_seq_pos] - seq_band[0, curr_seq_pos]}):"
                f"{seq_band[0, curr_seq_pos]}-{seq_band[1, curr_seq_pos]}"
                f"band_pos:{band_pos}"
            )
            raise MegaError(
                "Constrained basecall traceback outside band (stay)"
            )
        if band_pos >= seq_band[1, curr_seq_pos] - seq_band[0, curr_seq_pos]:
            LOGGER.debug(
                "Constrained basecall traceback outside band (step)  "
                f"block:{curr_block}/{nblocks-1}  "
                f"seq:{curr_seq_pos}/{seq_len-1}  "
                "band(len="
                f"{seq_band[1, curr_seq_pos] - seq_band[0, curr_seq_pos]}):"
                f"{seq_band[0, curr_seq_pos]}-{seq_band[1, curr_seq_pos]}"
                f"band_pos:{band_pos}"
            )
            raise MegaError(
                "Constrained basecall traceback outside band (step)"
            )
        tb_pos = base_offsets[curr_seq_pos] + band_pos
        tb_base = traceback[tb_pos, int_seq[curr_seq_pos]]
        if tb_base != int_seq[curr_seq_pos]:
            if tb_base == -1:
                reg_seq_band = seq_band[
                    :,
                    max(0, curr_seq_pos - 3)
                    : min(seq_band.shape[1] - 1, curr_seq_pos + 3)]
                LOGGER.debug(
                    f"Invalid traceback position: "
                    f"block:{curr_block}/{nblocks-1}  "
                    f"seq:{curr_seq_pos}/{seq_len-1}  "
                    "band(len="
                    f"{seq_band[1, curr_seq_pos] - seq_band[0, curr_seq_pos]}):"
                    f"{seq_band[0, curr_seq_pos]}-{seq_band[1, curr_seq_pos]}"
                    f"band_pos:{band_pos}  "
                    f"{tb_base}  "
                    f"{np.array(reg_seq_band)}"
                )
                raise MegaError("Invalid traceback position")
            if curr_seq_pos == 0:
                raise MegaError(
                    "Constrained basecall traceback outside band (seq)"
                )
            curr_seq_pos -= 1
            int_seq[curr_seq_pos] = tb_base
    if curr_seq_pos != 0:
        raise MegaError(
            "Constrained basecall traceback did not reach first base"
        )


cdef inline void constrain_decode_forward_move_step(
    float[::1] pos_scores,
    int[::1] pos_tb,
    const float[::1] pos_logprob,
    const short[::1] pos_allowed_bases,
    const float[::1] prev_base_scores,
    int nbase,
):
    """Compute move forward scores and traceback values for a step into a
    particular block by sequence position. Note that stay scores are computed
    in a separate function as these derive from the previous block in the
    *same* base.

    Args:
        pos_scores (float array): Forward score computed thus far at this
            position.
        pos_tb (int array): Traceback computed thus far at this position.
        pos_logprob (float array): Transition log probabilities into this
            position.
        pos_allowed_bases (int array): Allowed bases at this position.
        prev_base_scores (float array): Forward scores for each flip-flop base
            up to the previous base (from state) and previous block.
        nbase (int): Number of bases in alphabet
    """
    cdef int nff_state = nbase + nbase

    cdef size_t offset_to_state, to_flip_state, to_flop_state, from_state
    cdef float score
    cdef size_t offset_to_flop = nff_state * nbase
    for to_flip_state in range(nbase):
        if not pos_allowed_bases[to_flip_state]:
            # if this base is not allowed skip
            continue
        # move to flip base states
        offset_to_state = to_flip_state * nff_state
        # can move to flip state from any state
        for from_state in range(nff_state):
            # skip stay computation. Computed from prev_block_scores in
            # constrained decoding.
            if to_flip_state == from_state:
                continue
            score = (pos_logprob[offset_to_state + from_state]
                     + prev_base_scores[from_state])
            if score > pos_scores[to_flip_state]:
                pos_scores[to_flip_state] = score
                pos_tb[to_flip_state] = from_state

        # move to flop base state from flip state
        from_state = to_flip_state
        to_flop_state = to_flip_state + nbase
        # can only move from flip to flop state
        score = (pos_logprob[offset_to_flop + from_state]
                 + prev_base_scores[from_state])
        if score > pos_scores[to_flop_state]:
            pos_scores[to_flop_state] = score
            pos_tb[to_flop_state] = from_state


cdef inline void constrain_decode_forward_stay_step(
    float[::1] pos_scores,
    int[::1] pos_tb,
    const float[::1] pos_logprob,
    const short[::1] pos_allowed_bases,
    const float[::1] prev_block_scores,
    int nbase,
):
    """Compute stay forward scores and traceback values for a step into a
    particular block by sequence position. Note that move scores are computed
    in a separate function as these derive from the previous block in the
    *previous* base.

    Args:
        pos_scores (float array): Forward score computed thus far at this
            position.
        pos_tb (int array): Traceback computed thus far at this position.
        pos_logprob (float array): Transition log probabilities into this
            position.
        pos_allowed_bases (int array): Allowed bases at this position.
        prev_block_scores (float array): Forward scores for each flip-flop base
            up to the previous block in this same base(from state).
        nbase (int): Number of bases in alphabet
    """
    cdef size_t nff_state = nbase + nbase
    cdef size_t offset_to_flop = nff_state * nbase
    cdef size_t flip_base, flop_base
    cdef float score
    for flip_base in range(nbase):
        if not pos_allowed_bases[flip_base]:
            # if this base is not allowed retain previous score and pos_tb
            continue
        # stay in flip base state
        score = (pos_logprob[(flip_base * nff_state) + flip_base]
                 + prev_block_scores[flip_base])
        if score > pos_scores[flip_base]:
            pos_scores[flip_base] = score
            pos_tb[flip_base] = flip_base
        # stay in flop base state
        flop_base = nbase + flip_base
        score = (prev_block_scores[flop_base]
                 + pos_logprob[offset_to_flop + flop_base])
        if score > pos_scores[flop_base]:
            pos_scores[flop_base] = score
            pos_tb[flop_base] = flop_base


cdef void constrain_banded_forward_vit_step(
    float[:, ::1] band_scores,
    int[:, ::1] band_tb,
    const float[:, ::1] band_logprob,
    const float[:, ::1] prev_scores,
    const short[::1] pos_allowed_bases,
    int band_start_diff,
    int nbase
):
    """Process one base using standard Viterbi path flip-flop scoring.

    Args:
        band_logprob (2D float array): Flip-flop transition log probabilities
            corresponding to same band as band_scores.
        band_scores (2D float array): To be populated by this function. Elements
            will be forward scores at each position for each flip-flop state
        band_tb (2D int array): To be populated by this function. Elements will
            be flip-flop base to go to at the next step. Value equal to index
            indicates a stay in this base.
        prev_scores (const 2D float array): Forward scores for the previous base
        band_start_diff (int): Difference in starting coordinates between
            current and previous base
        pos_allowed_bases (const 1D int array): Bases allowed
    """
    cdef int nff_state = nbase + nbase
    cdef int band_pos
    cdef float base_score, move_score, stay_score
    # compute start position in band
    # if the previous band starts before this base's band compute band start
    # state; else invalid values remain at start of band
    if band_start_diff > 0:
        constrain_decode_forward_move_step(
            pos_scores=band_scores[0, :],
            pos_tb=band_tb[0, :],
            pos_logprob=band_logprob[0, :],
            pos_allowed_bases=pos_allowed_bases,
            prev_base_scores=prev_scores[band_start_diff - 1, :],
            nbase=nbase,
        )
        # clip prev_scores to start at same position as band_scores
        prev_scores = prev_scores[band_start_diff:, :]
    # if base bands are the same clip last element of prev_scores
    if prev_scores.shape[0] == band_scores.shape[0]:
        prev_scores = prev_scores[:prev_scores.shape[0] - 1]

    # compute scores where curr and prev base overlap
    for band_pos in range(1, prev_scores.shape[0] + 1):
        # compute max over move and stay scores
        constrain_decode_forward_move_step(
            pos_scores=band_scores[band_pos, :],
            pos_tb=band_tb[band_pos, :],
            pos_logprob=band_logprob[band_pos, :],
            pos_allowed_bases=pos_allowed_bases,
            prev_base_scores=prev_scores[band_pos - 1, :],
            nbase=nbase,
        )
        constrain_decode_forward_stay_step(
            pos_scores=band_scores[band_pos, :],
            pos_tb=band_tb[band_pos, :],
            pos_logprob=band_logprob[band_pos, :],
            pos_allowed_bases=pos_allowed_bases,
            prev_block_scores=band_scores[band_pos - 1, :],
            nbase=nbase,
        )

    # stay through rest of the band
    for band_pos in range(prev_scores.shape[0] + 1, band_scores.shape[0]):
        constrain_decode_forward_stay_step(
            pos_scores=band_scores[band_pos, :],
            pos_tb=band_tb[band_pos, :],
            pos_logprob=band_logprob[band_pos, :],
            pos_allowed_bases=pos_allowed_bases,
            prev_block_scores=band_scores[band_pos - 1, :],
            nbase=nbase,
        )


def flipflop_constrain_decode(
    const float[:, ::1] trans_logprob,
    const short[:, ::1] all_allowed_bases,
    const int[:, ::1] seq_band
):
    """Decode flip-flop transitions with constains on the produced sequence. In
    particular the sequence length is fixed and any number of bases within the
    fixed sequence may be constrained to any particular set of bases.

    This decoding is also constrain to a band through the sequence to
    transitions matrix output. This is to increase efficiency.

    Args:
        trans_logprob (np.array): 2D Float array containing flip-flop transition
            log probabilities. Shape should be num_blocks by num_transitions.
            num_blocks is signal // stride and num_transitions is the number of
            flip-flop transitions (40 for 4 canonical bases).
        all_allowed_bases (np.array): 2D Boolean array (np.uint8) containing
            allowed bases for each position within the produced basecalls.
            Shape should be seq_len by num_bases. True indicies indicate that
            this base is allowed at this position in the sequence.
        seq_band (np.array): 2D Int array containing band start and end
            positions for each base. Positions are in num_blocks dimension of
            trans_logprob. Shape should be 2 by seq_len

    Returns:
        Integer encoded contrained decoded sequence.
    """
    cdef int nblock, nstate, seq_len
    nblock = trans_logprob.shape[0]
    nstate = trans_logprob.shape[1]
    seq_len = all_allowed_bases.shape[0]
    if nblock < seq_band[1, seq_band.shape[1] - 1]:
        raise MegaError("Band extends beyond transitions posteriors array.")
    cdef int nbase = nstate_to_nbase(nstate)
    cdef int nff_state = nbase + nbase
    if nbase != all_allowed_bases.shape[1]:
        raise MegaError(
            "Number of bases implied by trans_logprob and bases included in "
            "all_allowed_bases do not agree."
        )
    if 2 != seq_band.shape[0]:
        raise MegaError("Second band dimension should be length 2.")
    if seq_len != seq_band.shape[1]:
        raise MegaError(
            "Sequence constraints length and band length do not agree."
        )
    if nblock < seq_len:
        raise MegaError(
            "num_blocks < seq_len results in no valid sequence decodings."
        )

    # compute offset within all_scores and traceback arrays to the start of
    # each base in the band
    cdef int[::1] base_offsets = np.empty(seq_len + 1, dtype=np.int32)
    base_offsets[0] = 0
    # cumsum
    cdef int base_len, base_pos
    for base_pos, base_len in enumerate(np.diff(seq_band, axis=0)[0]):
        base_offsets[base_pos + 1] = base_offsets[base_pos] + base_len

    # initialize ragged-banded score and traceback arrays. First dimension
    # represented ragged positions with block positions contiguous in memory.
    cdef int tot_band_len = base_offsets[base_offsets.shape[0] - 1]
    cdef float[:, ::1] all_scores = np.full(
        (tot_band_len, nff_state), -HUGE_VALF, dtype=np.float32)
    cdef int[:, ::1] traceback = np.full(
        (tot_band_len, nff_state), -1, dtype=np.int32)

    cdef int curr_bw, start_base, flip_stay_idx, flop_stay_idx, band_pos
    curr_bw = seq_band[1, 0]
    # initialize first base band with allowed bases
    for start_base in range(nbase):
        if not all_allowed_bases[0, start_base]:
            continue
        traceback[0, start_base] = start_base
        traceback[0, start_base + nbase] = start_base + nbase
        flip_stay_idx = (start_base * nff_state) + start_base
        flop_stay_idx = (nff_state * nbase) + nbase + start_base
        all_scores[0, start_base] = trans_logprob[0, flip_stay_idx]
        all_scores[0, start_base + nbase] = trans_logprob[0, flop_stay_idx]
        for band_pos in range(1, curr_bw):
            all_scores[band_pos, start_base] = (
                all_scores[band_pos - 1, start_base]
                + trans_logprob[band_pos, flip_stay_idx])
            all_scores[band_pos, start_base + nbase] = (
                all_scores[band_pos - 1, start_base + nbase]
                + trans_logprob[band_pos, flop_stay_idx])
            traceback[band_pos, start_base] = start_base
            traceback[band_pos, start_base + nbase] = start_base + nbase

    cdef int prev_offset, prev_band_st, prev_bw, seq_pos
    cdef int curr_offset, curr_band_st, curr_band_en
    # loop over sequence positions computing forward scores and traceback
    prev_bw = curr_bw
    prev_band_st = prev_offset = 0
    for seq_pos in range(1, seq_len):
        curr_band_st = seq_band[0, seq_pos]
        curr_band_en = seq_band[1, seq_pos]
        curr_bw = curr_band_en - curr_band_st
        curr_offset = base_offsets[seq_pos]

        # compute all scores and traceback for this base
        constrain_banded_forward_vit_step(
            all_scores[curr_offset:curr_offset + curr_bw, :],
            traceback[curr_offset:curr_offset + curr_bw, :],
            trans_logprob[curr_band_st : curr_band_en],
            all_scores[prev_offset:prev_offset + prev_bw, :],
            all_allowed_bases[seq_pos, :],
            curr_band_st - prev_band_st,
            nbase,
        )

        prev_band_st = curr_band_st
        prev_bw = curr_bw
        prev_offset = curr_offset

    # prepare output seq array for decoding
    int_seq = np.empty(seq_len, dtype=np.int32)
    # TODO also save and return the path
    flipflop_constrain_traceback(
        int_seq,
        traceback,
        seq_band,
        base_offsets,
        all_scores[all_scores.shape[0] - 1, :]
    )

    return int_seq


####################
# Sequence scoring #
####################

cdef float score_best_path(
        float[:, ::1] tpost, uintptr_t[::1] seq,
        size_t tpost_start, size_t tpost_end, size_t nseq, size_t nbase):
    cdef size_t ntrans_state = (nbase + nbase) * (nbase + 1)
    cdef size_t nblk = tpost_end - tpost_start
    cdef size_t window_width = nblk - nseq + 2

    cdef size_t * stay_indices = <size_t *> calloc(nseq, sizeof(size_t))
    cdef size_t * step_indices = <size_t *> calloc(nseq - 1, sizeof(size_t))
    # get_stay_step(seq, nseq, nbase, stay_indices, step_indices)
    cdef size_t * flop_mask_states = <size_t *> calloc(nseq, sizeof(size_t))
    stay_indices[0] = (seq[0] * (nbase + nbase) + seq[0]
                       if seq[0] < nbase else
                       nbase * (nbase + nbase) + seq[0])
    flop_mask_states[0] = seq[0]
    cdef size_t seq_pos, prev_fm_state, curr_fm_state
    for seq_pos in range(1, nseq):
        prev_fm_state = flop_mask_states[seq_pos - 1]
        if seq[seq_pos] == prev_fm_state:
            curr_fm_state = seq[seq_pos] + nbase
        else:
            curr_fm_state = seq[seq_pos]
        flop_mask_states[seq_pos] = curr_fm_state
        stay_indices[seq_pos] = (
            curr_fm_state * (nbase + nbase) + curr_fm_state
            if curr_fm_state < nbase else
            nbase * (nbase + nbase) + curr_fm_state)
        step_indices[seq_pos - 1] = (
            curr_fm_state * (nbase + nbase) + prev_fm_state
            if curr_fm_state < nbase else
            nbase * (nbase + nbase) + prev_fm_state)
    free(flop_mask_states)

    cdef float * curr = <float *> calloc(window_width, sizeof(float))
    cdef float * prev = <float *> calloc(window_width, sizeof(float))
    cdef float * tmp

    cdef size_t win_pos, trans_pos, stay_idx, step_idx
    cdef float stay_score, step_score, wp_score, pos_tpost
    # cumsum over stay in first seq pos to init prev
    prev[0] = 0.0
    stay_idx = stay_indices[0]
    for win_pos in range(1, window_width):
        wp_score = prev[win_pos - 1]
        trans_pos = tpost_start + win_pos - 1
        prev[win_pos] = wp_score + tpost[trans_pos, stay_idx]

    for seq_pos in range(1, nseq):
        wp_score = prev[0]
        stay_idx = stay_indices[seq_pos]
        step_idx = step_indices[seq_pos - 1]
        pos_tpost = tpost[tpost_start + seq_pos - 1, step_idx]
        curr[0] = wp_score + pos_tpost
        for win_pos in range(1, window_width):
            trans_pos = tpost_start + seq_pos + win_pos - 1

            # compute step score
            wp_score = prev[win_pos]
            pos_tpost = tpost[trans_pos, step_idx]
            step_score = wp_score + pos_tpost

            # compute stay score
            wp_score = curr[win_pos - 1]
            pos_tpost = tpost[trans_pos, stay_idx]
            stay_score = wp_score + pos_tpost

            # store best path score
            if step_score > stay_score:
                curr[win_pos] = step_score
            else:
                curr[win_pos] = stay_score
        # swap prev and curr scores
        tmp = prev
        prev = curr
        curr = tmp

    cdef float score = prev[window_width - 1]
    free(stay_indices)
    free(step_indices)
    free(curr)
    free(prev)

    return score


cdef float score_all_paths(
        float[:, ::1] tpost, uintptr_t[:] seq,
        size_t tpost_start, size_t tpost_end, size_t nseq, size_t nbase):
    cdef size_t ntrans_state = (nbase + nbase) * (nbase + 1)
    cdef size_t nblk = tpost_end - tpost_start
    cdef size_t window_width = nblk - nseq + 2

    cdef size_t * stay_indices = <size_t *> calloc(nseq, sizeof(size_t));
    cdef size_t * step_indices = <size_t *> calloc(nseq - 1, sizeof(size_t));
    # get_stay_step(seq, nseq, nbase, stay_indices, step_indices)
    cdef size_t * flop_mask_states = <size_t *> calloc(nseq, sizeof(size_t))
    stay_indices[0] = (seq[0] * (nbase + nbase) + seq[0]
                       if seq[0] < nbase else
                       nbase * (nbase + nbase) + seq[0])
    flop_mask_states[0] = seq[0]
    cdef size_t seq_pos, prev_fm_state, curr_fm_state
    for seq_pos in range(1, nseq):
        prev_fm_state = flop_mask_states[seq_pos - 1]
        if seq[seq_pos] == prev_fm_state:
            curr_fm_state = seq[seq_pos] + nbase
        else:
            curr_fm_state = seq[seq_pos]
        flop_mask_states[seq_pos] = curr_fm_state
        stay_indices[seq_pos] = (
            curr_fm_state * (nbase + nbase) + curr_fm_state
            if curr_fm_state < nbase else
            nbase * (nbase + nbase) + curr_fm_state)
        step_indices[seq_pos - 1] = (
            curr_fm_state * (nbase + nbase) + prev_fm_state
            if curr_fm_state < nbase else
            nbase * (nbase + nbase) + prev_fm_state)
    free(flop_mask_states)

    cdef float * curr = <float *> calloc(window_width, sizeof(float))
    cdef float * prev = <float *> calloc(window_width, sizeof(float))
    cdef float * tmp

    cdef size_t win_pos, trans_pos, stay_idx, step_idx
    cdef float stay_score, step_score, wp_score, pos_tpost
    # cumsum over stay in first seq pos to init prev
    prev[0] = 0.0
    stay_idx = stay_indices[0]
    for win_pos in range(1, window_width):
        wp_score = prev[win_pos - 1]
        trans_pos = tpost_start + win_pos - 1
        prev[win_pos] = wp_score + tpost[trans_pos, stay_idx]

    for seq_pos in range(1, nseq):
        wp_score = prev[0]
        stay_idx = stay_indices[seq_pos]
        step_idx = step_indices[seq_pos - 1]
        pos_tpost = tpost[tpost_start + seq_pos - 1, step_idx]
        curr[0] = wp_score + pos_tpost
        for win_pos in range(1, window_width):
            trans_pos = tpost_start + seq_pos + win_pos - 1

            # compute step score
            wp_score = prev[win_pos]
            pos_tpost = tpost[trans_pos, step_idx]
            step_score = wp_score + pos_tpost

            # compute stay score
            wp_score = curr[win_pos - 1]
            pos_tpost = tpost[trans_pos, stay_idx]
            stay_score = wp_score + pos_tpost

            # store Viterbi all paths score
            curr[win_pos] = fmaxf(stay_score, step_score) + log1pf(
                expf(-fabsf(stay_score - step_score)))
        # swap prev and curr scores
        tmp = prev
        prev = curr
        curr = tmp

    cdef float score = prev[window_width - 1]
    free(stay_indices)
    free(step_indices)
    free(curr)
    free(prev)

    return score


def score_seq(
        float[:, ::1] tpost, uintptr_t[::1] seq,
        tpost_start, tpost_end, all_paths):
    nseq = seq.shape[0]
    nbase = nstate_to_nbase(tpost.shape[1])
    if all_paths:
        return score_all_paths(tpost, seq, tpost_start, tpost_end, nseq, nbase)
    return score_best_path(tpost, seq, tpost_start, tpost_end, nseq, nbase)


########################
# Mod Sequence scoring #
########################

cdef float score_best_path_mod(
        float[:, ::1] tpost, uintptr_t[::1] seq,
        uintptr_t[::1] mod_cats,
        uintptr_t[::1] can_mods_offsets, size_t tpost_start,
        size_t tpost_end, size_t nseq, size_t nstate):
    cdef size_t nbase = nstate_to_nbase(nstate - can_mods_offsets[4])
    cdef size_t ntrans_state = (nbase + nbase) * (nbase + 1)
    cdef size_t nblk = tpost_end - tpost_start
    cdef size_t window_width = nblk - nseq + 2

    cdef size_t * stay_indices = <size_t *> calloc(nseq, sizeof(size_t));
    cdef size_t * step_indices = <size_t *> calloc(nseq - 1, sizeof(size_t));
    #get_stay_step(seq, nseq, nbase, stay_indices, step_indices)
    cdef size_t * flop_mask_states = <size_t *> calloc(nseq, sizeof(size_t))
    stay_indices[0] = (seq[0] * (nbase + nbase) + seq[0]
                       if seq[0] < nbase else
                       nbase * (nbase + nbase) + seq[0])
    flop_mask_states[0] = seq[0]
    cdef size_t seq_pos, prev_fm_state, curr_fm_state
    for seq_pos in range(1, nseq):
        prev_fm_state = flop_mask_states[seq_pos - 1]
        if seq[seq_pos] == prev_fm_state:
            curr_fm_state = seq[seq_pos] + nbase
        else:
            curr_fm_state = seq[seq_pos]
        flop_mask_states[seq_pos] = curr_fm_state
        stay_indices[seq_pos] = (
            curr_fm_state * (nbase + nbase) + curr_fm_state
            if curr_fm_state < nbase else
            nbase * (nbase + nbase) + curr_fm_state)
        step_indices[seq_pos - 1] = (
            curr_fm_state * (nbase + nbase) + prev_fm_state
            if curr_fm_state < nbase else
            nbase * (nbase + nbase) + prev_fm_state)
    free(flop_mask_states)

    cdef float * curr = <float *> calloc(window_width, sizeof(float))
    cdef float * prev = <float *> calloc(window_width, sizeof(float))
    cdef float * tmp

    cdef size_t win_pos, trans_pos, stay_idx, step_idx, mod_idx
    cdef float stay_score, step_score, wp_score, pos_tpost, pos_mod_tpost
    # cumsum over stay in first seq pos to init prev
    prev[0] = 0.0
    stay_idx = stay_indices[0]
    for win_pos in range(1, window_width):
        wp_score = prev[win_pos - 1]
        trans_pos = tpost_start + win_pos - 1
        prev[win_pos] = wp_score + tpost[trans_pos, stay_idx]

    for seq_pos in range(1, nseq):
        wp_score = prev[0]
        stay_idx = stay_indices[seq_pos]
        step_idx = step_indices[seq_pos - 1]
        pos_tpost = tpost[tpost_start + seq_pos - 1, step_idx]
        mod_idx = ntrans_state + can_mods_offsets[seq[seq_pos]] + \
                  mod_cats[seq_pos]
        pos_mod_tpost = tpost[tpost_start + seq_pos - 1, mod_idx]
        curr[0] = wp_score + pos_tpost + pos_mod_tpost
        for win_pos in range(1, window_width):
            trans_pos = tpost_start + seq_pos + win_pos - 1

            # compute step score
            wp_score = prev[win_pos]
            pos_tpost = tpost[trans_pos, step_idx]
            pos_mod_tpost = tpost[trans_pos, mod_idx]
            step_score = wp_score + pos_tpost + pos_mod_tpost

            # compute stay score
            wp_score = curr[win_pos - 1]
            pos_tpost = tpost[trans_pos, stay_idx]
            stay_score = wp_score + pos_tpost

            # store best path score
            if step_score > stay_score:
                curr[win_pos] = step_score
            else:
                curr[win_pos] = stay_score
        # swap prev and curr scores
        tmp = prev
        prev = curr
        curr = tmp

    cdef float score = prev[window_width - 1]
    free(stay_indices)
    free(step_indices)
    free(curr)
    free(prev)

    return score


cdef float score_all_paths_mod(
        float[:, ::1] tpost, uintptr_t[::1] seq,
        uintptr_t[::1] mod_cats, uintptr_t[::1] can_mods_offsets,
        size_t tpost_start, size_t tpost_end, size_t nseq, size_t nstate):
    cdef size_t nbase = nstate_to_nbase(nstate - can_mods_offsets[4])
    cdef size_t ntrans_state = (nbase + nbase) * (nbase + 1)
    cdef size_t nblk = tpost_end - tpost_start
    cdef size_t window_width = nblk - nseq + 2

    cdef size_t * stay_indices = <size_t *> calloc(nseq, sizeof(size_t));
    cdef size_t * step_indices = <size_t *> calloc(nseq - 1, sizeof(size_t));
    #get_stay_step(seq, nseq, nbase, stay_indices, step_indices)
    cdef size_t * flop_mask_states = <size_t *> calloc(nseq, sizeof(size_t))
    stay_indices[0] = (seq[0] * (nbase + nbase) + seq[0]
                       if seq[0] < nbase else
                       nbase * (nbase + nbase) + seq[0])
    flop_mask_states[0] = seq[0]
    cdef size_t seq_pos, prev_fm_state, curr_fm_state
    for seq_pos in range(1, nseq):
        prev_fm_state = flop_mask_states[seq_pos - 1]
        if seq[seq_pos] == prev_fm_state:
            curr_fm_state = seq[seq_pos] + nbase
        else:
            curr_fm_state = seq[seq_pos]
        flop_mask_states[seq_pos] = curr_fm_state
        stay_indices[seq_pos] = (
            curr_fm_state * (nbase + nbase) + curr_fm_state
            if curr_fm_state < nbase else
            nbase * (nbase + nbase) + curr_fm_state)
        step_indices[seq_pos - 1] = (
            curr_fm_state * (nbase + nbase) + prev_fm_state
            if curr_fm_state < nbase else
            nbase * (nbase + nbase) + prev_fm_state)
    free(flop_mask_states)

    cdef float * curr = <float *> calloc(window_width, sizeof(float))
    cdef float * prev = <float *> calloc(window_width, sizeof(float))
    cdef float * tmp

    cdef size_t win_pos, trans_pos, stay_idx, step_idx, mod_idx
    cdef float stay_score, step_score, wp_score, pos_tpost, pos_mod_tpost
    # cumsum over stay in first seq pos to init prev
    prev[0] = 0.0
    stay_idx = stay_indices[0]
    for win_pos in range(1, window_width):
        wp_score = prev[win_pos - 1]
        trans_pos = tpost_start + win_pos - 1
        prev[win_pos] = wp_score + tpost[trans_pos, stay_idx]

    for seq_pos in range(1, nseq):
        wp_score = prev[0]
        stay_idx = stay_indices[seq_pos]
        step_idx = step_indices[seq_pos - 1]
        pos_tpost = tpost[tpost_start + seq_pos - 1, step_idx]
        mod_idx = ntrans_state + can_mods_offsets[seq[seq_pos]] + \
                  mod_cats[seq_pos]
        pos_mod_tpost = tpost[tpost_start + seq_pos - 1, mod_idx]
        curr[0] = wp_score + pos_tpost + pos_mod_tpost
        for win_pos in range(1, window_width):
            trans_pos = tpost_start + seq_pos + win_pos - 1

            # compute step score
            wp_score = prev[win_pos]
            pos_tpost = tpost[trans_pos, step_idx]
            pos_mod_tpost = tpost[trans_pos, mod_idx]
            step_score = wp_score + pos_tpost + pos_mod_tpost

            # compute stay score
            wp_score = curr[win_pos - 1]
            pos_tpost = tpost[trans_pos, stay_idx]
            stay_score = wp_score + pos_tpost

            # store Viterbi all paths score
            curr[win_pos] = fmaxf(stay_score, step_score) + log1pf(
                expf(-fabsf(stay_score - step_score)))
        # swap prev and curr scores
        tmp = prev
        prev = curr
        curr = tmp

    cdef float score = prev[window_width - 1]
    free(stay_indices)
    free(step_indices)
    free(curr)
    free(prev)

    return score

cdef float get_best_flipflop_path_core(
    int[::1] path,
    const float[:, ::1] tpost,
    const uintptr_t[::1] seq,
    size_t nbase
):
    cdef int nseq = seq.shape[0]
    cdef int nblk = tpost.shape[0]
    cdef int window_width = nblk - nseq + 1

    cdef int[:, ::1] tb = -np.ones((nblk, nseq), dtype=np.int32)
    cdef float[::1] prev = np.empty(window_width, dtype=np.float32)
    cdef float[::1] curr = np.empty(window_width, dtype=np.float32)

    # cumsum over stay in first seq pos to init prev
    prev[0] = tpost[0, stay_idx]
    tb[0, 0] = 0
    # always start in flip base
    cdef int prev_ff_base = seq[0]
    cdef int stay_idx = prev_ff_base * (nbase + nbase) + prev_ff_base
    cdef int win_pos
    for win_pos in range(1, window_width):
        prev[win_pos] = prev[win_pos - 1] + tpost[win_pos, stay_idx]
        tb[win_pos, 0] = win_pos

    cdef int seq_pos, step_idx, trans_pos
    cdef float stay_score, step_score
    for seq_pos in range(1, nseq):
        # compute flip-flop stay and step indices
        if prev_ff_base < nbase and seq[seq_pos] == seq[seq_pos - 1]:
            # if this is a flop base
            stay_idx = nbase * (nbase + nbase) + seq[seq_pos] + nbase
            step_idx = nbase * (nbase + nbase) + seq[seq_pos]
            prev_ff_base = seq[seq_pos] + nbase
        else:
            # else this is a flip base
            stay_idx = seq[seq_pos] * (nbase + nbase) + seq[seq_pos]
            step_idx = seq[seq_pos] * (nbase + nbase) + prev_ff_base
            prev_ff_base = seq[seq_pos]

        # step into first position in window
        curr[0] = prev[0] + tpost[seq_pos, step_idx]
        tb[seq_pos, seq_pos] = 0
        for win_pos in range(1, window_width):
            trans_pos = seq_pos + win_pos
            # compute stay and step scores
            stay_score = curr[win_pos - 1] + tpost[trans_pos, stay_idx]
            step_score = prev[win_pos] + tpost[trans_pos, step_idx]
            # store best path score
            if step_score > stay_score:
                tb[trans_pos, seq_pos] = 0
                curr[win_pos] = step_score
            else:
                # add to stays from previous position in base
                tb[trans_pos, seq_pos] = tb[trans_pos - 1, seq_pos] + 1
                curr[win_pos] = stay_score
        # swap prev and curr
        curr, prev = prev, curr

    # perform traceback
    cdef int block_pos = nblk
    path[nseq] = nblk
    cdef int blocks_to_base_start
    for seq_pos in range(nseq - 1, -1, -1):
        # look one block back from block_pos to find the last position in the
        # path attributed to previous base. The value is the number of blocks
        # back to the start of this previous base
        blocks_to_base_start = tb[block_pos - 1, seq_pos]
        if blocks_to_base_start < 0:
            raise MegaError("Invalid traceback")
        block_pos -= blocks_to_base_start + 1
        path[seq_pos] = block_pos

    return prev[window_width - 1]


def get_best_flipflop_path(
    const float[:, ::1] tpost,
    const uintptr_t[::1] seq,
    size_t nbase,
):
    """Compute the path (blocks to seq mapping) with the best score.

    Args:
        tpost (np.array): 2D float32 array containing transition probabilities.
            Only canonical base transition probabilities are used, but modified
            base probabilities may be included in tpost and will be ignored.
        seq (np.array): 1D int32 array containing integer encoded sequence
        nbase (int): Number of bases included in alphabet for seq

    Returns:
        2-tuple containing float score for the path and a numpy array
        containing the path. The path will be one element longer than seq.
        Each element represents the index within tpost that the corresponding
        position in seq starts.
    """
    if tpost.shape[0] < seq.shape[0]:
        raise MegaError(
            "Cannot compute path. Transition probabilities than sequence "
            "length."
        )
    path = np.empty(seq.shape[0] + 1, dtype=np.int32)
    score = get_best_flipflop_path_core(path, tpost, seq, nbase)
    return score, path


def score_mod_seq(
    float[:, ::1] tpost,
    uintptr_t[::1] seq,
    uintptr_t[::1] mod_cats,
    uintptr_t[::1] can_mods_offsets,
    tpost_start,
    tpost_end,
    all_paths,
):
    nseq = seq.shape[0]
    nstate = tpost.shape[1]
    if all_paths:
        return score_all_paths_mod(tpost, seq, mod_cats, can_mods_offsets,
                                   tpost_start, tpost_end, nseq, nstate)
    return score_best_path_mod(tpost, seq, mod_cats, can_mods_offsets,
                               tpost_start, tpost_end, nseq, nstate)


################################################
# Categorical modification flip-flop functions #
################################################

@cython.wraparound(True)
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


@cython.wraparound(True)
def decode_post(
        float[:, ::1] r_post,
        can_alphabet=ALPHABET,
        float[:, ::1] mod_weights=None,
        int64_t[::1] can_nmods=None):
    """Decode a posterior using Viterbi algorithm for transducer.
    :param r_post: numpy array containing transducer posteriors.
    :param can_alphabet: canonical alphabet corresponding to flip-flop labels.
    :returns: tuple containing (base calls, score and raw block positions).
    """
    nblock, nstate = r_post.shape[:2]
    cdef size_t n_can_base = len(can_alphabet)
    if n_can_base != nstate_to_nbase(nstate):
        raise NotImplementedError(
            'Incompatible decoding alphabet and posterior states.')

    path = np.zeros(nblock + 1, dtype=np.uintp)
    qpath = np.zeros(nblock + 1, dtype=np.float32)

    score = crf_flipflop_viterbi(r_post, path, qpath)

    # only process positions "transitioned into" along the path
    # first position doesn't have a score anyways
    # This aligned the indices of path and the posterior matricies
    runval, runlen = rle(path)
    basecall = ''.join(can_alphabet[int(b) % n_can_base] for b in runval)
    rl_cumsum = np.cumsum(np.concatenate([[0], runlen]))

    mods_scores = None
    if mod_weights is not None:
        mods_scores = np.empty(
            (runval.shape[0], len(can_alphabet) + sum(can_nmods)),
            dtype=np.float32
        )
        mods_scores[0] = 0
        mods_scores[1:] = mod_weights[rl_cumsum[1:-1] - 1]

    return basecall, score, rl_cumsum, mods_scores
