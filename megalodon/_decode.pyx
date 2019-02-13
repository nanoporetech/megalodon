cimport _libdecode
import cython
import numpy as np
cimport numpy as np


@cython.boundscheck(False)
@cython.wraparound(False)
def nstate_to_nbase(size_t nstate):
    """Convert number of flip-flop states to number of bases and check for
    valid states number.
    """
    cdef np.float32_t nbase_f = np.sqrt(0.25 + (0.5 * nstate)) - 0.5
    assert np.mod(nbase_f, 1) == 0, (
        'Number of states not valid for flip-flop model. ' +
        'nstates: {}\tconverted nbases: {}').format(nstate, nbase_f)
    cdef size_t nbase = <size_t>nbase_f
    return nbase


##########################################
###### Standard flip-flop functions ######
##########################################

@cython.boundscheck(False)
@cython.wraparound(False)
def crf_flipflop_trans_post(np.ndarray[np.float32_t, ndim=2, mode="c"] trans,
                            log=True):
    """ Get posteriors from transition weights
    """
    cdef size_t nblock, nparam, nbase
    nblock, nparam = trans.shape[0], trans.shape[1]
    nbase = nstate_to_nbase(nparam)

    cdef np.ndarray[np.float32_t, ndim=2, mode="c"] tpost = np.zeros_like(trans)
    _libdecode.flipflop_trans_post(&trans[0,0], nbase, nblock, &tpost[0,0])
    if not log:
        # perform exp inplace
        np.exp(tpost, tpost)

    return tpost


@cython.boundscheck(False)
@cython.wraparound(False)
def crf_flipflop_viterbi(np.ndarray[np.float32_t, ndim=2, mode="c"] trans,
                         np.ndarray[np.uintp_t, ndim=1, mode="c"] path,
                         np.ndarray[np.float32_t, ndim=1, mode="c"] qpath):
    """ Fast flip-flop Viterbi decoding calling C implementation
    """
    cdef size_t nblock, nparam, nbase
    nblock, nparam = trans.shape[0], trans.shape[1]
    nbase = nstate_to_nbase(nparam)

    score = _libdecode.flipflop_decode_trans(
        &trans[0,0], nblock, nbase, &path[0], &qpath[0])

    return score


@cython.boundscheck(False)
@cython.wraparound(False)
def score_seq(
        np.ndarray[np.float32_t, ndim=2, mode="c"] tpost,
        np.ndarray[np.uintp_t, ndim=1, mode="c"] seq,
        tpost_start, tpost_end, all_paths):
    nseq = seq.shape[0]
    nbase = nstate_to_nbase(tpost.shape[1])
    if all_paths:
        return _libdecode.score_all_paths(
            &tpost[0,0], &seq[0], tpost_start, tpost_end, nseq, nbase)
    return _libdecode.score_best_path(
        &tpost[0,0], &seq[0], tpost_start, tpost_end, nseq, nbase)


##########################################################
###### Categorical modification flip-flop functions ######
##########################################################

@cython.boundscheck(False)
@cython.wraparound(False)
def score_mod_seq(
        np.ndarray[np.float32_t, ndim=2, mode="c"] tpost,
        np.ndarray[np.uintp_t, ndim=1, mode="c"] seq,
        np.ndarray[np.uintp_t, ndim=1, mode="c"] mod_cats,
        np.ndarray[np.uintp_t, ndim=1, mode="c"] can_mods_offsets,
        tpost_start, tpost_end, all_paths):
    nseq = seq.shape[0]
    nstate = tpost.shape[1]
    if all_paths:
        return _libdecode.score_all_paths_mod(
            &tpost[0,0], &seq[0], &mod_cats[0], &can_mods_offsets[0],
            tpost_start, tpost_end, nseq, nstate)
    return _libdecode.score_best_path_mod(
        &tpost[0,0], &seq[0], &mod_cats[0], &can_mods_offsets[0],
        tpost_start, tpost_end, nseq, nstate)
