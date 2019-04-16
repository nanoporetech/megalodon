cimport _libdecode
import cython
import numpy as np
cimport numpy as np


ALPHABET = 'ACGT'


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

def decode_post(r_post, alphabet=ALPHABET, mod_weights=None, can_nmods=None):
    """Decode a posterior using Viterbi algorithm for transducer.
    :param r_post: numpy array containing transducer posteriors.
    :param alphabet: alphabet corresponding to flip-flop labels.
    :returns: tuple containing (base calls, score and raw block positions).
    """
    nblock, nstate = r_post.shape[:2]
    nbase = len(set(alphabet))
    if nbase != nstate_to_nbase(nstate):
        raise NotImplementedError(
            'Incompatible decoding alphabet and posterior states.')

    path = np.zeros(nblock + 1, dtype=np.uintp)
    qpath = np.zeros(nblock + 1, dtype=np.float32)

    score = crf_flipflop_viterbi(r_post, path, qpath)

    # only process positions "transitioned into" along the path
    # first position doesn't have a score anyways
    # This aligned the indices of path and the posterior matricies
    runval, runlen = rle(path)
    basecall = ''.join(alphabet[int(b) % nbase] for b in runval)
    rl_cumsum = np.cumsum(np.concatenate([[0], runlen]))

    mods_scores = None
    if mod_weights is not None:
        # TODO cythonize this function
        # extract modified base probabilities for each modification included in
        # the input model
        # don't test first base since it is never "moved into"
        # and subtract 1 to align "moved into" indices
        bc_matched_mod_weights = mod_weights[rl_cumsum[1:-1] - 1]
        curr_can_pos = 0
        mods_scores = []
        for base_i, can_nmod in enumerate(can_nmods):
            if can_nmod > 0:
                base_poss = np.where(np.equal(np.mod(
                    runval[1:], len(can_nmods)), base_i))[0]
            for mod_i in range(can_nmod):
                mod_i_scores = np.full(runval.shape, np.NAN)
                mod_i_scores[base_poss + 1] = bc_matched_mod_weights[
                    base_poss, curr_can_pos + 1 + mod_i]
                mods_scores.append(mod_i_scores)
            curr_can_pos += 1 + can_nmod
        mods_scores = np.stack(mods_scores, axis=1)

    return basecall, score, rl_cumsum, mods_scores
