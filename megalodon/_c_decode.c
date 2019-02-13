#include <assert.h>
#include <immintrin.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "_c_decode.h"

#define LARGE_VAL 1e30f


/*
****************
Helper functions
****************
 */

static inline size_t nstate_to_nbase(size_t ntrans_state){
    double nbase_d = sqrt(0.25 + (0.5 * ntrans_state)) - 0.5;
    assert(fmod(nbase_d, 1.0) == 0.0);
    return (size_t) round(nbase_d);
}


static inline float logsumexpf(float x, float y, float a){
    return fmaxf(x, y) + log1pf(expf(-a * fabsf(x-y))) / a;
}


int argmaxf(float const * x, size_t n) {
    assert(n > 0);
    if (NULL == x) {
        return -1;
    }
    size_t imax = 0;
    float vmax = x[0];
    for (size_t i = 1; i < n; i++) {
        if (x[i] > vmax) {
            vmax = x[i];
            imax = i;
        }
    }
    return imax;
}

float valmaxf(float const * x, size_t n) {
    assert(n > 0);
    if (NULL == x) {
        return NAN;
    }
    float vmax = x[0];
    for (size_t i = 1; i < n; i++) {
        if (x[i] > vmax) {
            vmax = x[i];
        }
    }
    return vmax;
}

inline size_t trans_lookup(size_t from, size_t to, size_t nbase){
    assert(nbase >= 0);
    assert(from >= 0 && from < nbase + nbase);
    assert(to >= 0 && to < nbase + nbase);
    assert(to < nbase || ((to % nbase) == (from % nbase)));

    const size_t nff_state = nbase + nbase;
    const size_t offset = nbase * nff_state;

    return (to < nbase) ? (to * nff_state + from) : (offset + from);
}


void log_row_normalise_inplace(float * mat, size_t ncol, size_t nrow){
    if(NULL == mat){
        return;
    }

    for (size_t col=0 ; col < ncol; col++) {
        const size_t offset = col * nrow;
        float row_logsum = mat[offset];
        for(size_t row=1 ; row < nrow ; row++){
            row_logsum = logsumexpf(row_logsum, mat[offset + row], 1.0);
        }

        for(size_t row=0 ; row < nrow ; row++){
            mat[offset + row] -= row_logsum;
        }
    }
}


/*
******************
Decoding functions
******************
 */

void decode_forward_step(float const * curr_scores, size_t nbase,
                         float * prev_fwd, float * curr_fwd, size_t * tb){
    if(NULL == curr_scores || NULL == prev_fwd || NULL == curr_fwd ||
       NULL == tb){
        return;
    }

    const size_t nff_state = nbase + nbase;

    // ib -- fl[i]p to base states
    for(size_t ib=0 ; ib < nbase ; ib++){
        const size_t offset_to_state = ib * nff_state;
        // from 0 base state
        curr_fwd[ib] = curr_scores[offset_to_state + 0] + prev_fwd[0];
        tb[ib] = 0;
        for(size_t from_state=1 ; from_state < nff_state ; from_state++){
            // rest of from base states (either flip or flop)
            const float score = curr_scores[offset_to_state + from_state] +
                prev_fwd[from_state];
            if(score > curr_fwd[ib]){
                curr_fwd[ib] = score;
                tb[ib] = from_state;
            }
        }
    }
    const size_t offset_to_flop = nff_state * nbase;
    // ob -- fl[o]p to base state
    for(size_t ob=nbase ; ob < nff_state ; ob++){
        // stay state
        curr_fwd[ob] = prev_fwd[ob] + curr_scores[offset_to_flop + ob];
        tb[ob] = ob;
        // move from flip to flop state
        const size_t from_base = ob - nbase;
        const float score = curr_scores[offset_to_flop + from_base] +
            prev_fwd[from_base];
        if(score > curr_fwd[ob]){
            curr_fwd[ob] = score;
            tb[ob] = from_base;
        }
    }
}


float flipflop_decode_trans(float const * trans, size_t nblk, size_t nbase,
                            size_t * path, float * qpath){
    if(NULL == trans || NULL == path || NULL == qpath){
        return NAN;
    }

    const size_t nff_state = nbase + nbase;
    const size_t ntrans_state = nff_state * (nbase + 1);

    float * mem = calloc(2 * nff_state, sizeof(float));
    size_t * tb = calloc(nff_state * nblk, sizeof(size_t));
    if(NULL == mem || NULL == tb){
        free(tb);
        free(mem);
        return NAN;
    }

    float * curr = mem;
    float * prev = mem + nff_state;
    for(size_t init_state=0 ; init_state < nff_state ; init_state++){
        curr[init_state] = 0;
    }

    //  Forwards Viterbi pass
    for(size_t blk=0 ; blk < nblk ; blk++){
        {   // Swap
            float * tmp = curr;
            curr = prev;
            prev = tmp;
        }
        float const * curr_trans = trans + blk * ntrans_state;
        decode_forward_step(curr_trans, nbase, prev, curr,
                            tb + blk * nff_state);
    }

    //  Traceback
    const float score = valmaxf(curr, nff_state);
    path[nblk] = argmaxf(curr, nff_state);
    for(size_t blk=nblk ; blk > 0 ; blk--){
        path[blk - 1] = tb[nff_state * (blk - 1) + path[blk]];
        qpath[blk] =
            trans[ntrans_state * (blk - 1) +
                  trans_lookup(path[blk-1], path[blk], nbase)];
    }
    qpath[0] = NAN;

    free(tb);
    free(mem);

    return score;
}


void decode_forward(float const * scores, size_t nbase, size_t nblk,
                    float * fwd){
    const size_t nff_state = nbase + nbase;
    const size_t ntrans_state = nff_state * (nbase + 1);

    for(size_t init_state=0 ; init_state < nff_state ; init_state++){
        fwd[init_state] = 0;
    }

    size_t * tmp_tb = calloc(nff_state, sizeof(size_t));
    if(NULL == tmp_tb){
        free(tmp_tb);
        return;
    }
    //  Forwards Viterbi pass
    for(size_t blk=0 ; blk < nblk ; blk++){
        float * prev_fwd = fwd + blk * nff_state;
        float * curr_fwd = fwd + (blk + 1) * nff_state;
        float const * curr_scores = scores + blk * ntrans_state;
        decode_forward_step(curr_scores, nbase, prev_fwd, curr_fwd, tmp_tb);
    }

    free(tmp_tb);
}


void flipflop_trans_post(float const * logprob, size_t nbase, size_t nblk,
                         float * tpost){
    const size_t nff_state = nbase + nbase;
    const size_t ntrans_state = nff_state * (nbase + 1);

    float * fwd = calloc(nff_state * (nblk + 1), sizeof(float));
    float * mem = calloc(2 * nff_state, sizeof(float));
    if(NULL == fwd || NULL == mem){
        free(fwd);
        free(mem);
        return;
    }

    decode_forward(logprob, nbase, nblk, fwd);
    float * prev = mem;
    float * curr = mem + nff_state;
    for(size_t init_state=0 ; init_state < nff_state ; init_state++){
        curr[init_state] = 0;
    }
    for(size_t init_state=0 ; init_state < nff_state ; init_state++){
        prev[init_state] = -LARGE_VAL;
    }

    //  Backwards pass
    for(size_t blk=nblk ; blk > 0 ; blk--){
        const size_t offset_fwd = (blk - 1) * nff_state;
        const size_t offset = (blk - 1) * ntrans_state;
        const size_t offset_flop = offset + nff_state * nbase;

        {  // Swap
           float * tmp = prev;
           prev = curr;
           curr = tmp;
        }

        for(size_t ib=0 ; ib < nbase ; ib++){
            // End up in flip state
            const size_t offset_state = offset + ib * nff_state;
            for(size_t st=0 ; st < nff_state ; st++){
                tpost[offset_state + st] = fwd[offset_fwd + st] + prev[ib]
                    + logprob[offset_state + st];
            }
        }
        for(size_t ob=nbase ; ob < nff_state ; ob++){
            // End up in flop state
            const size_t ib = ob - nbase;
            tpost[offset_flop + ob] = fwd[offset_fwd + ob] + prev[ob]
                + logprob[offset_flop + ob];
            tpost[offset_flop + ib] = fwd[offset_fwd + ib] + prev[ob]
                + logprob[offset_flop + ib];
        }

        // Update backwards vector
        // ob -- fl[o]p to base state
        for(size_t ob=nbase ; ob < nff_state ; ob++){
            const size_t from_base = ob - nbase;
            // Stay in flop state
            curr[ob] = logprob[offset_flop + ob] + prev[ob];
            // Move from flip to flop state
            curr[from_base] = logprob[offset_flop + from_base] + prev[ob];
        }
        // ib -- fl[i]p to base states
        for(size_t ib=0 ; ib < nbase ; ib++){
            const size_t offset_state = offset + ib * nff_state;
            for(size_t from_state=0 ; from_state < nff_state ; from_state++){
                // rest of from base states (either flip or flop)
                const float score = logprob[offset_state + from_state] +
                    prev[ib];
                if(score > curr[from_state]){
                    curr[from_state] = score;
                }
            }
        }
    }

    free(mem);
    free(fwd);

    log_row_normalise_inplace(tpost, nblk, ntrans_state);
}


/*
****************
Sequence scoring
****************
 */


float score_best_path(float const * tpost, size_t * seq, size_t tpost_start,
                      size_t tpost_end, size_t nseq, size_t nbase){
    if(NULL == tpost){
        return NAN;
    }

    size_t ntrans_state = (nbase + nbase) * (nbase + 1);
    float score = NAN;
    size_t * flop_mask_states = calloc(nseq, sizeof(size_t));
    size_t * stay_indices = calloc(nseq, sizeof(size_t));
    size_t * step_indices = calloc(nseq - 1, sizeof(size_t));
    if(NULL == flop_mask_states || NULL == step_indices ||
       NULL == stay_indices){
        goto cleanup1;
    }
    stay_indices[0] = trans_lookup(seq[0], seq[0], nbase);
    flop_mask_states[0] = seq[0];
    for(size_t seq_pos=1 ; seq_pos < nseq ; seq_pos++){
        if(seq[seq_pos] == flop_mask_states[seq_pos - 1]){
            flop_mask_states[seq_pos] = seq[seq_pos] + nbase;
        } else {
            flop_mask_states[seq_pos] = seq[seq_pos];
        }
        stay_indices[seq_pos] =
            trans_lookup(flop_mask_states[seq_pos],
                         flop_mask_states[seq_pos], nbase);
        step_indices[seq_pos - 1] =
            trans_lookup(flop_mask_states[seq_pos - 1],
                         flop_mask_states[seq_pos], nbase);
    }

    const size_t nblk = tpost_end - tpost_start;
    const size_t window_width = nblk - nseq + 2;

    float * curr_scores = calloc(window_width, sizeof(float));
    float * prev_scores = calloc(window_width, sizeof(float));
    if(NULL == curr_scores || NULL == prev_scores){
        goto cleanup2;
    }

    // cumsum over stay in first seq pos to init prev_scores
    prev_scores[0] = 0.0f;
    for(size_t win_blk=1 ; win_blk < window_width ; win_blk++){
        prev_scores[win_blk] =
            prev_scores[win_blk - 1] +
            tpost[((tpost_start + win_blk - 1) * ntrans_state) +
                  stay_indices[0]];
    }

    for(size_t seq_pos=1 ; seq_pos < nseq ; seq_pos++){
        curr_scores[0] = prev_scores[0] +
            tpost[((tpost_start + seq_pos - 1) *
                   ntrans_state) + step_indices[seq_pos - 1]];
        for(size_t win_pos=1 ; win_pos < window_width ; win_pos++){
            const size_t blk_offset = (tpost_start + seq_pos + win_pos - 1) *
                ntrans_state;
            // step score
            curr_scores[win_pos] = prev_scores[win_pos] +
                tpost[blk_offset + step_indices[seq_pos - 1]];
            // stay score
            const float pos_score = curr_scores[win_pos - 1] +
                tpost[blk_offset + stay_indices[seq_pos]];
            // store best path score
            if(pos_score > curr_scores[win_pos]){
                curr_scores[win_pos] = pos_score;
            }
        }
        {   // Swap curr and prev vectors
            float * tmp = curr_scores;
            curr_scores = prev_scores;
            prev_scores = tmp;
        }
    }

    score = prev_scores[window_width - 1];

 cleanup2:
    free(prev_scores);
    free(curr_scores);
 cleanup1:
    free(flop_mask_states);
    free(stay_indices);
    free(step_indices);

    return score;
}


float score_all_paths(float const * tpost, size_t * seq, size_t tpost_start,
                      size_t tpost_end, size_t nseq, size_t nbase){
    if(NULL == tpost){
        return NAN;
    }

    size_t ntrans_state = (nbase + nbase) * (nbase + 1);
    float score = NAN;
    size_t * flop_mask_states = calloc(nseq, sizeof(size_t));
    size_t * stay_indices = calloc(nseq, sizeof(size_t));
    size_t * step_indices = calloc(nseq - 1, sizeof(size_t));
    if(NULL == flop_mask_states || NULL == step_indices ||
       NULL == stay_indices){
        goto cleanup1;
    }
    stay_indices[0] = trans_lookup(seq[0], seq[0], nbase);
    flop_mask_states[0] = seq[0];
    for(size_t seq_pos=1 ; seq_pos < nseq ; seq_pos++){
        if(seq[seq_pos] == flop_mask_states[seq_pos - 1]){
            flop_mask_states[seq_pos] = seq[seq_pos] + nbase;
        } else {
            flop_mask_states[seq_pos] = seq[seq_pos];
        }
        stay_indices[seq_pos] =
            trans_lookup(flop_mask_states[seq_pos],
                         flop_mask_states[seq_pos], nbase);
        step_indices[seq_pos - 1] =
            trans_lookup(flop_mask_states[seq_pos - 1],
                         flop_mask_states[seq_pos], nbase);
    }

    const size_t nblk = tpost_end - tpost_start;
    const size_t window_width = nblk - nseq + 2;

    float * curr_scores = calloc(window_width, sizeof(float));
    float * prev_scores = calloc(window_width, sizeof(float));
    if(NULL == curr_scores || NULL == prev_scores){
        goto cleanup2;
    }

    // cumsum over stay in first seq pos to init prev_scores
    prev_scores[0] = 0.0f;
    for(size_t win_blk=1 ; win_blk < window_width ; win_blk++){
        prev_scores[win_blk] =
            prev_scores[win_blk - 1] +
            tpost[((tpost_start + win_blk - 1) * ntrans_state) +
                  stay_indices[0]];
    }

    for(size_t seq_pos=1 ; seq_pos < nseq ; seq_pos++){
        curr_scores[0] = prev_scores[0] +
            tpost[((tpost_start + seq_pos - 1) *
                   ntrans_state) + step_indices[seq_pos - 1]];
        for(size_t win_pos=1 ; win_pos < window_width ; win_pos++){
            const size_t blk_offset = (tpost_start + seq_pos + win_pos - 1) *
                ntrans_state;
            // step score
            curr_scores[win_pos] = prev_scores[win_pos] +
                tpost[blk_offset + step_indices[seq_pos - 1]];
            // stay score
            const float pos_score = curr_scores[win_pos - 1] +
                tpost[blk_offset + stay_indices[seq_pos]];
            // store Viterbi all paths score
            curr_scores[win_pos] =
                logsumexpf(curr_scores[win_pos], pos_score, 1.0);
        }
        {   // Swap curr and prev vectors
            float * tmp = curr_scores;
            curr_scores = prev_scores;
            prev_scores = tmp;
        }
    }

    score = prev_scores[window_width - 1];

 cleanup2:
    free(prev_scores);
    free(curr_scores);
 cleanup1:
    free(flop_mask_states);
    free(stay_indices);
    free(step_indices);

    return score;
}


/*
********************
Mod Sequence scoring
********************
 */

float score_best_path_mod(float const * tpost, size_t * seq, size_t * mod_cats,
                          size_t * can_mods_offsets, size_t tpost_start,
                          size_t tpost_end, size_t nseq, size_t nstate){
    if(NULL == tpost || NULL == seq || NULL == mod_cats ||
       NULL == can_mods_offsets){
        return NAN;
    }

    size_t nbase = nstate_to_nbase(nstate - can_mods_offsets[4]);
    size_t ntrans_state = (nbase + nbase) * (nbase + 1);
    float score = NAN;
    size_t * flop_mask_states = calloc(nseq, sizeof(size_t));
    size_t * stay_indices = calloc(nseq, sizeof(size_t));
    size_t * step_indices = calloc(nseq - 1, sizeof(size_t));
    if(NULL == flop_mask_states || NULL == step_indices ||
       NULL == stay_indices){
        goto cleanup1;
    }
    stay_indices[0] = trans_lookup(seq[0], seq[0], nbase);
    flop_mask_states[0] = seq[0];
    for(size_t seq_pos=1 ; seq_pos < nseq ; seq_pos++){
        if(seq[seq_pos] == flop_mask_states[seq_pos - 1]){
            flop_mask_states[seq_pos] = seq[seq_pos] + nbase;
        } else {
            flop_mask_states[seq_pos] = seq[seq_pos];
        }
        stay_indices[seq_pos] =
            trans_lookup(flop_mask_states[seq_pos],
                         flop_mask_states[seq_pos], nbase);
        step_indices[seq_pos - 1] =
            trans_lookup(flop_mask_states[seq_pos - 1],
                         flop_mask_states[seq_pos], nbase);
    }

    const size_t nblk = tpost_end - tpost_start;
    const size_t window_width = nblk - nseq + 2;

    float * curr_scores = calloc(window_width, sizeof(float));
    float * prev_scores = calloc(window_width, sizeof(float));
    if(NULL == curr_scores || NULL == prev_scores){
        goto cleanup2;
    }

    // cumsum over stay in first seq pos to init prev_scores
    prev_scores[0] = 0.0f;
    for(size_t win_blk=1 ; win_blk < window_width ; win_blk++){
        prev_scores[win_blk] =
            prev_scores[win_blk - 1] +
            tpost[((tpost_start + win_blk - 1) * nstate) +
                  stay_indices[0]];
    }

    for(size_t seq_pos=1 ; seq_pos < nseq ; seq_pos++){
        const size_t blk0_offset = (tpost_start + seq_pos - 1) * nstate;
        curr_scores[0] = prev_scores[0] +
            tpost[blk0_offset + step_indices[seq_pos - 1]] +
            tpost[blk0_offset + ntrans_state +
                  can_mods_offsets[seq[seq_pos]] + mod_cats[seq_pos]];
        for(size_t win_pos=1 ; win_pos < window_width ; win_pos++){
            const size_t blk_offset = (tpost_start + seq_pos + win_pos - 1) *
                nstate;
            // step score
            curr_scores[win_pos] = prev_scores[win_pos] +
                tpost[blk_offset + step_indices[seq_pos - 1]] +
                tpost[blk_offset + ntrans_state +
                      can_mods_offsets[seq[seq_pos]] + mod_cats[seq_pos]];
            // stay score
            const float pos_score = curr_scores[win_pos - 1] +
                tpost[blk_offset + stay_indices[seq_pos]];
            // store best path score
            if(pos_score > curr_scores[win_pos]){
                curr_scores[win_pos] = pos_score;
            }
        }
        {   // Swap curr and prev vectors
            float * tmp = curr_scores;
            curr_scores = prev_scores;
            prev_scores = tmp;
        }
    }

    score = prev_scores[window_width - 1];

 cleanup2:
    free(prev_scores);
    free(curr_scores);
 cleanup1:
    free(flop_mask_states);
    free(stay_indices);
    free(step_indices);

    return score;
}


float score_all_paths_mod(float const * tpost, size_t * seq, size_t * mod_cats,
                          size_t * can_mods_offsets, size_t tpost_start,
                          size_t tpost_end, size_t nseq, size_t nstate){
    if(NULL == tpost || NULL == seq || NULL == mod_cats ||
       NULL == can_mods_offsets){
        return NAN;
    }

    size_t nbase = nstate_to_nbase(nstate - can_mods_offsets[4]);
    size_t ntrans_state = (nbase + nbase) * (nbase + 1);
    float score = NAN;
    size_t * flop_mask_states = calloc(nseq, sizeof(size_t));
    size_t * stay_indices = calloc(nseq, sizeof(size_t));
    size_t * step_indices = calloc(nseq - 1, sizeof(size_t));
    if(NULL == flop_mask_states || NULL == step_indices ||
       NULL == stay_indices){
        goto cleanup1;
    }
    stay_indices[0] = trans_lookup(seq[0], seq[0], nbase);
    flop_mask_states[0] = seq[0];
    for(size_t seq_pos=1 ; seq_pos < nseq ; seq_pos++){
        if(seq[seq_pos] == flop_mask_states[seq_pos - 1]){
            flop_mask_states[seq_pos] = seq[seq_pos] + nbase;
        } else {
            flop_mask_states[seq_pos] = seq[seq_pos];
        }
        stay_indices[seq_pos] =
            trans_lookup(flop_mask_states[seq_pos],
                         flop_mask_states[seq_pos], nbase);
        step_indices[seq_pos - 1] =
            trans_lookup(flop_mask_states[seq_pos - 1],
                         flop_mask_states[seq_pos], nbase);
    }

    const size_t nblk = tpost_end - tpost_start;
    const size_t window_width = nblk - nseq + 2;

    float * curr_scores = calloc(window_width, sizeof(float));
    float * prev_scores = calloc(window_width, sizeof(float));
    if(NULL == curr_scores || NULL == prev_scores){
        goto cleanup2;
    }

    // cumsum over stay in first seq pos to init prev_scores
    prev_scores[0] = 0.0f;
    for(size_t win_blk=1 ; win_blk < window_width ; win_blk++){
        prev_scores[win_blk] =
            prev_scores[win_blk - 1] +
            tpost[((tpost_start + win_blk - 1) * nstate) +
                  stay_indices[0]];
    }

    for(size_t seq_pos=1 ; seq_pos < nseq ; seq_pos++){
        const size_t blk0_offset = (tpost_start + seq_pos - 1) * nstate;
        curr_scores[0] = prev_scores[0] +
            tpost[blk0_offset + step_indices[seq_pos - 1]] +
            tpost[blk0_offset + ntrans_state +
                  can_mods_offsets[seq[seq_pos]] + mod_cats[seq_pos]];
        for(size_t win_pos=1 ; win_pos < window_width ; win_pos++){
            const size_t blk_offset = (tpost_start + seq_pos + win_pos - 1) *
                nstate;
            // step score
            curr_scores[win_pos] = prev_scores[win_pos] +
                tpost[blk_offset + step_indices[seq_pos - 1]] +
                tpost[blk_offset + ntrans_state +
                      can_mods_offsets[seq[seq_pos]] + mod_cats[seq_pos]];
            // stay score
            const float pos_score = curr_scores[win_pos - 1] +
                tpost[blk_offset + stay_indices[seq_pos]];
            // store Viterbi all paths score
            curr_scores[win_pos] =
                logsumexpf(curr_scores[win_pos], pos_score, 1.0);
        }
        {   // Swap curr and prev vectors
            float * tmp = curr_scores;
            curr_scores = prev_scores;
            prev_scores = tmp;
        }
    }

    score = prev_scores[window_width - 1];

 cleanup2:
    free(prev_scores);
    free(curr_scores);
 cleanup1:
    free(flop_mask_states);
    free(stay_indices);
    free(step_indices);

    return score;
}
