from libc.stdint cimport int32_t

cdef extern from "_c_decode.h":
    float flipflop_decode_trans(
        const float * trans, size_t nblk, size_t nbase, size_t * path,
        float * qpath);

    void flipflop_trans_post(
        const float * logprob, size_t nbase, size_t nblk, float * tpost);

    float score_best_path(
        const float * tpost, size_t * seq, size_t tpost_start,
        size_t tpost_end, size_t nseq, size_t nbase);

    float score_all_paths(
        const float * tpost, size_t * seq, size_t tpost_start,
        size_t tpost_end, size_t nseq, size_t nbase);

    float score_best_path_mod(
        const float * tpost, size_t * seq, const size_t * mod_cats,
        const size_t * can_mods_offsets, size_t tpost_start,
        size_t tpost_end, size_t nseq, size_t nstate);

    float score_all_paths_mod(
        const float * tpost, size_t * seq, const size_t * mod_cats,
        const size_t * can_mods_offsets, size_t tpost_start,
        size_t tpost_end, size_t nseq, size_t nstate);
