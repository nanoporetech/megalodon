#include <stdint.h>

float flipflop_decode_trans(float const * trans, size_t nblk, size_t nbase,
                            size_t * path, float * qpath);

void flipflop_trans_post(float const * logprob, size_t nbase, size_t nblk,
                         float * tpost);

float score_best_path(float const * tpost, size_t * seq, size_t tpost_start,
                      size_t tpost_end, size_t nseq, size_t nbase);

float score_all_paths(float const * tpost, size_t * seq, size_t tpost_start,
                      size_t tpost_end, size_t nseq, size_t nbase);

float score_best_path_mod(float const * tpost, size_t * seq, size_t * mod_cats,
                          size_t * can_mods_offsets, size_t tpost_start,
                          size_t tpost_end, size_t nseq, size_t nstate);

float score_all_paths_mod(float const * tpost, size_t * seq, size_t * mod_cats,
                          size_t * can_mods_offsets, size_t tpost_start,
                          size_t tpost_end, size_t nseq, size_t nstate);
