import queue
import sqlite3
from time import sleep

import numpy as np

from megalodon import decode, megalodon_helper as mh


CREATE_MODS_TBLS = """
CREATE TABLE mods (
    read_id TEXT,
    chrm TEXT,
    strand INTEGER,
    pos INTEGER,
    score FLOAT,
    mod_base TEXT,
    motif TEXT
)"""
ADDMANY_MODS = "INSERT INTO mods VALUES (?,?,?,?,?,?,?)"
CREATE_MODS_IDX = "CREATE INDEX mod_pos ON mods (chrm, strand, pos)"


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

def _get_mods_queue(mods_q, mods_conn, mods_db_fn, mods_txt_fn):
    mods_db = sqlite3.connect(mods_db_fn)
    mods_db_c = mods_db.cursor()
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


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
