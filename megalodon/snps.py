import sys
import queue
import sqlite3
import numpy as np
from time import sleep
import multiprocessing as mp
from collections import defaultdict, namedtuple

from megalodon import decode, megalodon_helper as mh


FIELD_NAMES = ('read_id', 'chrm', 'strand', 'pos', 'score',
               'ref_seq', 'alt_seq', 'snp_id')
SNP_DATA = namedtuple('SNP_DATA', FIELD_NAMES)
CREATE_snps_TBLS = """
CREATE TABLE snps (
    {} TEXT,
    {} TEXT,
    {} INTEGER,
    {} INTEGER,
    {} FLOAT,
    {} TEXT,
    {} TEXT,
    {} TEXT
)""".format(*FIELD_NAMES)
ADDMANY_SNPS = "INSERT INTO snps VALUES (?,?,?,?,?,?,?,?)"
CREATE_SNPS_IDX = "CREATE INDEX snp_pos ON snps (chrm, strand, pos)"

COUNT_UNIQ_POS = """
SELECT COUNT(*) FROM (
SELECT DISTINCT chrm, strand, pos FROM snps)"""
SEL_UNIQ_POS = 'SELECT DISTINCT chrm, strand, pos FROM snps'
SEL_POS_STATS = 'SELECT * FROM snps WHERE chrm=? AND strand=? AND pos=?'


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
            snp_ref_seq = mh.revcomp(snp_ref_seq)
            snp_alt_seq = mh.revcomp(snp_alt_seq)
        ovlp_snps.append((read_pos, snp_ref_seq, snp_alt_seq, snp_id))

    return ovlp_snps

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
                              mh.revcomp(snp_ref_seq))
        fwd_strand_alt_seq = (snp_alt_seq if r_ref_pos.strand == 1 else
                              mh.revcomp(snp_alt_seq))
        if len(fwd_strand_ref_seq) == 0 or len(fwd_strand_alt_seq) == 0:
            fwd_ref_base = (
                mh.ALPHABET[r_ref_seq[pos_bb - 1]] if r_ref_pos.strand == 1 else
                mh.revcomp(mh.ALPHABET[r_ref_seq[pos_bb + len(snp_ref_seq)]]))
            fwd_strand_ref_seq = fwd_ref_base + fwd_strand_ref_seq
            fwd_strand_alt_seq = fwd_ref_base + fwd_strand_alt_seq
        r_snp_calls.append((
            snp_ref_pos, loc_ref_score - loc_alt_score, fwd_strand_ref_seq,
            fwd_strand_alt_seq, snp_id))

    return r_snp_calls

def _get_snps_queue(snps_q, snps_conn, snp_id_tbl, snps_db_fn, snps_txt_fn):
    snps_db = sqlite3.connect(snps_db_fn)
    snps_db_c = snps_db.cursor()
    snps_db_c.execute(CREATE_SNPS_TBLS)
    snps_txt_fp = None if snps_txt_fn is None else open(snps_txt_fn, 'w')

    while True:
        try:
            # note strand is +1 for fwd or -1 for rev
            r_snp_calls, (read_id, chrm, strand) = snps_q.get(block=False)
            snps_db_c.executemany(ADDMANY_SNPS, [
                (read_id, chrm, strand, pos, score, ref_seq, alt_seq,
                 snp_id_tbl[chrm][snp_i])
                for pos, score, ref_seq, alt_seq, snp_i in r_snp_calls])
            if snps_txt_fp is not None:
                # would involve batching and creating several conversion tables
                # for var strings (read_if and chrms).
                snps_txt_fp.write('\n'.join((
                    '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                        read_id, chrm, strand, pos, score, ref_seq, alt_seq,
                        snp_id_tbl[chrm][snp_i])
                    for pos, score, ref_seq, alt_seq, snp_i in r_snp_calls)) +
                              '\n')
                snps_txt_fp.flush()
        except queue.Empty:
            if snps_conn.poll():
                break
            sleep(0.1)
            continue

    while not snps_q.empty():
        r_snp_calls, (read_id, chrm, strand) = snps_q.get(block=False)
        snps_db_c.executemany(ADDMANY_SNPS, [
            (read_id, chrm, strand, pos, score, ref_seq, alt_seq,
             snp_id_tbl[chrm][snp_i])
            for pos, score, ref_seq, alt_seq, snp_i in r_snp_calls])
        if snps_txt_fp is not None:
            snps_txt_fp.write('\n'.join((
                '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                    read_id, chrm, strand, pos, score, ref_seq, alt_seq,
                    snp_id_tbl[chrm][snp_i])
                for pos, score, ref_seq, alt_seq, snp_i in r_snp_calls)) + '\n')
            snps_txt_fp.flush()
    if snps_txt_fp is not None: snps_txt_fp.close()
    snps_db_c.execute(CREATE_SNPS_IDX)
    snps_db.commit()
    snps_db.close()

    return

class Snps(object):
    def __init__(
            self, snp_fn, do_prepend_chr_vcf, max_snp_size, all_paths,
            write_snps_txt):
        self.all_paths = all_paths
        self.write_snps_txt = write_snps_txt
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


class AggSnps(object):
    def __init__(self, snps_db_fn):
        # open as read only database
        self.snps_db = sqlite3.connect(snps_db_fn, uri=True)
        self.snps_db_c = snps_db.cursor()
        self.n_uniq_pos = self.snps_db_c.execute(COUNT_UNIQ_POS).fetchone()[0]

    def iter_uniq_pos(self):
        for q_val in self.snps_db_c.execute(SEL_UNIQ_POS):
            yield q_val

    def get_pos_stats(self, chrm, strand, pos):
        return [SNP_DATA(snp_stats) for snp_stats in self.snps_db_c.execute(
            SEL_POS_STATS, (chrm, strand, pos))]

    def close(self):
        self.snps_db.close()


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
