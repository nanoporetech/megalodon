import sys
import queue
import sqlite3
import datetime
from time import sleep
import multiprocessing as mp
from collections import defaultdict, namedtuple, OrderedDict

import numpy as np

from megalodon import decode, megalodon_helper as mh
from megalodon._version import MEGALODON_VERSION


DIPLOID_MODE = 'diploid'
HAPLIOD_MODE = 'haploid'
AMBIG_LLRS_THRESH = None

FIELD_NAMES = ('read_id', 'chrm', 'strand', 'pos', 'score',
               'ref_seq', 'alt_seq', 'snp_id')
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
    {} TEXT
)""".format(*FIELD_NAMES)

SET_NO_ROLLBACK_MODE='PRAGMA journal_mode = OFF'
SET_ASYNC_MODE='PRAGMA synchronous = OFF'

ADDMANY_SNPS = "INSERT INTO snps VALUES (?,?,?,?,?,?,?,?)"
CREATE_SNPS_IDX = '''
CREATE INDEX snp_pos ON snps (chrm, pos, ref_seq, alt_seq, snp_id)'''

COUNT_UNIQ_SNPS = """
SELECT COUNT(*) FROM (
SELECT DISTINCT chrm, pos, snp_id, ref_seq, alt_seq FROM snps)"""
SEL_UNIQ_SNP_ID = '''
SELECT DISTINCT chrm, pos, snp_id, ref_seq, alt_seq FROM snps'''
SEL_SNP_STATS = '''
SELECT * FROM snps WHERE chrm=? AND pos=? AND
snp_id=? AND ref_seq=? AND alt_seq=?'''

FIXED_VCF_MI = [
    'phasing=none',
    'INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
    'FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    'FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
    'FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
]


########################
##### SNP Encoding #####
########################

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

def simplify_and_encode_snp(
        snp_ref_seq, snp_alt_seq, ref_pos, max_snp_size, trim_matching=False):
    """ Simplify SNP when extra bases are included (for indels)
    """
    snp_ref_seq = snp_ref_seq.upper()
    snp_alt_seq = snp_alt_seq.upper()
    if any(b == ',' for b in snp_alt_seq):
        raise mh.MegaError('Multiple alt SNPs not currently supported.')
    # handle cases containing non-base values (e.g. dash for deletion)
    # assume ref or alt deletion, but only if it is a single character
    if any(rb not in mh.ALPHABET for rb in snp_ref_seq):
        if len(snp_ref_seq) > 1:
            raise mh.MegaError('Invalid SNP')
        if any(ab not in mh.ALPHABET for ab in snp_alt_seq):
            raise mh.MegaError('Invalid SNP')
        if max_snp_size is not None and len(snp_alt_seq) > max_snp_size:
            raise mh.MegaError('SNP too long')
        return 0, encode_snp_seq(snp_alt_seq), ref_pos
    elif any(ab not in mh.ALPHABET for ab in snp_alt_seq):
        if len(snp_alt_seq) > 1:
            raise mh.MegaError('Invalid SNP')
        if max_snp_size is not None and len(snp_ref_seq) > max_snp_size:
            raise mh.MegaError('SNP too long')
        return encode_snp_seq(snp_ref_seq), 0, ref_pos

    if trim_matching:
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

    if (max_snp_size is not None and
        np.abs(len(snp_ref_seq) - len(snp_alt_seq)) > max_snp_size):
        raise mh.MegaError('SNP too long')
    return encode_snp_seq(snp_ref_seq), encode_snp_seq(snp_alt_seq), ref_pos


################################
##### Per-read SNP Scoring #####
################################

def iter_overlapping_snps(r_ref_pos, snps_to_test, edge_buffer):
    """Iterator over SNPs overlapping the read mapped position.

    SNPs within edge buffer of the end of the mapping will be ignored.
    """
    try:
        chrm_poss, chrm_ref_es, chrm_alt_es = snps_to_test[r_ref_pos.chrm]
    except KeyError:
        raise mh.MegaError(
            'No SNPs on mapped chromosome/record (see --prepend-chr-*)')
    start_idx, end_idx = np.searchsorted(
        chrm_poss, (r_ref_pos.start + edge_buffer, r_ref_pos.end - edge_buffer))
    if start_idx >= end_idx:
        raise mh.MegaError('No overlapping SNPs')

    for ref_pos, ref_es, alt_es, snp_id in zip(chrm_poss[start_idx:end_idx],
                                               chrm_ref_es[start_idx:end_idx],
                                               chrm_alt_es[start_idx:end_idx],
                                               range(start_idx, end_idx)):
        snp_ref_seq = decode_snp_seq(ref_es)
        snp_alt_seq = decode_snp_seq(alt_es)
        if r_ref_pos.strand == 1:
            read_pos = ref_pos - r_ref_pos.start
        else:
            read_pos = r_ref_pos.end - ref_pos - len(snp_ref_seq)
            snp_ref_seq = mh.revcomp(snp_ref_seq)
            snp_alt_seq = mh.revcomp(snp_alt_seq)
        yield read_pos, snp_ref_seq, snp_alt_seq, snp_id, ref_pos

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

    return decode.score_seq(tpost, seq, tpost_start, tpost_end, all_paths)

def call_read_snps(
        r_ref_pos, snps_to_test, edge_buffer, context_bases, r_ref_seq,
        rl_cumsum, r_to_q_poss, r_post, post_mapped_start, all_paths,
        snp_calib_tbl):
    # call all snps overlapping this read
    r_snp_calls = []
    for (r_snp_pos, snp_ref_seq, snp_alt_seq, snp_id,
         snp_ref_pos) in iter_overlapping_snps(
             r_ref_pos, snps_to_test, edge_buffer):
        # select single base SNP or indel context width
        snp_context_bases = (context_bases[0]
                             if len(snp_ref_seq) == len(snp_alt_seq) else
                             context_bases[1])
        pos_bb = min(snp_context_bases, r_snp_pos)
        pos_ab = min(snp_context_bases,
                     r_ref_seq.shape[0] - r_snp_pos - len(snp_ref_seq))
        pos_ref_seq = r_ref_seq[r_snp_pos - pos_bb:
                                r_snp_pos + pos_ab + len(snp_ref_seq)]
        if any(pos_ref_seq[pos_bb:pos_bb + len(snp_ref_seq)] !=
               np.array([mh.ALPHABET.find(b) for b in snp_ref_seq])):
            """print('*'*10+'Refernce seq at {} expected "{}" got "{}"'.format(
                snp_ref_pos,
                ''.join(mh.ALPHABET[b] for b in
                        pos_ref_seq[pos_bb:pos_bb + len(snp_ref_seq)]),
                snp_ref_seq) + '*' * 10)"""
            raise mh.MegaError(
                'Reference SNP sequence does not match reference FASTA.')
        pos_alt_seq = np.concatenate([
            pos_ref_seq[:pos_bb],
            np.array([mh.ALPHABET.find(b) for b in snp_alt_seq],
                     dtype=np.uintp),
            pos_ref_seq[pos_bb + len(snp_ref_seq):]])
        blk_start  = rl_cumsum[r_to_q_poss[r_snp_pos - pos_bb]]
        blk_end = rl_cumsum[r_to_q_poss[r_snp_pos + pos_ab] + 1]

        if blk_end - blk_start < max(len(pos_ref_seq), len(pos_alt_seq)):
            # no valid mapping over large inserted query bases
            # i.e. need as many "events/strides" as bases for valid mapping
            continue
        loc_ref_score = score_seq(
            r_post, pos_ref_seq, post_mapped_start + blk_start,
            post_mapped_start + blk_end, all_paths)
        loc_alt_score = score_seq(
            r_post, pos_alt_seq, post_mapped_start + blk_start,
            post_mapped_start + blk_end, all_paths)
        if np.isnan(loc_ref_score) or np.isnan(loc_alt_score):
            raise mh.MegaError(
                'Score computation error (mapped signal too short for ' +
                'proposed sequence or memory error)')

        fwd_strand_ref_seq = (snp_ref_seq if r_ref_pos.strand == 1 else
                              mh.revcomp(snp_ref_seq))
        fwd_strand_alt_seq = (snp_alt_seq if r_ref_pos.strand == 1 else
                              mh.revcomp(snp_alt_seq))
        if len(fwd_strand_ref_seq) == 0 or len(fwd_strand_alt_seq) == 0:
            # add forward strand reference base to both ref and alt snp seqs
            if r_ref_pos.strand == 1:
                fwd_ref_base = mh.ALPHABET[pos_ref_seq[pos_bb - 1]]
            else:
                fwd_ref_base = mh.comp(mh.ALPHABET[
                    pos_ref_seq[::-1][pos_ab - 1]])
            fwd_strand_ref_seq = fwd_ref_base + fwd_strand_ref_seq
            fwd_strand_alt_seq = fwd_ref_base + fwd_strand_alt_seq
            snp_ref_pos -= 1

        # calibrate llr
        calib_llr = snp_calib_tbl.calibrate_llr(
            loc_ref_score - loc_alt_score, snp_ref_seq, snp_alt_seq)
        r_snp_calls.append((
            snp_ref_pos, calib_llr, fwd_strand_ref_seq, fwd_strand_alt_seq,
            snp_id))

    return r_snp_calls


###############################
##### Per-read SNP Output #####
###############################

def _get_snps_queue(
        snps_q, snps_conn, snp_id_tbl, snps_db_fn, snps_txt_fn, db_safety):
    snps_db = sqlite3.connect(snps_db_fn)
    if db_safety < 2:
        snps_db.execute(SET_ASYNC_MODE)
    if db_safety < 1:
        snps_db.execute(SET_NO_ROLLBACK_MODE)
    snps_db.execute(CREATE_SNPS_TBLS)
    snps_txt_fp = None if snps_txt_fn is None else open(snps_txt_fn, 'w')

    while True:
        try:
            # note strand is +1 for fwd or -1 for rev
            r_snp_calls, (read_id, chrm, strand) = snps_q.get(block=False)
            snps_db.executemany(ADDMANY_SNPS, [
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
        snps_db.executemany(ADDMANY_SNPS, [
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
    snps_db.execute(CREATE_SNPS_IDX)
    snps_db.commit()
    snps_db.close()

    return


######################
##### VCF Reader #####
######################

class SnpCalibrator(object):
    def _load_calibration(self):
        calib_data = np.load(self.fn)
        self.stratify_type = str(calib_data['stratify_type'])
        assert self.stratify_type == 'snp_ins_del'

        self.num_calib_vals = np.int(calib_data['smooth_nvals'])

        self.snp_llr_range = calib_data['snp_llr_range'].copy()
        self.snp_step = (self.snp_llr_range[1] - self.snp_llr_range[0]) / (
            self.num_calib_vals - 1)
        self.snp_calib_table = calib_data['snp_calibration_table'].copy()

        self.del_llr_range = calib_data['deletion_llr_range'].copy()
        self.del_step = (self.del_llr_range[1] - self.del_llr_range[0]) / (
            self.num_calib_vals - 1)
        self.del_calib_table = calib_data['deletion_calibration_table'].copy()

        self.ins_llr_range = calib_data['insertion_llr_range'].copy()
        self.ins_step = (self.ins_llr_range[1] - self.ins_llr_range[0]) / (
            self.num_calib_vals - 1)
        self.ins_calib_table = calib_data['insertion_calibration_table'].copy()

        return

    def __init__(self, snps_calib_fn):
        self.fn = snps_calib_fn
        if self.fn is not None:
            self._load_calibration()
        self.calib_loaded = self.fn is not None
        return

    def calibrate_llr(self, llr, ref_seq, alt_seq):
        if not self.calib_loaded:
            return llr
        if len(ref_seq) == len(alt_seq):
            return self.snp_calib_table[np.around((
                np.clip(llr, self.snp_llr_range[0], self.snp_llr_range[1]) -
                self.snp_llr_range[0]) / self.snp_step).astype(int)]
        elif len(ref_seq) > len(alt_seq):
            return self.del_calib_table[np.around((
                np.clip(llr, self.del_llr_range[0], self.del_llr_range[1]) -
                self.del_llr_range[0]) / self.del_step).astype(int)]
        return self.ins_calib_table[np.around((
            np.clip(llr, self.ins_llr_range[0], self.ins_llr_range[1]) -
            self.ins_llr_range[0]) / self.ins_step).astype(int)]


class SnpData(object):
    def __init__(
            self, snp_fn, do_prepend_chr_vcf, max_snp_size, all_paths,
            write_snps_txt, context_bases, snps_calib_fn=None):
        self.all_paths = all_paths
        self.write_snps_txt = write_snps_txt
        self.calib_table = SnpCalibrator(snps_calib_fn)
        self.context_bases = context_bases
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


######################
##### VCF Writer #####
######################

class Variant(object):
    """ Variant for entry into VcfWriter.
    Currently only handles a single sample.
    """
    def __init__(
            self, chrom, pos, ref, alt, id='.', qual='.', filter='.',
            info=None, sample_dict=None):
        self.chrom = chrom
        self.pos = int(pos)
        self.ref = ref.upper()
        self.alt = alt.upper()
        self.id = str(id)
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
            sorted_keys = ['GT'] + [k for k in sorted_keys if k != 'GT']
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

    def add_haploid_probs(self, probs):
        assert len(probs) == 2, 'Haploid probabilities must be length 2.'
        gt_call = np.argmax(probs)
        # phred scaled likelihoods
        with np.errstate(divide='ignore'):
            raw_pl = -10 * np.log10(probs)
        # "normalized" PL values stored as decsribed by VCF format
        # abs to remove negative 0 from file
        pl = np.abs(np.minimum(raw_pl - raw_pl.min(), mh.MAX_PL_VALUE))
        s_pl = np.sort(pl)

        # add sample tags
        self.qual = '{:.0f}'.format(
            np.around(np.minimum(raw_pl[0], mh.MAX_PL_VALUE)))
        self.add_sample_field('GQ', '{:.0f}'.format(np.around(s_pl[1])))
        self.add_sample_field('PL', '{:.0f},{:.0f}'.format(*np.around(pl)))
        if gt_call == 0:
            # ref
            self.add_sample_field('GT', '0')
        else:
            # alt
            self.add_sample_field('GT', '1')
        return

    def add_diploid_probs(self, probs):
        assert len(probs) == 3, 'Diploid probabilities must be length 3.'
        gt_call = np.argmax(probs)
        # phred scaled likelihoods
        with np.errstate(divide='ignore'):
            raw_pl = -10 * np.log10(probs)
        # "normalized" PL values stored as decsribed by VCF format
        # abs to remove negative 0 from file
        pl = np.abs(np.minimum(raw_pl - raw_pl.min(), mh.MAX_PL_VALUE))
        s_pl = np.sort(pl)

        # add sample tags
        self.qual = '{:.0f}'.format(
            np.around(np.minimum(raw_pl[0], mh.MAX_PL_VALUE)))
        self.add_sample_field('GQ', '{:.0f}'.format(np.around(s_pl[1])))
        self.add_sample_field('PL', '{:.0f},{:.0f},{:.0f}'.format(
            *np.around(pl)))
        if gt_call == 0:
            # homo-ref
            self.add_sample_field('GT', '0/0')
        elif gt_call == 1:
            # het
            self.add_sample_field('GT', '0/1')
        else:
            # homo-alt
            self.add_sample_field('GT', '1/1')
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
    version_options = set(['4.2',])
    def __init__(
            self, filename, mode='w',
            header=('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
                    'INFO', 'FORMAT', 'SAMPLE'),
            extra_meta_info=FIXED_VCF_MI, version='4.2', ref_fn=None):
        self.filename = filename
        self.mode = mode
        self.header = header
        if version not in self.version_options:
            raise ValueError('version must be one of {}'.format(
                self.version_options))
        self.version = version
        self.meta = [
            mh.VCF_VERSION_MI.format(self.version),
            mh.FILE_DATE_MI.format(datetime.date.today().strftime("%Y%m%d")),
            mh.SOURCE_MI.format(MEGALODON_VERSION),
            mh.REF_MI.format(ref_fn)] + extra_meta_info

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

def logsumexp(x):
    x_max = x.max()
    return np.log(np.sum(np.exp(x - x_max))) + x_max

def binom_pmf(k, n, p):
    return (np.math.factorial(n) / (
        np.math.factorial(k) * np.math.factorial(n - k))) * (
            p ** k) * ((1 - p) ** (n - k))

class AggSnps(mh.AbstractAggregationClass):
    """ Class to assist in database queries for per-site aggregation of
    SNP calls over reads.
    """
    def __init__(self, snps_db_fn, write_vcf_llr=False):
        # open as read only database
        self.snps_db = sqlite3.connect(snps_db_fn, uri=True)
        self.n_uniq_snps = None
        self.write_vcf_llr = write_vcf_llr
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

    def compute_diploid_probs(self, llrs, het_factor=1.0):
        if np.errstate(over='ignore'):
            exp_llrs = np.exp(llrs)
        lp_alt = np.sort(np.log(1) - np.log1p(exp_llrs))[::-1]
        lp_ref = np.sort(llrs - np.log1p(exp_llrs))
        lp_hom_alt = np.sum(lp_alt)
        lp_hom_ref = np.sum(lp_ref)
        lp_het = np.log(sum(
            np.exp(np.log(binom_pmf(i, len(llrs), 0.5)) +
                   np.sum(lp_alt[:i]) + np.sum(lp_ref[i:]))
            for i in range(len(llrs) + 1)))
        prior_weights = np.array([1.0, het_factor ** len(llrs), 1.0])
        prior_weights /= prior_weights.sum()
        snp_lps = np.array([lp_hom_ref, lp_het, lp_hom_alt]) + np.log(
            prior_weights)
        post_snp_lps = snp_lps - logsumexp(snp_lps)
        return np.exp(post_snp_lps)

    def compute_haploid_probs(self, llrs):
        if np.errstate(over='ignore'):
            exp_llrs = np.exp(llrs)
        lp_ref = np.sum(llrs - np.log1p(exp_llrs))
        lp_alt = np.sum(np.log(1) - np.log1p(exp_llrs))
        snp_lps = np.array([lp_ref, lp_alt])
        post_snp_lps = snp_lps - logsumexp(snp_lps)
        return np.exp(post_snp_lps)

    def compute_snp_stats(self, snp_loc, het_factors, call_mode=DIPLOID_MODE):
        assert call_mode in (HAPLIOD_MODE, DIPLOID_MODE), (
            'Invalid SNP aggregation call mode.')

        pr_snp_stats = self.get_per_read_snp_stats(snp_loc)
        llrs = np.array([r_stats.score for r_stats in pr_snp_stats])
        if AMBIG_LLRS_THRESH is not None:
            llrs = llrs[np.logical_or(llrs < -AMBIG_LLRS_THRESH,
                                      llrs > AMBIG_LLRS_THRESH)]
            if len(llrs) == 0:
                return None
        r0_stats = pr_snp_stats[0]
        snp_var = Variant(
            chrom=r0_stats.chrm, pos=r0_stats.pos, ref=r0_stats.ref_seq,
            alt=r0_stats.alt_seq, id=r0_stats.snp_id)
        snp_var.add_tag('DP', '{}'.format(len(pr_snp_stats)))
        snp_var.add_sample_field('DP', '{}'.format(len(pr_snp_stats)))

        if self.write_vcf_llr:
            snp_var.add_sample_field('LLRS', ','.join(
                '{:.2f}'.format(llr) for llr in  sorted(llrs)))

        if call_mode == DIPLOID_MODE:
            het_factor = (
                het_factors[0] if len(r0_stats.ref_seq) == 1 and
                len(r0_stats.alt_seq) == 1 else
                het_factors[1])
            diploid_probs = self.compute_diploid_probs(llrs, het_factor)
            snp_var.add_diploid_probs(diploid_probs)
        elif call_mode == HAPLIOD_MODE:
            haploid_probs = self.compute_haploid_probs(llrs)
            snp_var.add_haploid_probs(haploid_probs)

        return snp_var

    def close(self):
        self.snps_db.close()
        return


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
