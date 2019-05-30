import sys
import queue
import sqlite3
import datetime
from time import sleep
import multiprocessing as mp
from collections import defaultdict, namedtuple, OrderedDict

import pysam
import numpy as np

from megalodon import (calibration, decode, logging, mapping,
                       megalodon_helper as mh)
from megalodon._version import MEGALODON_VERSION


DEBUG = False

DIPLOID_MODE = 'diploid'
HAPLIOD_MODE = 'haploid'

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
SELECT DISTINCT chrm, pos, snp_id, ref_seq FROM snps)"""
SEL_UNIQ_SNP_ID = '''
SELECT DISTINCT chrm, pos, snp_id, ref_seq FROM snps'''
SEL_SNP_STATS = '''
SELECT * FROM snps WHERE chrm=? AND pos=? AND
snp_id=? AND ref_seq=?'''

FIXED_VCF_MI = [
    'phasing=none',
    'INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
    'FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    'FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
    'FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
    'FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">'
]


################################
##### Per-read SNP Scoring #####
################################

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

    score = decode.score_seq(tpost, seq, tpost_start, tpost_end, all_paths)
    if np.isnan(score):
        raise mh.MegaError(
            'Score computation error (mapped signal too short for ' +
            'proposed sequence or memory error)')

    return score

def call_read_snps(
        snps_data, r_ref_pos, edge_buffer, r_ref_seq, rl_cumsum, r_to_q_poss,
        r_post, post_mapped_start):
    # call all snps overlapping this read
    r_snp_calls = []
    for (snp_ref_seq, snp_alt_seqs, snp_id,
         snp_ref_pos) in snps_data.iter_overlapping_snps(
             r_ref_pos, edge_buffer):

        if r_ref_pos.strand == 1:
            read_pos = snp_ref_pos - 1 - r_ref_pos.start
            read_ref_seq = snp_ref_seq
            read_alt_seqs = snp_alt_seqs
        else:
            read_pos = r_ref_pos.end - snp_ref_pos + 1 - len(snp_ref_seq)
            read_ref_seq = mh.revcomp(snp_ref_seq)
            read_alt_seqs = [mh.revcomp(alt_seq) for alt_seq in snp_alt_seqs]

        # select single base SNP or indel context width
        is_snp = all(len(snp_ref_seq) == len(snp_alt_seq)
                     for snp_alt_seq in snp_alt_seqs)
        is_del = not is_snp and len(snp_ref_seq) > len(snp_alt_seqs[0])
        snp_context_bases = snps_data.indel_context if is_snp else \
                            snps_data.snp_context
        pos_bb = min(snp_context_bases, read_pos)
        pos_ab = min(snp_context_bases,
                     r_ref_seq.shape[0] - read_pos - len(read_ref_seq))
        pos_ref_seq = r_ref_seq[read_pos - pos_bb:
                                read_pos + pos_ab + len(read_ref_seq)]
        if DEBUG and any(
                pos_ref_seq[pos_bb:pos_bb + len(snp_ref_seq)] !=
                np.array([mh.ALPHABET.find(b) for b in read_ref_seq])):
            logger = logging.get_logger()
            logger.debug(
                '*'*10 + 'Refernce seq at {} expected {}[{}]{} got "{}"'.format(
                    snp_ref_pos,
                    ''.join(mh.ALPHABET[b] for b in
                            pos_ref_seq[pos_bb-3:pos_bb]),
                    ''.join(mh.ALPHABET[b] for b in
                            pos_ref_seq[pos_bb:pos_bb + len(snp_ref_seq)]),
                    ''.join(mh.ALPHABET[b] for b in
                            pos_ref_seq[pos_bb + len(snp_ref_seq):
                                        pos_bb + len(snp_ref_seq) + 3]),
                    read_ref_seq, ) + '*' * 10)
            raise mh.MegaError(
                'Reference SNP sequence does not match reference FASTA.')
        blk_start  = rl_cumsum[r_to_q_poss[read_pos - pos_bb]]
        blk_end = rl_cumsum[r_to_q_poss[read_pos + pos_ab] + 1]
        if blk_end - blk_start < max(len(pos_ref_seq), max(
                len(read_alt_seq) for read_alt_seq in read_alt_seqs)):
            # no valid mapping over large inserted query bases
            # i.e. need as many "events/strides" as bases for valid mapping
            continue

        loc_ref_score = score_seq(
            r_post, pos_ref_seq, post_mapped_start + blk_start,
            post_mapped_start + blk_end, snps_data.all_paths)
        loc_alt_llrs = []
        for read_alt_seq in read_alt_seqs:
            pos_alt_seq = np.concatenate([
                pos_ref_seq[:pos_bb],
                np.array([mh.ALPHABET.find(b) for b in read_alt_seq],
                         dtype=np.uintp),
                pos_ref_seq[pos_bb + len(snp_ref_seq):]])
            loc_alt_score = score_seq(
                r_post, pos_alt_seq, post_mapped_start + blk_start,
                post_mapped_start + blk_end, snps_data.all_paths)
            # calibrate log probs
            loc_alt_llrs.append(snps_data.calibrate_llr(
                loc_ref_score - loc_alt_score, snp_ref_seq, pos_alt_seq))

        # due to calibration mutli-allelic log likelihoods could result in
        # inferred negative reference likelihood, so re-normalize her
        loc_alt_log_ps = calibration.compute_alt_log_probs(
            np.array(loc_alt_llrs))

        r_snp_calls.append((
            snp_ref_pos, loc_alt_log_ps, snp_ref_seq, snp_alt_seqs, snp_id))

    return r_snp_calls


###############################
##### Per-read SNP Output #####
###############################

def annotate_snps(r_start, ref_seq, r_snp_calls, strand):
    """ Annotate reference sequence with called snps.

    Note: Reference sequence is in read orientation and snp calls are in
    genome coordiates.
    """
    snp_seqs = []
    prev_pos = 0
    if strand == -1:
        ref_seq = ref_seq[::-1]
    for snp_pos, alt_lps, snp_ref_seq, snp_alt_seqs, _ in sorted(r_snp_calls):
        ref_lp = np.log1p(-np.exp(alt_lps).sum())
        # called canonical
        if ref_lp >= min(alt_lps): continue
        alt_seq = snp_alt_seqs[np.argmax(alt_lps)]
        if strand == -1:
            alt_seq = mh.revcomp(alt_seq)
        snp_seqs.append(ref_seq[prev_pos:snp_pos - r_start] + alt_seq)
        prev_pos = snp_pos - r_start + len(snp_ref_seq)
    snp_seqs.append(ref_seq[prev_pos:])
    snp_seq = ''.join(snp_seqs)
    if strand == -1:
        snp_seq = snp_seq[::-1]

    return snp_seq

def _get_snps_queue(
        snps_q, snps_conn, snps_db_fn, snps_txt_fn, db_safety,
        pr_refs_fn, pr_ref_filts):
    def get_snp_call():
        # note strand is +1 for fwd or -1 for rev
        r_snp_calls, (read_id, chrm, strand, r_start, ref_seq, read_len,
                      q_st, q_en, cigar) = snps_q.get(block=False)
        snps_db.executemany(ADDMANY_SNPS, [
            (read_id, chrm, strand, pos, alt_lp,
             snp_ref_seq, snp_alt_seq, snp_id)
            for pos, alt_lps, snp_ref_seq, snp_alt_seqs, snp_id in r_snp_calls
            for alt_lp, snp_alt_seq in zip(alt_lps, snp_alt_seqs)])
        if snps_txt_fp is not None:
            # would involve batching and creating several conversion tables
            # for var strings (read_if and chrms).
            snps_txt_fp.write('\n'.join((
                '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                    read_id, chrm, strand, pos, score,
                    snp_ref_seq, snp_alt_seq, snp_id)
                for pos, score, snp_ref_seq, snp_alt_seq, snp_id in
                r_snp_calls)) + '\n')
            snps_txt_fp.flush()
        if pr_refs_fn is not None:
            if not mapping.read_passes_filters(
                    pr_ref_filts, read_len, q_st, q_en, cigar):
                return
            pr_refs_fp.write('>{}\n{}\n'.format(read_id, annotate_snps(
                r_start, ref_seq, r_snp_calls, strand)))
            pr_refs_fp.flush()
        return


    snps_db = sqlite3.connect(snps_db_fn)
    if db_safety < 2:
        snps_db.execute(SET_ASYNC_MODE)
    if db_safety < 1:
        snps_db.execute(SET_NO_ROLLBACK_MODE)
    snps_db.execute(CREATE_SNPS_TBLS)
    snps_txt_fp = None if snps_txt_fn is None else open(snps_txt_fn, 'w')

    if pr_refs_fn is not None:
        pr_refs_fp = open(pr_refs_fn, 'w')

    while True:
        try:
            get_snp_call()
        except queue.Empty:
            if snps_conn.poll():
                break
            sleep(0.1)
            continue

    while not snps_q.empty():
        get_snp_call()
    if snps_txt_fp is not None: snps_txt_fp.close()
    if pr_refs_fn is not None: pr_refs_fp.close()
    snps_db.execute(CREATE_SNPS_IDX)
    snps_db.commit()
    snps_db.close()

    return


######################
##### VCF Reader #####
######################

class SnpData(object):
    def __init__(
            self, variant_fn, max_indel_size, all_paths,
            write_snps_txt, context_bases, snps_calib_fn=None,
            call_mode=DIPLOID_MODE, do_pr_ref_snps=False, aligner=None,
            keep_snp_fp_open=False):
        logger = logging.get_logger('snps')
        self.max_indel_size = max_indel_size
        self.all_paths = all_paths
        self.write_snps_txt = write_snps_txt
        self.snps_calib_fn = snps_calib_fn
        self.calib_table = calibration.SnpCalibrator(self.snps_calib_fn)
        self.context_bases = context_bases
        if len(self.context_bases) != 2:
            raise mh.MegaError(
                'Must provide 2 context bases values (for single base SNPs ' +
                'and indels).')
        self.call_mode = call_mode
        self.do_pr_ref_snps = do_pr_ref_snps
        self.variant_fn = variant_fn
        self.variants_idx = None
        if self.variant_fn is None:
            return

        logger.info('Loading variants.')
        vars_idx = pysam.VariantFile(variant_fn)
        try:
            contigs = list(vars_idx.header.contigs.keys())
            vars_idx.fetch(next(iter(contigs)), 0, 0)
        except ValueError:
            raise mh.MegaError(
                'Variants file must be indexed. Use bgzip and tabix.')
        if keep_snp_fp_open:
            self.variants_idx = vars_idx
        else:
            vars_idx.close()
            self.variants_idx = None

        if aligner is None:
            raise mh.MegaError(
                'Must provide aligner if SNP filename is provided')
        if len(set(aligner.ref_names_and_lens[0]).intersection(contigs)) == 0:
            raise mh.MegaError((
                'Reference and variant files contain no chromosomes/contigs ' +
                'in common.\n\t\tFirst 3 reference contigs:\t{}\n\t\tFirst 3 ' +
                'variant file contigs:\t{}').format(
                    ', '.join(aligner.ref_names_and_lens[0][:3]),
                    ', '.join(contigs[:3])))

        return

    @property
    def snp_context(self):
        return self.context_bases[0]
    @property
    def indel_context(self):
        return self.context_bases[1]

    def reopen_variant_index(self):
        self.variants_idx = pysam.VariantFile(self.variant_fn)
        return

    def iter_overlapping_snps(self, r_ref_pos, edge_buffer):
        """Iterator over SNPs overlapping the read mapped position.

        SNPs within edge buffer of the end of the mapping will be ignored.
        """
        if r_ref_pos.end - r_ref_pos.start <= 2 * edge_buffer:
            mh.MegaError('Mapped region too short for SNP calling.')
        for variant in self.variants_idx.fetch(
                r_ref_pos.chrm, r_ref_pos.start + edge_buffer,
                r_ref_pos.end - edge_buffer):
            snp_ref_seq = variant.ref
            snp_alt_seqs = variant.alts
            # skip SNPs larger than specified limit
            if self.max_indel_size is not None and max(
                    np.abs(len(snp_ref_seq) - len(snp_alt_seq))
                    for snp_alt_seq in snp_alt_seqs) > self.max_indel_size:
                continue
            yield snp_ref_seq, snp_alt_seqs, variant.id, variant.pos

        return

    def calibrate_llr(self, llr, snp_ref_seq, snp_alt_seq):
        return self.calib_table.calibrate_llr(llr, snp_ref_seq, snp_alt_seq)


######################
##### VCF Writer #####
######################

class Variant(object):
    """ Variant for entry into VcfWriter.
    Currently only handles a single sample.
    """
    def __init__(
            self, chrom, pos, ref, alts, id='.', qual='.', filter='.',
            info=None, sample_dict=None):
        self.chrom = chrom
        self.pos = int(pos)
        self.ref = ref.upper()
        self.alts = [alt.upper() for alt in alts]
        self.alt = ','.join(self.alts)
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
        # phred scaled likelihoods
        with np.errstate(divide='ignore'):
            raw_pl = -10 * np.log10(probs)
        # "normalized" PL values stored as decsribed by VCF format
        # abs to remove negative 0 from file
        pl = np.abs(np.minimum(raw_pl - raw_pl.min(), mh.MAX_PL_VALUE))
        s_pl = np.sort(pl)

        # add sample tags
        self.add_sample_field('GT', gts[np.argmax(probs)])
        self.qual = '{:.0f}'.format(
            np.around(np.minimum(raw_pl[0], mh.MAX_PL_VALUE)))
        self.add_sample_field('GQ', '{:.0f}'.format(np.around(s_pl[1])))
        self.add_sample_field(
            'PL', ','.join('{:.0f}' for _ in range(probs.shape[0])).format(
                *np.around(pl)))
        return

    def add_diploid_probs(self, probs, gts):
        # phred scaled likelihoods
        with np.errstate(divide='ignore'):
            raw_pl = -10 * np.log10(probs)
        # "normalized" PL values stored as decsribed by VCF format
        # abs to remove negative 0 from file
        pl = np.abs(np.minimum(raw_pl - raw_pl.min(), mh.MAX_PL_VALUE))
        s_pl = np.sort(pl)

        # add sample tags
        self.add_sample_field('GT', gts[np.argmax(probs)])
        self.qual = '{:.0f}'.format(
            np.around(np.minimum(raw_pl[0], mh.MAX_PL_VALUE)))
        self.add_sample_field('GQ', '{:.0f}'.format(np.around(s_pl[1])))
        self.add_sample_field(
            'PL', ','.join('{:.0f}' for _ in range(probs.shape[0])).format(
                *np.around(pl)))
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
            extra_meta_info=FIXED_VCF_MI, version='4.2', ref_fn=None,
            ref_names_and_lens=None):
        self.filename = filename
        self.mode = mode
        self.header = header
        if version not in self.version_options:
            raise ValueError('version must be one of {}'.format(
                self.version_options))
        self.version = version
        contig_mis = [] if ref_names_and_lens is None else [
            mh.CONTIG_MI.format(ref_name, ref_len)
            for ref_name, ref_len in zip(*ref_names_and_lens)]
        self.meta = [
            mh.VCF_VERSION_MI.format(self.version),
            mh.FILE_DATE_MI.format(datetime.date.today().strftime("%Y%m%d")),
            mh.SOURCE_MI.format(MEGALODON_VERSION),
            mh.REF_MI.format(ref_fn)] + contig_mis + extra_meta_info
        # TODO add meta info for LOG_PROBS to make field valid VCF

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
    def __init__(self, snps_db_fn, write_vcf_log_probs=False):
        # open as read only database
        self.snps_db = sqlite3.connect(snps_db_fn, uri=True)
        self.n_uniq_snps = None
        self.write_vcf_log_probs = write_vcf_log_probs
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

    def compute_diploid_probs(self, ref_lps, alts_lps, het_factor=1.0):
        def compute_het_lp(a1, a2):
            a1_lps = all_lps[a1]
            a2_lps = all_lps[a2]
            lp_ord = np.argsort(a1_lps - a2_lps)
            s_a1_lps = a1_lps[lp_ord]
            s_a2_lps = a2_lps[lp_ord]
            return np.log(np.sum(
                np.exp(np.log(binom_pmf(i, all_lps.shape[1], 0.5)) +
                       np.sum(s_a1_lps[:i]) + np.sum(s_a2_lps[i:]))
                for i in range(all_lps.shape[1] + 1)))


        all_lps = np.concatenate([ref_lps.reshape(1, -1), alts_lps], axis=0)
        genotype_lps, het_gts, gts = [], [], []
        for a2 in range(all_lps.shape[0]):
            for a1 in range(a2 + 1):
                gts.append('{}/{}'.format(a1, a2))
                if a1 == a2:
                    genotype_lps.append(np.sum(all_lps[a1]))
                    het_gts.append(False)
                else:
                    genotype_lps.append(compute_het_lp(a1, a2))
                    het_gts.append(True)

        prior_weights = np.array([
            1.0 if het_gt else het_factor ** all_lps.shape[1]
            for het_gt in het_gts])
        prior_weights /= prior_weights.sum()
        snp_lps = np.array(genotype_lps) + np.log(prior_weights)
        post_snp_lps = snp_lps - logsumexp(snp_lps)
        return np.exp(post_snp_lps), gts

    def compute_haploid_probs(self, ref_lps, alts_lps):
        snp_lps = np.concatenate([[ref_lps.sum()], alts_lps.sum(axis=0)])
        post_snp_lps = snp_lps - logsumexp(snp_lps)
        return np.exp(post_snp_lps), list(map(str, range(snp_lps.shape[0])))

    def compute_snp_stats(self, snp_loc, het_factors, call_mode=DIPLOID_MODE):
        assert call_mode in (HAPLIOD_MODE, DIPLOID_MODE), (
            'Invalid SNP aggregation ploidy call mode: {}.'.format(call_mode))

        pr_snp_stats = self.get_per_read_snp_stats(snp_loc)
        alt_seqs = sorted(set(r_stats.alt_seq for r_stats in pr_snp_stats))
        pr_alt_lps = defaultdict(dict)
        for r_stats in pr_snp_stats:
            pr_alt_lps[r_stats.read_id][r_stats.alt_seq] = r_stats.score
        alt_seq_lps = [[] for _ in range(len(alt_seqs))]
        for read_lps in pr_alt_lps.values():
            for i, alt_seq in enumerate(alt_seqs):
                try:
                    alt_seq_lps[i].append(read_lps[alt_seq])
                except KeyError:
                    mh.MegaError(
                        'Alternative SNP seqence must exist for all reads.')
        alts_lps = np.stack(alt_seq_lps, axis=0)
        with np.errstate(all='ignore'):
            ref_lps = np.log1p(-np.exp(alts_lps).sum(axis=0))

        r0_stats = pr_snp_stats[0]
        snp_var = Variant(
            chrom=r0_stats.chrm, pos=r0_stats.pos, ref=r0_stats.ref_seq,
            alts=alt_seqs, id=r0_stats.snp_id)
        snp_var.add_tag('DP', '{}'.format(ref_lps.shape[0]))
        snp_var.add_sample_field('DP', '{}'.format(ref_lps.shape[0]))

        if self.write_vcf_log_probs:
            snp_var.add_sample_field('LOG_PROBS', ';'.join(
                ','.join('{:.2f}'.format(lp) for lp in alt_i_lps)
                for alt_i_lps in alts_lps))

        if call_mode == DIPLOID_MODE:
            het_factor = (
                het_factors[0] if len(r0_stats.ref_seq) == 1 and
                len(r0_stats.alt_seq) == 1 else
                het_factors[1])
            diploid_probs, gts = self.compute_diploid_probs(
                ref_lps, alts_lps, het_factor)
            snp_var.add_diploid_probs(diploid_probs, gts)
        elif call_mode == HAPLIOD_MODE:
            haploid_probs, gts = self.compute_haploid_probs(ref_lps, alts_lps)
            snp_var.add_haploid_probs(haploid_probs, gts)

        return snp_var

    def close(self):
        self.snps_db.close()
        return


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
