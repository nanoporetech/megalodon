import os
import sys
import queue
import sqlite3
import datetime
from time import sleep
from array import array
import multiprocessing as mp
from operator import itemgetter
from collections import defaultdict, namedtuple, OrderedDict

import pysam
import numpy as np

from megalodon import (calibration, decode, logging, mapping,
                       megalodon_helper as mh)
from megalodon._version import MEGALODON_VERSION


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
SELECT * FROM snps WHERE chrm IS ? AND pos IS ? AND
snp_id IS ? AND ref_seq IS ?'''

SAMPLE_NAME = 'SAMPLE'
# specified by sam format spec
WHATSHAP_MAX_QUAL = 40
WHATSHAP_RG_ID = '1'
FIXED_VCF_MI = [
    'phasing=none',
    'INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
    'FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    'FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
    'FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
    ('FORMAT=<ID=GL,Number=G,Type=Float,' +
     'Description="Log10 likelihoods for genotypes">'),
    ('FORMAT=<ID=PL,Number=G,Type=Integer,' +
     'Description="Normalized, Phred-scaled likelihoods for genotypes">')
]
FORMAT_LOG_PROB_MI = (
    'FORMAT=<ID=LOG_PROBS,Number=A,Type=String,' +
    'Description="Per-read log10 likelihoods for alternative ' +
    'alleles (semi-colon separated)">')


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
        r_post, post_mapped_start, snp_merge_gap=2):
    # call all snps overlapping this read
    r_snp_calls = []
    for (snp_ref_seq, snp_alt_seqs, snp_id,
         snp_ref_pos) in snps_data.iter_nearby_snp_grps(
             r_ref_pos, edge_buffer, snp_merge_gap):
        # TODO test all combinations of snps in the group and then
        # take marginal probabilities for each individual SNP

        if r_ref_pos.strand == 1:
            read_pos = snp_ref_pos - r_ref_pos.start
            read_ref_seq = snp_ref_seq
            read_alt_seqs = snp_alt_seqs
        else:
            read_pos = r_ref_pos.end - snp_ref_pos - len(snp_ref_seq)
            read_ref_seq = mh.revcomp(snp_ref_seq)
            read_alt_seqs = [mh.revcomp(alt_seq) for alt_seq in snp_alt_seqs]

        # select single base SNP or indel context width
        snp_context_bases = snps_data.indel_context if all(
            len(snp_ref_seq) == len(snp_alt_seq)
            for snp_alt_seq in snp_alt_seqs) else snps_data.snp_context
        pos_bb = min(snp_context_bases, read_pos)
        pos_ab = min(snp_context_bases,
                     r_ref_seq.shape[0] - read_pos - len(read_ref_seq))
        pos_ref_seq = r_ref_seq[read_pos - pos_bb:
                                read_pos + pos_ab + len(read_ref_seq)]
        # TODO move this to an initial check of a small number of variants
        # against the reference
        if any(pos_ref_seq[pos_bb:pos_bb + len(snp_ref_seq)] !=
               np.array([mh.ALPHABET.find(b) for b in read_ref_seq])):
            # variant reference sequence does not match fasta reference
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
            continue
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
                loc_ref_score - loc_alt_score, read_ref_seq, read_alt_seq))

        # due to calibration mutli-allelic log likelihoods could result in
        # inferred negative reference likelihood, so re-normalize here
        loc_alt_log_ps = calibration.compute_log_probs(np.array(loc_alt_llrs))

        r_snp_calls.append((
            snp_ref_pos, loc_alt_log_ps, snp_ref_seq, snp_alt_seqs, snp_id))

    return r_snp_calls


###############################
##### Per-read SNP Output #####
###############################

def log_prob_to_phred(log_prob):
    return -10 * np.log10(1 - np.exp(log_prob))

def simplify_snp_seq(ref_seq, alt_seq):
    trim_before = trim_after = 0
    while (len(ref_seq) > 0 and len(alt_seq) > 0 and
           ref_seq[0] == alt_seq[0]):
        trim_before += 1
        ref_seq = ref_seq[1:]
        alt_seq = alt_seq[1:]
    while (len(ref_seq) > 0 and len(alt_seq) > 0 and
           ref_seq[-1] == alt_seq[-1]):
        trim_after += 1
        ref_seq = ref_seq[:-1]
        alt_seq = alt_seq[:-1]

    return ref_seq, alt_seq, trim_before, trim_after

def iter_non_overlapping_snps(snp_calls):
    def get_max_prob_allele_snp(snp_grp):
        """ For overlapping SNPs return the snp with the highest probability
        single allele as this one will be added to the reference sequence.

        More complex chained SNPs could be handled, but are not here.
        For example, a 5 base deletion covering 2 single base swap SNPs could
        validly result in 2 alternative single base swap alleles, but the
        logic here would only allow one of those alternatives since they
        are covered by the same reference deletion. There are certainly many
        more edge cases than this and each one would require specific logic.
        This likely covers the majority of valid cases and limiting to
        50 base indels by default limits the scope of this issue.
        """
        most_prob_snp = None
        for snp_data in snp_grp:
            ref_lp = np.log1p(-np.exp(snp_data[1]).sum())
            snp_max_lp = max(ref_lp, snp_data[1].max())
            if most_prob_snp is None or snp_max_lp > most_prob_snp[0]:
                most_prob_snp = (snp_max_lp, ref_lp, snp_data)

        _, ref_lp, (snp_pos, alt_lps, snp_ref_seq,
                    snp_alt_seqs, _) = most_prob_snp
        return snp_pos, alt_lps, snp_ref_seq, snp_alt_seqs, ref_lp


    snp_calls_iter = iter(snp_calls)
    # initialize snp_grp with first snp
    snp_data = next(snp_calls_iter)
    prev_snp_end = snp_data[0] + len(snp_data[2])
    snp_grp = [snp_data]
    for snp_data in sorted(snp_calls, key=itemgetter(0)):
        if snp_data[0] < prev_snp_end:
            prev_snp_end = max(snp_data[0] + len(snp_data[2]), prev_snp_end)
            snp_grp.append(snp_data)
        else:
            yield get_max_prob_allele_snp(snp_grp)
            prev_snp_end = snp_data[0] + len(snp_data[2])
            snp_grp = [snp_data]

    # yeild last snp grp data
    yield get_max_prob_allele_snp(snp_grp)
    return

def annotate_snps(r_start, ref_seq, r_snp_calls, strand):
    """ Annotate reference sequence with called snps.

    Note: Reference sequence is in read orientation and snp calls are in
    genome coordiates.
    """
    snp_seqs, snp_quals, snp_cigar = [], [], []
    prev_pos, curr_match = 0, 0
    # ref_seq is read-centric so flop order to process snps in genomic order
    if strand == -1:
        ref_seq = ref_seq[::-1]
    for (snp_pos, alt_lps, snp_ref_seq, snp_alt_seqs,
         ref_lp) in iter_non_overlapping_snps(r_snp_calls):
        prev_len = snp_pos - r_start - prev_pos
        # called canonical
        if ref_lp >= max(alt_lps):
            snp_seqs.append(
                ref_seq[prev_pos:snp_pos - r_start + len(snp_ref_seq)])
            snp_quals.extend(
                ([WHATSHAP_MAX_QUAL] * prev_len) +
                ([min(log_prob_to_phred(ref_lp), WHATSHAP_MAX_QUAL)] *
                 len(snp_ref_seq)))
            curr_match += prev_len + len(snp_ref_seq)
        else:
            alt_seq = snp_alt_seqs[np.argmax(alt_lps)]
            # complement since ref_seq is complement seq
            # (not reversed; see loop init)
            read_alt_seq = alt_seq if strand == 1 else mh.comp(alt_seq)
            snp_seqs.append(ref_seq[prev_pos:snp_pos - r_start] + read_alt_seq)
            snp_quals.extend(
                ([WHATSHAP_MAX_QUAL] * prev_len) +
                ([min(log_prob_to_phred(max(alt_lps)), WHATSHAP_MAX_QUAL)] *
                 len(alt_seq)))

            # add cigar information for snp or indel
            t_ref_seq, t_alt_seq, t_before, t_after = simplify_snp_seq(
                snp_ref_seq, alt_seq)
            curr_match += t_before
            snp_cigar.append((7, curr_match + prev_len))
            if len(t_alt_seq) == len(t_ref_seq):
                snp_cigar.append((8, len(t_alt_seq)))
            elif len(t_alt_seq) > len(t_ref_seq):
                # left justify mismatch bases in complex insertion
                if len(t_ref_seq) != 0:
                    snp_cigar.append((8, len(t_ref_seq)))
                snp_cigar.append((1, len(t_alt_seq) - len(t_ref_seq)))
            else:
                # left justify mismatch bases in complex deletion
                if len(t_alt_seq) != 0:
                    snp_cigar.append((8, len(t_alt_seq)))
                snp_cigar.append((2, len(t_ref_seq) - len(t_alt_seq)))
            curr_match = t_after
        prev_pos = snp_pos - r_start + len(snp_ref_seq)

    snp_seqs.append(ref_seq[prev_pos:])
    snp_seq = ''.join(snp_seqs)
    if strand == -1:
        snp_seq = snp_seq[::-1]
    len_remain = len(ref_seq) - prev_pos
    snp_quals.extend([WHATSHAP_MAX_QUAL] * len_remain)
    if strand == -1:
        snp_quals = snp_quals[::-1]
    snp_quals = list(map(int, snp_quals))
    snp_cigar.append((7, len_remain + curr_match))
    if strand == -1:
        snp_cigar = snp_cigar[::-1]

    return snp_seq, snp_quals, snp_cigar

def _get_snps_queue(
        snps_q, snps_conn, snps_db_fn, snps_txt_fn, db_safety, pr_refs_fn,
        pr_ref_filts, whatshap_map_fn, ref_names_and_lens, ref_fn):
    def write_whatshap_alignment(
            read_id, snp_seq, snp_quals, chrm, strand, r_st, snp_cigar):
        a = pysam.AlignedSegment()
        a.query_name = read_id
        a.flag = 0 if strand == 1 else 16
        a.reference_id = whatshap_map_fp.get_tid(chrm)
        a.reference_start = r_st
        a.template_length = len(snp_seq)
        a.mapping_quality = WHATSHAP_MAX_QUAL
        a.set_tags([('RG', WHATSHAP_RG_ID)])

        # convert to reference based sequence
        if strand == -1:
            snp_seq = mh.revcomp(snp_seq)
            snp_quals = snp_quals[::-1]
            snp_cigar = snp_cigar[::-1]
        a.query_sequence = snp_seq
        a.query_qualities = array('B', snp_quals)
        a.cigartuples = snp_cigar
        whatshap_map_fp.write(a)

        return

    def get_snp_call():
        # note strand is +1 for fwd or -1 for rev
        r_snp_calls, (read_id, chrm, strand, r_start, ref_seq, read_len,
                      q_st, q_en, cigar) = snps_q.get(block=False)
        snps_db.executemany(ADDMANY_SNPS, [
            (read_id, chrm, strand, pos, alt_lp,
             snp_ref_seq, snp_alt_seq, snp_id)
            for pos, alt_lps, snp_ref_seq, snp_alt_seqs, snp_id in r_snp_calls
            for alt_lp, snp_alt_seq in zip(alt_lps, snp_alt_seqs)])
        if snps_txt_fp is not None and len(r_snp_calls) > 0:
            snps_txt_fp.write('\n'.join((
                ('\t'.join('{}' for _ in field_names)).format(
                    read_id, chrm, strand, pos,
                    np.log1p(-np.exp(alt_lps).sum()), alt_lp,
                    snp_ref_seq, snp_alt_seq, snp_id)
                for pos, alt_lps, snp_ref_seq, snp_alt_seqs, snp_id
                in r_snp_calls
                for alt_lp, snp_alt_seq in zip(alt_lps, snp_alt_seqs))) + '\n')
            snps_txt_fp.flush()
        if do_ann_snps:
            if not mapping.read_passes_filters(
                    pr_ref_filts, read_len, q_st, q_en, cigar):
                return
            snp_seq, snp_quals, snp_cigar = annotate_snps(
                r_start, ref_seq, r_snp_calls, strand)
            if pr_refs_fn is not None:
                pr_refs_fp.write('>{}\n{}\n'.format(read_id, snp_seq))
                pr_refs_fp.flush()
            if whatshap_map_fn is not None:
                write_whatshap_alignment(
                    read_id, snp_seq, snp_quals, chrm, strand, r_start,
                    snp_cigar)

        return


    snps_db = sqlite3.connect(snps_db_fn)
    if db_safety < 2:
        snps_db.execute(SET_ASYNC_MODE)
    if db_safety < 1:
        snps_db.execute(SET_NO_ROLLBACK_MODE)
    snps_db.execute(CREATE_SNPS_TBLS)
    if snps_txt_fn is None:
        snps_txt_fp = None
    else:
        snps_txt_fp = open(snps_txt_fn, 'w')
        field_names = (
            'read_id', 'chrm', 'strand', 'pos', 'ref_log_prob', 'alt_log_prob',
            'ref_seq', 'alt_seq', 'snp_id')
        snps_txt_fp.write('\t'.join(field_names) + '\n')

    if pr_refs_fn is not None:
        pr_refs_fp = open(pr_refs_fn, 'w')

    if whatshap_map_fn is not None:
        _, map_fmt = os.path.splitext(whatshap_map_fn)
        if map_fmt == '.bam': w_mode = 'wb'
        elif map_fmt == '.cram': w_mode = 'wc'
        elif map_fmt == '.sam': w_mode = 'w'
        else:
            raise mh.MegaError('Invalid mapping output format')
        header = {
            'HD': {'VN': '1.4'},
            'SQ': [{'LN': ref_len, 'SN': ref_name}
                   for ref_name, ref_len in sorted(zip(*ref_names_and_lens))],
            'RG': [{'ID':WHATSHAP_RG_ID, 'SM':SAMPLE_NAME},]}
        whatshap_map_fp = pysam.AlignmentFile(
            whatshap_map_fn, w_mode, header=header, reference_filename=ref_fn)

    do_ann_snps = whatshap_map_fn is not None or pr_refs_fn is not None

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
    if whatshap_map_fn is not None: whatshap_map_fp.close()
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
        vars_idx = pysam.VariantFile(self.variant_fn)
        try:
            contigs = list(vars_idx.header.contigs.keys())
            vars_idx.fetch(next(iter(contigs)), 0, 0)
        except ValueError:
            logger.warn(
                'Variants file must be indexed. Performing indexing now.')
            vars_idx.close()
            self.variant_fn = index_variants(self.variant_fn)
            vars_idx = pysam.VariantFile(self.variant_fn)
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
        if self.variant_fn is not None:
            self.variants_idx = pysam.VariantFile(self.variant_fn)
        return

    def iter_nearby_snp_grps(self, r_ref_pos, edge_buffer, snp_merge_gap):
        """Iterator over SNPs overlapping the read mapped position.

        SNPs within edge buffer of the end of the mapping will be ignored.
        """
        if r_ref_pos.end - r_ref_pos.start <= 2 * edge_buffer:
            raise mh.MegaError('Mapped region too short for SNP calling.')
        try:
            fetch_res = self.variants_idx.fetch(
                r_ref_pos.chrm, r_ref_pos.start + edge_buffer,
                r_ref_pos.end - edge_buffer)
        except ValueError:
            raise mh.MegaError('Mapped location not valid for variants file.')

        snp_calls = sorted(fetch_res, key=itemgetter(0))
        # initialize snp_grp with first snp
        snp_data = snp_calls[0]
        prev_snp_end = snp_data[0] + len(snp_data[2])
        snp_grp = [snp_data]
        for variant in snp_calls[1:]:
            if self.max_indel_size is not None and max(
                    np.abs(len(variant.ref) - len(snp_alt_seq))
                    for snp_alt_seq in variant.alts) > self.max_indel_size:
                continue
            if snp_data[0] < prev_snp_end + snp_merge_gap:
                prev_snp_end = max(snp_data[0] + len(snp_data[2]), prev_snp_end)
                snp_grp.append(snp_data)
            else:
                yield filter_overlapping_dels(snp_grp)
                prev_snp_end = snp_data[0] + len(snp_data[2])
                snp_grp = [snp_data]

        # yeild last snp grp data
        yield filter_overlapping_dels(snp_grp)
        return

    def calibrate_llr(self, llr, snp_ref_seq, snp_alt_seq):
        return self.calib_table.calibrate_llr(
            llr, snp_ref_seq, snp_alt_seq)


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
        self.id = '' if id is None else str(id)
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
            sorted_keys.insert(0, sorted_keys.pop(sorted_keys.index('GT')))
        if 'LOG_PROBS' in sorted_keys:
            # move log probs to end of format field for easier human readability
            sorted_keys.append(sorted_keys.pop(sorted_keys.index('LOG_PROBS')))
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

    def add_haploid_probs(self, probs, gts):
        # phred scaled likelihoods
        with np.errstate(divide='ignore'):
            gl = np.log10(probs)
        raw_pl = -10 * gl
        # "normalized" PL values stored as decsribed by VCF format
        # abs to remove negative 0 from file
        pl = np.abs(np.minimum(raw_pl - raw_pl.min(), mh.MAX_PL_VALUE))
        s_pl = np.sort(pl)

        # add sample tags
        self.add_sample_field('GT', gts[np.argmax(probs)])
        qual = int(np.around(np.minimum(raw_pl[0], mh.MAX_PL_VALUE)))
        self.qual = '{:.0f}'.format(qual) if qual > 0 else '.'
        self.add_sample_field('GQ', '{:.0f}'.format(np.around(s_pl[1])))
        self.add_sample_field(
            'GL', ','.join('{:.2f}' for _ in range(probs.shape[0])).format(*gl))
        self.add_sample_field(
            'PL', ','.join('{:.0f}' for _ in range(probs.shape[0])).format(
                *np.around(pl)))
        return

    def add_diploid_probs(self, probs, gts):
        # phred scaled likelihoods
        with np.errstate(divide='ignore'):
            gl = np.log10(probs)
        raw_pl = -10 * gl
        # "normalized" PL values stored as decsribed by VCF format
        pl = np.minimum(raw_pl - raw_pl.min(), mh.MAX_PL_VALUE)
        s_pl = np.sort(pl)

        # add sample tags
        self.add_sample_field('GT', gts[np.argmax(probs)])
        try:
            qual = int(np.minimum(np.around(raw_pl[0]), mh.MAX_PL_VALUE))
        except ValueError:
            logger = logging.get_logger()
            logger.debug(
                'NAN quality value encountered. gts:{}, probs:{}'.format(
                    str(gts), str(probs)))
            qual = mg.MAX_PL_VALUE
        self.qual = '{:.0f}'.format(qual) if qual > 0 else '.'
        self.add_sample_field('GQ', '{:.0f}'.format(np.around(s_pl[1])))
        self.add_sample_field(
            'GL', ','.join('{:.2f}' for _ in range(probs.shape[0])).format(*gl))
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
    version_options = set(['4.1',])
    def __init__(
            self, filename, mode='w',
            header=('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
                    'INFO', 'FORMAT', SAMPLE_NAME),
            extra_meta_info=FIXED_VCF_MI, version='4.1', ref_fn=None,
            ref_names_and_lens=None, write_vcf_lp=False):
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
        if write_vcf_lp:
            self.meta.append(FORMAT_LOG_PROB_MI)

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
            # order by the inverse log likelihood ratio
            llr_ord = np.argsort(all_lps[a2] - all_lps[a1])
            s_a1_lps = all_lps[a1][llr_ord]
            s_a2_lps = all_lps[a2][llr_ord]
            # compute log probability of heterozygous genotype by binomial
            # weighted sum of maximum likelihoods
            return logsumexp(np.array([
                np.log(binom_pmf(i, all_lps.shape[1], 0.5)) +
                np.sum(s_a1_lps[:i]) + np.sum(s_a2_lps[i:])
                for i in range(all_lps.shape[1] + 1)]))


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
        snp_lps = np.concatenate([[ref_lps.sum()], alts_lps.sum(axis=1)])
        post_snp_lps = snp_lps - logsumexp(snp_lps)
        return np.exp(post_snp_lps), list(map(str, range(snp_lps.shape[0])))

    def compute_snp_stats(
            self, snp_loc, het_factors, call_mode=DIPLOID_MODE,
            valid_read_ids=None):
        assert call_mode in (HAPLIOD_MODE, DIPLOID_MODE), (
            'Invalid SNP aggregation ploidy call mode: {}.'.format(call_mode))

        pr_snp_stats = self.get_per_read_snp_stats(snp_loc)
        alt_seqs = sorted(set(r_stats.alt_seq for r_stats in pr_snp_stats))
        pr_alt_lps = defaultdict(dict)
        for r_stats in pr_snp_stats:
            if (valid_read_ids is not None and
                r_stats.read_id not in valid_read_ids):
                continue
            pr_alt_lps[r_stats.read_id][r_stats.alt_seq] = r_stats.score
        if len(pr_alt_lps) == 0:
            raise mh.MegaError('No valid reads cover SNP')

        alt_seq_lps = [[] for _ in range(len(alt_seqs))]
        for read_lps in pr_alt_lps.values():
            for i, alt_seq in enumerate(alt_seqs):
                try:
                    alt_seq_lps[i].append(read_lps[alt_seq])
                except KeyError:
                    raise mh.MegaError(
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
            snp_var.add_sample_field('LOG_PROBS', ','.join(
                ';'.join('{:.2f}'.format(lp) for lp in alt_i_lps)
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


#########################################
##### Whatshap Mapping Post-process #####
#########################################

def get_whatshap_command(index_variant_fn, whatshap_sort_fn, phase_fn):
    return ('Run following command to obtain phased variants:\n\t\t' +
            'whatshap phase --indels --distrust-genotypes -o {} {} {}').format(
                phase_fn, index_variant_fn, whatshap_sort_fn)

def sort_variants(in_vcf_fn, out_vcf_fn):
    in_vcf_fp = pysam.VariantFile(in_vcf_fn)
    with pysam.VariantFile(
            out_vcf_fn, 'w', header=in_vcf_fp.header) as out_vcf_fp:
        for rec in sorted(in_vcf_fp.fetch(), key=lambda r: (r.chrom, r.start)):
            out_vcf_fp.write(rec)
    return

def index_variants(variant_fn):
    try:
        return pysam.tabix_index(
            variant_fn, force=True, preset='vcf', keep_original=True)
    except OSError:
        # file likely not sorted
        sort_variant_fn = mh.add_fn_suffix(variant_fn, 'sorted')
        sort_variants(variant_fn, sort_variant_fn)
        return pysam.tabix_index(
            sort_variant_fn, force=True, preset='vcf', keep_original=True)


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
