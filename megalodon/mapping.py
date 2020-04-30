import sys
import queue
from time import sleep
from collections import namedtuple

import mappy
import pysam
import numpy as np

from megalodon import megalodon_helper as mh, logging


MAP_POS = namedtuple('MAP_POS', (
    'chrm', 'strand', 'start', 'end', 'q_trim_start', 'q_trim_end'))
MAP_SUMM = namedtuple('MAP_SUMM', (
    'read_id', 'pct_identity', 'num_align', 'num_match',
    'num_del', 'num_ins', 'read_pct_coverage', 'chrom', 'strand',
    'start', 'end'))
# Defaults for backwards compatibility when reading
MAP_SUMM.__new__.__defaults__ = (None, None, None, None, None)
MAP_SUMM_TMPLT = (
    '{0.read_id}\t{0.pct_identity:.2f}\t{0.num_align}\t{0.num_match}\t' +
    '{0.num_del}\t{0.num_ins}\t{0.read_pct_coverage:.2f}\t{0.chrom}\t' +
    '{0.strand}\t{0.start}\t{0.end}\n')
MAP_SUMM_TYPES = dict(zip(
    MAP_SUMM._fields,
    (str, float, int, int, int, int, float, str, str, int, int)))

LOGGER = logging.get_logger()


# mappy aligner with extra attrs
class alignerPlus(mappy.Aligner):
    ref_names_and_lens = None
    out_fmt = None
    ref_fn = None

    def add_ref_lens(self):
        ref_names, ref_lens = [], []
        for ref_name in self.seq_names:
            ref_names.append(ref_name)
            ref_lens.append(len(self.seq(ref_name)))
        self.ref_names_and_lens = (ref_names, ref_lens)

        return


def align_read(q_seq, aligner, map_thr_buf, read_id=None):
    try:
        # enumerate all alignments to avoid memory leak from mappy
        r_algn = list(aligner.map(str(q_seq), buf=map_thr_buf))[0]
    except IndexError:
        # alignment not produced
        return [None, None], None

    ref_seq = aligner.seq(r_algn.ctg, r_algn.r_st, r_algn.r_en)
    if r_algn.strand == -1:
        ref_seq = mh.revcomp(ref_seq)
    r_algn_data = [
        r_algn.ctg, r_algn.strand, r_algn.r_st, r_algn.r_en,
        r_algn.q_st, r_algn.q_en, r_algn.cigar]
    return [ref_seq, r_algn_data], (
        read_id, q_seq, r_algn.ctg, r_algn.strand, r_algn.r_st,
        r_algn.q_st, r_algn.q_en, r_algn.cigar)


def _map_read_worker(aligner, map_conn, mo_q):
    # get mappy aligner thread buffer
    map_thr_buf = mappy.ThreadBuffer()

    while True:
        try:
            q_seq, read_id = map_conn.recv()
        except Exception:
            # exit gracefully
            return
        if q_seq is None:
            break
        map_res, full_res = align_read(
            q_seq, aligner, map_thr_buf, read_id)
        map_conn.send(map_res)

        if mo_q is not None and full_res is not None:
            mo_q.put((map_res[0], full_res))

    return


def parse_cigar(r_cigar, strand, ref_len):
    fill_invalid = -1
    # get each base calls genomic position
    r_to_q_poss = np.full(ref_len + 1, fill_invalid, dtype=np.int32)
    # process cigar ops in read direction
    curr_r_pos, curr_q_pos = 0, 0
    cigar_ops = r_cigar if strand == 1 else r_cigar[::-1]
    for op_len, op in cigar_ops:
        if op == 1:
            # inserted bases into ref
            curr_q_pos += op_len
        elif op in (2, 3):
            # deleted ref bases
            for r_pos in range(curr_r_pos, curr_r_pos + op_len):
                r_to_q_poss[r_pos] = curr_q_pos
            curr_r_pos += op_len
        elif op in (0, 7, 8):
            # aligned bases
            for op_offset in range(op_len):
                r_to_q_poss[curr_r_pos + op_offset] = curr_q_pos + op_offset
            curr_q_pos += op_len
            curr_r_pos += op_len
        elif op == 6:
            # padding (shouldn't happen in mappy)
            pass
    r_to_q_poss[curr_r_pos] = curr_q_pos
    if r_to_q_poss[-1] == fill_invalid:
        raise mh.MegaError((
            'Invalid cigar string encountered. Reference length: {}  Cigar ' +
            'implied reference length: {}').format(ref_len, curr_r_pos))

    return r_to_q_poss


def map_read(q_seq, read_id, caller_conn, signal_reversed=False):
    """Map read (query) sequence and return:
    1) reference sequence (endcoded as int labels)
    2) mapping from reference to read positions (after trimming)
    3) reference mapping position (including read trimming positions)
    4) cigar as produced by mappy
    """
    # send seq to _map_read_worker and receive mapped seq and pos
    if signal_reversed:
        q_seq = q_seq[::-1]
    caller_conn.send((q_seq, read_id))
    r_ref_seq, r_algn = caller_conn.recv()
    if r_ref_seq is None:
        raise mh.MegaError('No alignment')
    chrm, strand, r_st, r_en, q_st, q_en, r_cigar = r_algn
    if signal_reversed:
        q_st, q_en = len(q_seq) - q_en, len(q_seq) - q_st
        r_ref_seq = r_ref_seq[::-1]
        r_cigar = r_cigar[::-1]

    try:
        r_to_q_poss = parse_cigar(r_cigar, strand, r_en - r_st)
    except mh.MegaError as e:
        LOGGER.debug('Read {} ::: '.format(read_id) + str(e))
        raise mh.MegaError('Invalid cigar string encountered.')
    r_pos = MAP_POS(
        chrm=chrm, strand=strand, start=r_st, end=r_en,
        q_trim_start=q_st, q_trim_end=q_en)

    return r_ref_seq, r_to_q_poss, r_pos, r_cigar


def compute_pct_identity(cigar):
    nalign, nmatch = 0, 0
    for op_len, op in cigar:
        if op not in (4, 5):
            nalign += op_len
        if op in (0, 7):
            nmatch += op_len
    return 100 * nmatch / float(nalign)


def read_passes_filters(ref_out_info, read_len, q_st, q_en, cigar):
    if ref_out_info.min_len is not None and read_len < ref_out_info.min_len:
        return False
    if ref_out_info.max_len is not None and read_len > ref_out_info.max_len:
        return False
    if ref_out_info.pct_cov is not None and \
       100 * (q_en - q_st) / read_len < ref_out_info.pct_cov:
        return False
    if ref_out_info.pct_idnt is not None and \
       compute_pct_identity(cigar) < ref_out_info.pct_idnt:
        return False
    return True


def open_alignment_out_file(out_dir, map_fmt, ref_names_and_lens, ref_fn):
    map_fn = mh.get_megalodon_fn(out_dir, mh.MAP_NAME) + '.' + map_fmt
    if map_fmt == 'bam':
        w_mode = 'wb'
    elif map_fmt == 'cram':
        w_mode = 'wc'
    elif map_fmt == 'sam':
        w_mode = 'w'
    else:
        raise mh.MegaError('Invalid mapping output format')
    return pysam.AlignmentFile(
        map_fn, w_mode, reference_names=ref_names_and_lens[0],
        reference_lengths=ref_names_and_lens[1], reference_filename=ref_fn)


def test_open_alignment_out_file(out_dir, map_fmt, ref_names_and_lens, ref_fn):
    try:
        map_fp = open_alignment_out_file(
            out_dir, map_fmt, ref_names_and_lens, ref_fn)
    except ValueError:
        raise mh.MegaError(
            'Failed to open alignment file for writing. Check that ' +
            'reference file is compressed with bgzip for CRAM output.')
    map_fp.close()
    return


def _get_map_queue(
        mo_q, map_conn, out_dir, ref_names_and_lens, map_fmt, ref_fn,
        do_output_pr_refs, ref_out_info):
    def write_alignment(
            read_id, q_seq, chrm, strand, r_st, q_st, q_en, cigar):
        bc_len = len(q_seq)
        q_seq = q_seq[q_st:q_en]

        a = pysam.AlignedSegment()
        a.query_name = read_id
        a.query_sequence = q_seq if strand == 1 else mh.revcomp(q_seq)
        a.flag = 0 if strand == 1 else 16
        a.reference_id = map_fp.get_tid(chrm)
        a.reference_start = r_st
        a.cigartuples = [(op, op_l) for op_l, op in cigar]
        a.template_length = q_en - q_st
        map_fp.write(a)

        nalign, nmatch, ndel, nins = [0, ] * 4
        for op_len, op in cigar:
            if op not in (4, 5):
                nalign += op_len
            if op in (0, 7):
                nmatch += op_len
            elif op in (2, 3):
                ndel += op_len
            elif op == 1:
                nins += op_len
        # compute alignment stats
        r_map_summ = MAP_SUMM(
            read_id=read_id, pct_identity=100 * nmatch / float(nalign),
            num_align=nalign, num_match=nmatch, num_del=ndel, num_ins=nins,
            read_pct_coverage=(q_en - q_st) * 100 / float(bc_len), chrom=chrm,
            strand=mh.int_strand_to_str(strand), start=r_st,
            end=r_st + nalign - nins)
        summ_fp.write(MAP_SUMM_TMPLT.format(r_map_summ))
        summ_fp.flush()

        return

    def write_pr_ref(read_id, ref_seq):
        pr_ref_fp.write('>{}\n{}\n'.format(read_id, ref_seq))
        pr_ref_fp.flush()
        return

    def get_alignment():
        ref_seq, (read_id, q_seq, chrm, strand, r_st, q_st, q_en,
                  cigar) = mo_q.get(block=False)
        write_alignment(read_id, q_seq, chrm, strand, r_st, q_st, q_en, cigar)
        if do_output_pr_refs and read_passes_filters(
                ref_out_info, len(q_seq), q_st, q_en, cigar):
            write_pr_ref(read_id, ref_seq)
        return

    summ_fp = open(mh.get_megalodon_fn(out_dir, mh.MAP_SUMM_NAME), 'w')
    summ_fp.write('\t'.join(MAP_SUMM._fields) + '\n')

    map_fp = open_alignment_out_file(
        out_dir, map_fmt, ref_names_and_lens, ref_fn)

    if do_output_pr_refs:
        pr_ref_fp = open(mh.get_megalodon_fn(out_dir, mh.PR_REF_NAME), 'w')

    try:
        while True:
            try:
                get_alignment()
            except queue.Empty:
                if map_conn.poll():
                    break
                sleep(0.001)
                continue

        while not mo_q.empty():
            get_alignment()
    finally:
        map_fp.close()
        summ_fp.close()
        if do_output_pr_refs:
            pr_ref_fp.close()

    return


####################
# Samtools wrapper #
####################

def sort_and_index_mapping(map_fn, out_fn, ref_fn=None, do_index=False):
    sort_args = ['-O', 'BAM', '-o', out_fn, map_fn]
    if ref_fn is not None:
        sort_args.extend(('--reference', ref_fn))
    try:
        sleep(0.1)
        pysam.sort(*sort_args)
        if do_index:
            sleep(0.1)
            pysam.index(out_fn, catch_stdout=False, catch_stderr=False)
    except pysam.utils.SamtoolsError:
        LOGGER.warning('Sorting and/or indexing mapping failed.')

    return


##########################
# Mapping summary parser #
##########################

def parse_map_summary_file(map_summ_fn):
    def parse_line(line):
        return MAP_SUMM(*(None if v is None else MAP_SUMM_TYPES[fn](v)
                          for fn, v in zip(header, line.split())))

    with open(map_summ_fn) as map_summ_fp:
        header = map_summ_fp.readline().split()
        map_summ = [parse_line(line) for line in map_summ_fp]
    return map_summ


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
