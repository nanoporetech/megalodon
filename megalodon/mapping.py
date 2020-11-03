import os
import sys
import array
import queue
import traceback
import subprocess
from collections import namedtuple, OrderedDict

import mappy
import pysam
import numpy as np

from megalodon import megalodon_helper as mh, logging
from megalodon._version import MEGALODON_VERSION


LOGGER = logging.get_logger()

MAP_POS = namedtuple('MAP_POS', (
    'chrm', 'strand', 'start', 'end', 'q_trim_start', 'q_trim_end'))
MAP_RES = namedtuple('MAP_RES', (
    'read_id', 'q_seq', 'ref_seq', 'ctg', 'strand', 'r_st', 'r_en',
    'q_st', 'q_en', 'cigar', 'map_sig_start', 'map_sig_end', 'sig_len'))
MAP_RES.__new__.__defaults__ = (None, None, None)
MAP_SUMM = namedtuple('MAP_SUMM', (
    'read_id', 'pct_identity', 'num_align', 'num_match',
    'num_del', 'num_ins', 'read_pct_coverage', 'chrom', 'strand',
    'start', 'end', 'map_sig_start', 'map_sig_end', 'sig_len'))
# Defaults for backwards compatibility when reading
MAP_SUMM.__new__.__defaults__ = (
    None, None, None, None, None, None, None, None)
MAP_SUMM_TMPLT = (
    '{0.read_id}\t{0.pct_identity:.2f}\t{0.num_align}\t{0.num_match}\t' +
    '{0.num_del}\t{0.num_ins}\t{0.read_pct_coverage:.2f}\t{0.chrom}\t' +
    '{0.strand}\t{0.start}\t{0.end}\t{0.map_sig_start}\t{0.map_sig_end}\t' +
    '{0.sig_len}\n')
MAP_SUMM_TYPES = dict(zip(
    MAP_SUMM._fields,
    (str, float, int, int, int, int, float, str, str, int, int,
     int, int, int)))

MOD_POS_TAG = 'Mm'
MOD_PROB_TAG = 'Ml'


def get_mapping_mode(map_fmt):
    if map_fmt == 'bam':
        return 'wb'
    elif map_fmt == 'cram':
        return 'wc'
    elif map_fmt == 'sam':
        return 'w'
    raise mh.MegaError('Invalid mapping output format: {}'.format(map_fmt))


def open_unaligned_alignment_file(basename, map_fmt, mod_long_names=None):
    fn = '{}.{}'.format(basename, map_fmt)
    header_dict = OrderedDict([('PG', [OrderedDict([
        ('ID', 'megalodon'), ('PN', 'megalodon'), ('VN', MEGALODON_VERSION),
        ('CL', ' '.join(sys.argv))])])])
    if mod_long_names is not None:
        header_dict['CO'] = [
            'Modified base "{}" encoded as "{}"'.format(
                mln, mh.convert_legacy_mods(mod_base))
            for mod_base, mln in mod_long_names]
    header = pysam.AlignmentHeader.from_dict(header_dict)
    return pysam.AlignmentFile(fn, get_mapping_mode(map_fmt), header=header,
                               add_sq_text=False)


def prepare_mapping(
        read_id, seq, flag=0, ref_id=None, ref_st=None, qual=None,
        map_qual=None, mods_scores=None, cigartuples=None, tags=None):
    a = pysam.AlignedSegment()
    a.query_name = read_id
    a.query_sequence = seq
    a.template_length = len(seq)
    a.flag = flag
    if ref_id is not None:
        a.reference_id = ref_id
    if ref_st is not None:
        a.reference_start = ref_st
    if map_qual is not None:
        a.mapping_quality = map_qual
    if qual is None:
        qual = array.array('B', [255] * len(seq))
    a.query_qualities = qual
    if cigartuples is None:
        cigartuples = [(0, len(seq)), ]
    a.cigartuples = cigartuples

    if tags is None:
        tags = []
    if mods_scores is not None:
        # Add modified base tags
        #  see https://github.com/samtools/hts-specs/pull/418
        tags.append((MOD_POS_TAG, mods_scores[0], 'Z'))
        if len(mods_scores[1]) > 0:
            tags.append((MOD_PROB_TAG, mods_scores[1]))
    a.set_tags(tags)

    return a


def get_map_pos_from_res(map_res):
    return MAP_POS(
        chrm=map_res.ctg, strand=map_res.strand, start=map_res.r_st,
        end=map_res.r_en, q_trim_start=map_res.q_st, q_trim_end=map_res.q_en)


class MapInfo:
    def __init__(
            self, aligner, map_fmt, ref_fn, out_dir, do_output_mappings,
            samtools_exec, do_sort_mappings, cram_ref_fn):
        if aligner is None:
            self.ref_names_and_lens = None
        else:
            ref_names, ref_lens = [], []
            for ref_name in aligner.seq_names:
                ref_names.append(ref_name)
                ref_lens.append(len(aligner.seq(ref_name)))
            self.ref_names_and_lens = (ref_names, ref_lens)

        self.map_fmt = map_fmt
        self.ref_fn = mh.resolve_path(ref_fn)
        self.cram_ref_fn = self.ref_fn if cram_ref_fn is None else \
            mh.resolve_path(cram_ref_fn)
        self.out_dir = out_dir
        self.do_output_mappings = do_output_mappings
        self.samtools_exec = samtools_exec
        self.do_sort_mappings = do_sort_mappings

    def open_alignment_out_file(self):
        map_fn = '{}.{}'.format(
            mh.get_megalodon_fn(self.out_dir, mh.MAP_NAME), self.map_fmt)
        w_mode = get_mapping_mode(self.map_fmt)
        try:
            align_file = pysam.AlignmentFile(
                map_fn, w_mode, reference_names=self.ref_names_and_lens[0],
                reference_lengths=self.ref_names_and_lens[1],
                reference_filename=self.cram_ref_fn)
        except ValueError:
            LOGGER.error(
                'Failed to open alignment file for writing.\n\t\tFor CRAM ' +
                'output, if FASTA is compressed ensure it is with bgzip or ' +
                'if --reference is a minimap2 index see --cram-reference.')
            raise mh.MegaError('Reference loading error.')
        return align_file

    def test_open_alignment_out_file(self):
        map_fp = self.open_alignment_out_file()
        map_fp.close()
        os.remove(map_fp.filename)

    def test_samtools(self):
        try:
            test_sort_res = subprocess.run(
                [self.samtools_exec, 'sort'],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            test_index_res = subprocess.run(
                [self.samtools_exec, 'index'],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except FileNotFoundError:
            LOGGER.warning('Samtools executable not found. Mappings will ' +
                           'not be sorted or indexed.')
            self.do_sort_mappings = False
        # Note index returns non-zero exit status
        if test_sort_res.returncode != 0:
            LOGGER.warning('Samtools test commands return non-zero exit ' +
                           'status. Mappings will not be sorted or indexed.')
            LOGGER.debug(
                ('MappingTestFail:   sort_returncode: {}   ' +
                 'index_returncode: {}\nsort_call_stdout:\n{}\n' +
                 'sort_call_stderr:\n{}\nindex_call_stdout:\n{}' +
                 '\nindex_call_stderr:\n{}').format(
                     test_sort_res.returncode, test_index_res.returncode,
                     test_sort_res.stdout.decode(),
                     test_sort_res.stderr.decode(),
                     test_index_res.stdout.decode(),
                     test_index_res.stderr.decode()))
            self.do_sort_mappings = False


def align_read(q_seq, aligner, map_thr_buf, read_id=None):
    try:
        # enumerate all alignments to avoid memory leak from mappy
        r_algn = list(aligner.map(str(q_seq), buf=map_thr_buf))[0]
    except IndexError:
        # alignment not produced
        return None

    ref_seq = aligner.seq(r_algn.ctg, r_algn.r_st, r_algn.r_en)
    if r_algn.strand == -1:
        ref_seq = mh.revcomp(ref_seq)
    return MAP_RES(
        read_id=read_id, q_seq=q_seq, ref_seq=ref_seq, ctg=r_algn.ctg,
        strand=r_algn.strand, r_st=r_algn.r_st, r_en=r_algn.r_en,
        q_st=r_algn.q_st, q_en=r_algn.q_en, cigar=r_algn.cigar)


def _map_read_worker(aligner, map_conn):
    LOGGER.debug('MappingWorkerStarting')
    # get mappy aligner thread buffer
    map_thr_buf = mappy.ThreadBuffer()

    LOGGER.debug('MappingWorkerInitComplete')
    while True:
        try:
            q_seq, read_id = map_conn.recv()
        except EOFError:
            LOGGER.debug('MappingWorkerClosing')
            break
        map_res = align_read(q_seq, aligner, map_thr_buf, read_id)
        if map_res is not None:
            # only convert to tuple if result is valid
            map_res = tuple(map_res)
        map_conn.send(map_res)


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


def map_read(
        caller_conn, called_read, sig_info, mo_q=None, signal_reversed=False,
        rl_cumsum=None):
    """ Map read (query) sequence

    Returns:
        Tuple containing
            1) reference sequence (endcoded as int labels)
            2) mapping from reference to read positions (after trimming)
            3) reference mapping position (including read trimming positions)
            4) cigar as produced by mappy
    """
    # send seq to _map_read_worker and receive mapped seq and pos
    q_seq = called_read.seq[::-1] if signal_reversed else called_read.seq
    caller_conn.send((q_seq, sig_info.read_id))
    map_res = caller_conn.recv()
    if map_res is None:
        raise mh.MegaError('No alignment')
    map_res = MAP_RES(*map_res)
    # add signal coordinates to mapping output if run-length cumsum provided
    if rl_cumsum is not None:
        # convert query start and end to signal-anchored locations
        # Note that for signal_reversed reads, the start will be larger than
        # the end
        q_st = len(map_res.q_seq) - map_res.q_st if signal_reversed else \
            map_res.q_st
        q_en = len(map_res.q_seq) - map_res.q_en if signal_reversed else \
            map_res.q_en
        map_res = map_res._replace(
            map_sig_start=called_read.trimmed_samples +
            rl_cumsum[q_st] * sig_info.stride,
            map_sig_end=called_read.trimmed_samples +
            rl_cumsum[q_en] * sig_info.stride,
            sig_len=called_read.trimmed_samples +
            rl_cumsum[-1] * sig_info.stride)
    if mo_q is not None:
        mo_q.put(tuple(map_res))
    if signal_reversed:
        # if signal is reversed compared to mapping, reverse coordinates so
        # they are relative to signal/state_data
        map_res = map_res._replace(
            q_st=len(map_res.q_seq) - map_res.q_en,
            q_en=len(map_res.q_seq) - map_res.q_st,
            ref_seq=map_res.ref_seq[::-1],
            cigar=map_res.cigar[::-1])

    try:
        r_to_q_poss = parse_cigar(
            map_res.cigar, map_res.strand, map_res.r_en - map_res.r_st)
    except mh.MegaError as e:
        LOGGER.debug('{} CigarParsingError'.format(sig_info.read_id) + str(e))
        raise mh.MegaError('Invalid cigar string encountered.')
    map_pos = get_map_pos_from_res(map_res)

    return map_res.ref_seq, r_to_q_poss, map_pos, map_res.cigar


def compute_pct_identity(cigar):
    nalign, nmatch = 0, 0
    for op_len, op in cigar:
        if op not in (4, 5):
            nalign += op_len
        if op in (0, 7):
            nmatch += op_len
    return 100 * nmatch / float(nalign)


def read_passes_filters(filt_params, read_len, q_st, q_en, cigar):
    if filt_params.min_len is not None and read_len < filt_params.min_len:
        return False
    if filt_params.max_len is not None and read_len > filt_params.max_len:
        return False
    if filt_params.pct_cov is not None and \
       100 * (q_en - q_st) / read_len < filt_params.pct_cov:
        return False
    if filt_params.pct_idnt is not None and \
       compute_pct_identity(cigar) < filt_params.pct_idnt:
        return False
    return True


def _get_map_queue(mo_q, mo_conn, map_info, ref_out_info, aux_failed_q):
    def write_alignment(map_res):
        # convert tuple back to namedtuple
        map_res = MAP_RES(*map_res)
        nalign, nmatch, ndel, nins = [0, ] * 4
        for op_len, op in map_res.cigar:
            if op not in (4, 5):
                nalign += op_len
            if op in (0, 7):
                nmatch += op_len
            elif op in (2, 3):
                ndel += op_len
            elif op == 1:
                nins += op_len
        bc_len = len(map_res.q_seq)
        q_seq = map_res.q_seq[map_res.q_st:map_res.q_en]

        a = prepare_mapping(
            map_res.read_id,
            q_seq if map_res.strand == 1 else mh.revcomp(q_seq),
            flag=0 if map_res.strand == 1 else 16,
            ref_id=map_fp.get_tid(map_res.ctg), ref_st=map_res.r_st,
            cigartuples=[(op, op_l) for op_l, op in map_res.cigar],
            tags=[('NM', nalign - nmatch)])
        map_fp.write(a)

        # compute alignment stats
        r_map_summ = MAP_SUMM(
            read_id=map_res.read_id, pct_identity=100 * nmatch / float(nalign),
            num_align=nalign, num_match=nmatch, num_del=ndel, num_ins=nins,
            read_pct_coverage=((map_res.q_en - map_res.q_st) * 100 /
                               float(bc_len)), chrom=map_res.ctg,
            strand=mh.int_strand_to_str(map_res.strand), start=map_res.r_st,
            end=map_res.r_st + nalign - nins,
            map_sig_start=map_res.map_sig_start,
            map_sig_end=map_res.map_sig_end, sig_len=map_res.sig_len)
        summ_fp.write(MAP_SUMM_TMPLT.format(r_map_summ))

        if ref_out_info.do_output.pr_refs and read_passes_filters(
                ref_out_info.filt_params, len(map_res.q_seq), map_res.q_st,
                map_res.q_en, map_res.cigar):
            pr_ref_fp.write('>{}\n{}\n'.format(
                map_res.read_id, map_res.ref_seq))

    try:
        LOGGER.debug('GetterStarting')
        # initialize file pointers
        summ_fp = open(mh.get_megalodon_fn(
            map_info.out_dir, mh.MAP_SUMM_NAME), 'w')
        summ_fp.write('\t'.join(MAP_SUMM._fields) + '\n')
        map_fp = map_info.open_alignment_out_file()
        if ref_out_info.do_output.pr_refs:
            pr_ref_fp = open(mh.get_megalodon_fn(
                map_info.out_dir, mh.PR_REF_NAME), 'w')
        workers_active = True
        LOGGER.debug('GetterInitComplete')
    except Exception as e:
        aux_failed_q.put(('MappingsInitError', str(e), traceback.format_exc()))
        return

    # loop to get alignments and write to requested files
    try:
        while workers_active or not mo_q.empty():
            try:
                map_res = mo_q.get(timeout=0.1)
                mh.log_errors(write_alignment, map_res)
            except queue.Empty:
                if mo_conn.poll():
                    workers_active = False
        LOGGER.debug('GetterClosing')
    except Exception as e:
        aux_failed_q.put((
            'MappingsProcessingError', str(e), traceback.format_exc()))
    finally:
        map_fp.close()
        summ_fp.close()
        if ref_out_info.do_output.pr_refs:
            pr_ref_fp.close()


####################
# Samtools wrapper #
####################

def sort_and_index_mapping(
        samtools_exec, map_fn, out_fn, map_fmt, ref_fn=None):
    sort_args = [
        samtools_exec, 'sort', '-O', map_fmt.upper(), '-o', out_fn, map_fn]
    if map_fmt == mh.MAP_OUT_CRAM:
        sort_args.extend(('--reference', ref_fn))
    try:
        sort_res = subprocess.run(
            sort_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if sort_res.returncode != 0:
            LOGGER.warning('Mapping sort failed. Full error text in log.')
            LOGGER.debug(
                'MappingSortFail:\ncall_stdout:\n{}\ncall_stderr:\n{}'.format(
                    sort_res.stdout.decode(), sort_res.stderr.decode()))
        if map_fmt in (mh.MAP_OUT_BAM, mh.MAP_OUT_CRAM):
            index_args = [samtools_exec, 'index', out_fn]
            index_res = subprocess.run(
                index_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if sort_res.returncode != 0:
                LOGGER.warning('Mapping index failed. Full error text in log.')
                LOGGER.debug(
                    ('MappingIndexFail:\ncall_stdout:\n{}\ncall_stderr:' +
                     '\n{}').format(
                         index_res.stdout.decode(), index_res.stderr.decode()))
    except Exception as e:
        LOGGER.warning('Sorting and/or indexing mapping failed. Full error ' +
                       'text in log.')
        LOGGER.debug('MappingSortFail:\n{}'.format(str(e)))


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
