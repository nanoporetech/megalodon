import os
import re
import sys
import shutil
import traceback
import pkg_resources
from tqdm import tqdm
from abc import ABC, abstractmethod
from collections import defaultdict, namedtuple, OrderedDict

import numpy as np

from megalodon import logging


LOGGER = logging.get_logger()

# TODO move these values into model specific files as they may change from
# model to model. Need to create script to automate HET_FACTOR optimization
# process first
# Also potentially use sepatate insertion and deletion factor
# determine necessity from larger validation run
DEFAULT_SNV_HET_FACTOR = 1.0
DEFAULT_INDEL_HET_FACTOR = 1.0

DEFAULT_EDGE_BUFFER = 30
CONTEXT_MAX_DIST = 5
DEFAULT_MAX_INDEL_SIZE = 50
DEFUALT_MAX_VAR_CNTXTS = 16
DEFAULT_SNV_CONTEXT = 15
DEFAULT_INDEL_CONTEXT = 30
DEFAULT_VAR_CONTEXT_BASES = [DEFAULT_SNV_CONTEXT, DEFAULT_INDEL_CONTEXT]
DEFAULT_MOD_CONTEXT = 15
DEFAULT_CONTEXT_MIN_ALT_PROB = 0.0
DEFAULT_MOD_MIN_PROB = 0.01
MOD_BIN_THRESH_NAME = 'binary_threshold'
MOD_EM_NAME = 'expectation_maximization'
MOD_AGG_METHOD_NAMES = set((MOD_BIN_THRESH_NAME, MOD_EM_NAME))
DEFAULT_MOD_BINARY_THRESH = 0.8
DEFAULT_READ_ENUM_TS = 8
DEFAULT_EXTRACT_SIG_PROC = 2

MED_NORM_FACTOR = 1.4826

ALPHABET = 'ACGT'
# set RNA alphabet for use in reading guppy posterior output
# requiring assumed canonical alphabet
RNA_ALPHABET = 'ACGU'
VALID_ALPHABETS = [ALPHABET, RNA_ALPHABET]
COMP_BASES = dict(zip(map(ord, 'ACGT'), map(ord, 'TGCA')))
NP_COMP_BASES = np.array([3, 2, 1, 0], dtype=np.uintp)
SEQ_MIN = np.array(['A'], dtype='S1').view(np.uint8)[0]
SEQ_TO_INT_ARR = np.full(26, -1, dtype=np.uintp)
SEQ_TO_INT_ARR[0] = 0
SEQ_TO_INT_ARR[2] = 1
SEQ_TO_INT_ARR[6] = 2
SEQ_TO_INT_ARR[19] = 3
SINGLE_LETTER_CODE = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'B': 'CGT', 'D': 'AGT', 'H': 'ACT',
    'K': 'GT', 'M': 'AC', 'N': 'ACGT', 'R': 'AG', 'S': 'CG', 'V': 'ACG',
    'W': 'AT', 'Y': 'CT'}
PHRED_BASE = 33

CHAN_INFO_CHANNEL_SLOT = 'channel_number'
CHAN_INFO_OFFSET = 'offset'
CHAN_INFO_RANGE = 'range'
CHAN_INFO_DIGI = 'digitisation'
CHAN_INFO_SAMP_RATE = 'sampling_rate'

# allow 64GB for memory mapped sqlite file access
MEMORY_MAP_LIMIT = 64000000000
DEFAULT_VAR_DATABASE_TIMEOUT = 5
DEFAULT_MOD_DATABASE_TIMEOUT = 5
# default cache size in kilobytes
SQLITE_CACHE_SIZE = 10000
SQLITE_PAGE_SIZE = 65536
SQLITE_MAX_PAGE_COUNT = 2147483646
SQLITE_THREADS = 8

# VCF spec text
MIN_GL_VALUE = -999
MAX_PL_VALUE = 999
VCF_VERSION_MI = 'fileformat=VCFv{}'
FILE_DATE_MI = 'fileDate={}'
SOURCE_MI = 'source=ont-megalodon.v{}'
REF_MI = "reference={}"
CONTIG_MI = "contig=<ID={},length={}>"
STRAND_FIELD_NAME = 'STRAND'

# outputs specification
BC_NAME = 'basecalls'
BC_FMT_FQ = 'fastq'
BC_FMT_FA = 'fasta'
BC_OUT_FMTS = (BC_FMT_FQ, BC_FMT_FA)
SEQ_SUMM_NAME = 'seq_summary'
BC_MODS_NAME = 'mod_basecalls'
MAP_NAME = 'mappings'
MAP_SUMM_NAME = 'mappings_summary'
MAP_OUT_BAM = 'bam'
MAP_OUT_CRAM = 'cram'
MAP_OUT_SAM = 'sam'
MAP_OUT_FMTS = (MAP_OUT_BAM, MAP_OUT_CRAM, MAP_OUT_SAM)
MAP_OUT_WRITE_MODES = {MAP_OUT_BAM: 'wb', MAP_OUT_CRAM: 'wc', MAP_OUT_SAM: 'w'}
PR_VAR_NAME = 'per_read_variants'
PR_VAR_TXT_NAME = 'per_read_variants_text'
VAR_MAP_NAME = 'variant_mappings'
VAR_NAME = 'variants'
PR_MOD_NAME = 'per_read_mods'
PR_MOD_TXT_NAME = 'per_read_mods_text'
MOD_MAP_NAME = 'mod_mappings'
MOD_NAME = 'mods'
SIG_MAP_NAME = 'signal_mappings'
PR_REF_NAME = 'per_read_refs'
OUTPUT_FNS = {
    BC_NAME: 'basecalls',
    SEQ_SUMM_NAME: 'sequencing_summary.txt',
    BC_MODS_NAME: 'mod_basecalls',
    MAP_NAME: 'mappings',
    MAP_SUMM_NAME: 'mappings.summary.txt',
    PR_VAR_NAME: 'per_read_variant_calls.db',
    PR_VAR_TXT_NAME: 'per_read_variant_calls.txt',
    VAR_NAME: 'variants.vcf',
    VAR_MAP_NAME: 'variant_mappings',
    PR_MOD_NAME: 'per_read_modified_base_calls.db',
    PR_MOD_TXT_NAME: 'per_read_modified_base_calls.txt',
    MOD_MAP_NAME: 'mod_mappings',
    MOD_NAME: 'modified_bases',
    SIG_MAP_NAME: 'signal_mappings.hdf5',
    PR_REF_NAME: 'per_read_references.fasta'
}
# outputs to be selected with command line --outputs argument
OUTPUT_DESCS = OrderedDict([
    (BC_NAME, 'Called bases (FASTA/Q)'),
    (BC_MODS_NAME, 'Basecall-anchored modified base scores (HDF5)'),
    (MAP_NAME, 'Mapped reads (BAM/CRAM/SAM)'),
    (PR_VAR_NAME, 'Per-read, per-site sequence variant scores database'),
    (VAR_NAME, 'Sample-level aggregated sequence variant calls (VCF)'),
    (VAR_MAP_NAME, 'Per-read mappings annotated with variant calls'),
    (PR_MOD_NAME, 'Per-read, per-site modified base scores database'),
    (MOD_NAME, 'Sample-level aggregated modified base calls (modVCF)'),
    (MOD_MAP_NAME, 'Per-read mappings annotated with modified base calls'),
    (SIG_MAP_NAME, 'Signal mappings for taiyaki model training (HDF5)'),
    (PR_REF_NAME, 'Per-read reference sequence for model training (FASTA)')
])

# output formats for modified bases and file extensions
MOD_BEDMETHYL_NAME = 'bedmethyl'
MOD_VCF_NAME = 'modvcf'
MOD_WIG_NAME = 'wiggle'
MOD_OUTPUT_FMTS = {
    MOD_BEDMETHYL_NAME: 'bedMethyl',
    MOD_VCF_NAME: 'modVcf',
    MOD_WIG_NAME: 'wig'
}
MOD_OUTPUT_EXTNS = {
    MOD_BEDMETHYL_NAME: 'bed',
    MOD_VCF_NAME: 'vcf',
    MOD_WIG_NAME: 'wig'
}

ALIGN_OUTPUTS = set((MAP_NAME, PR_REF_NAME, SIG_MAP_NAME,
                     PR_VAR_NAME, VAR_NAME, VAR_MAP_NAME,
                     PR_MOD_NAME, MOD_NAME, MOD_MAP_NAME))
MOD_OUTPUTS = set((MOD_NAME, PR_MOD_NAME, BC_MODS_NAME, MOD_MAP_NAME))

_MAX_QUEUE_SIZE = 10000
GETTER_INFO = namedtuple('GETTER_INFO', (
    'name', 'do_output', 'func', 'args',  'max_size'))
GETTER_INFO.__new__.__defaults__ = (_MAX_QUEUE_SIZE, )
STATUS_INFO = namedtuple('STATUS_INFO', (
    'suppress_prog_bar', 'suppress_queues', 'num_prog_errs'))
INPUT_INFO = namedtuple('INPUT_INFO', (
    'fast5s_dir', 'recursive', 'num_reads', 'read_ids_fn', 'num_ps',
    'do_it_live', 'num_read_enum_ts', 'num_extract_sig_proc'))
INPUT_INFO.__new__.__defaults__ = (
    False, DEFAULT_READ_ENUM_TS, DEFAULT_EXTRACT_SIG_PROC)
BASECALL_DO_OUTPUT = namedtuple('BASECALL_DO_OUTPUT', (
    'any', 'basecalls', 'mod_basecalls'))
BASECALL_INFO = namedtuple('BASECALL_INFO', (
    'do_output', 'out_dir', 'bc_fmt', 'mod_bc_fmt', 'mod_bc_min_prob',
    'mod_long_names', 'rev_sig', 'reads_per_batch'))
REF_DO_OUTPUT = namedtuple('REF_DO_OUTPUT', (
    'pr_refs', 'can_pr_refs', 'mod_pr_refs', 'var_pr_refs',
    'sig_maps', 'can_sig_maps', 'mod_sig_maps', 'var_sig_maps'))
REF_OUT_FILTER_PARAMS = namedtuple('REF_OUT_FILTER_PARAMS', (
    'pct_idnt', 'pct_cov', 'min_len', 'max_len'))
REF_OUT_INFO = namedtuple('ref_out_info', (
    'do_output', 'filt_params', 'ref_mods_all_motifs', 'alphabet_info',
    'out_dir', 'get_sig_map_func', 'per_site_threshs'))
VAR_DO_OUTPUT = namedtuple('VAR_DO_OUTPUT', (
    'db', 'text', 'var_map'))
VAR_DO_OUTPUT.__new__.__defaults__ = (False, False, False, False)
MOD_DO_OUTPUT = namedtuple('MOD_DO_OUTPUT', (
    'db', 'text', 'mod_map'))
MOD_DO_OUTPUT.__new__.__defaults__ = (False, False, False, False)

# directory names define model preset string
# currently only one model trained
MODEL_DATA_DIR_NAME = 'model_data'
VAR_CALIBRATION_FN = 'megalodon_variant_calibration.npz'
MOD_CALIBRATION_FN = 'megalodon_mod_calibration.npz'
DEFAULT_CALIB_SMOOTH_BW = 0.8
DEFAULT_CALIB_SMOOTH_MAX = 200
DEFAULT_CALIB_SMOOTH_NVALS = 5001
DEFAULT_CALIB_MIN_DENSITY = 5e-8
DEFAULT_CALIB_DIFF_EPS = 1e-6
DEFAULT_CALIB_LLR_CLIP_BUFFER = 1

SEQ_SUMM_INFO = namedtuple('seq_summ_info', (
    'filename', 'read_id', 'run_id', 'batch_id', 'channel', 'mux',
    'start_time', 'duration', 'num_events', 'passes_filtering',
    'template_start', 'num_events_template', 'template_duration',
    'sequence_length_template', 'mean_qscore_template',
    'strand_score_template', 'median_template', 'mad_template',
    'scaling_median_template', 'scaling_mad_template'))
# set default value of None for ref, alts, ref_start and strand;
# false for has_context_base
SEQ_SUMM_INFO.__new__.__defaults__ = tuple(['NA', ] * 12)

# default guppy settings
DEFAULT_GUPPY_SERVER_PATH = './ont-guppy/bin/guppy_basecall_server'
DEFAULT_GUPPY_CFG = 'dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg'
DEFAULT_GUPPY_TIMEOUT = 30.0
DEFAULT_GUPPY_BATCH_SIZE = 50

# completed read information
READ_STATUS = namedtuple('READ_STATUS', (
    'is_err', 'do_update_prog', 'err_type', 'fast5_fn', 'err_tb', 'n_sig'))
READ_STATUS.__new__.__defaults__ = (False, True, None, None, None, 0)

TRUE_TEXT_VALUES = set(('y', 'yes', 't', 'true', 'on', '1'))
FALSE_TEXT_VALUES = set(('n', 'no', 'f', 'false', 'off', '0'))

# original modified base models used 5mC=Z and 6mA=Y. This has now been
# standardized here https://github.com/samtools/hts-specs/pull/418
LEGACY_MOD_BASES = dict(zip(map(ord, 'ZY'), map(ord, 'ma')))


class MegaError(Exception):
    """ Custom megalodon error for more graceful error handling
    """
    pass


####################
# Helper Functions #
####################

def nstate_to_nbase(nstate):
    return int(np.sqrt(0.25 + (0.5 * nstate)) - 0.5)


def comp(seq):
    return seq.translate(COMP_BASES)


def revcomp(seq):
    return seq.translate(COMP_BASES)[::-1]


def comp_np(np_seq):
    return NP_COMP_BASES[np_seq]


def revcomp_np(np_seq):
    return NP_COMP_BASES[np_seq][::-1]


def seq_to_int(seq, error_on_invalid=True):
    try:
        np_seq = SEQ_TO_INT_ARR[
            np.array(list(seq), dtype='c').view(np.uint8) - SEQ_MIN]
    except IndexError:
        if error_on_invalid:
            raise MegaError('Invalid character in sequence')
        else:
            # use slower string find method to convert seq with
            # invalid characters
            np_seq = np.array([ALPHABET.find(b) for b in seq], dtype=np.uintp)
    # if error_on_invalid and np_seq.shape[0] > 0 and np_seq.max() >= 4:
    #    raise MegaError('Invalid character in sequence')
    return np_seq


def int_to_seq(np_seq, alphabet=ALPHABET):
    if np_seq.shape[0] == 0:
        return ''
    if np_seq.max() >= len(alphabet):
        raise MegaError('Invalid character in sequence')
    return ''.join(alphabet[b] for b in np_seq)


def get_mean_q_score(read_q):
    """ Extract mean q-score from FASTQ quality string
    """
    return np.mean([q_val - PHRED_BASE
                    for q_val in read_q.encode('ASCII')])


def rolling_window(a, size):
    shape = a.shape[:-1] + (a.shape[-1] - size + 1, size)
    strides = a.strides + (a. strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


def log_prob_to_phred(log_prob, ignore_np_divide=True):
    if ignore_np_divide:
        with np.errstate(divide='ignore'):
            return -10 * np.log10(1 - np.exp(log_prob))
    return -10 * np.log10(1 - np.exp(log_prob))


def log_errors(func, *args, **kwargs):
    try:
        return func(*args, **kwargs)
    except KeyboardInterrupt:
        raise
    except MegaError as e:
        LOGGER.debug('MegaError {}'.format(str(e)))
    except Exception as e:
        LOGGER.debug('UnexpectedError {}\nFull traceback:\n{}'.format(
            str(e), traceback.format_exc()))


def compile_motif_pat(raw_motif):
    return re.compile(''.join(
        '[{}]'.format(SINGLE_LETTER_CODE[letter])
        for letter in raw_motif))


def compile_rev_comp_motif_pat(raw_motif):
    return re.compile(''.join(
        '[{}]'.format(''.join(comp(b) for b in SINGLE_LETTER_CODE[letter]))
        for letter in raw_motif[::-1]))


def convert_legacy_mods(mod_base):
    return mod_base .translate(LEGACY_MOD_BASES)


#######################
# Filename Extraction #
#######################

def resolve_path(fn_path):
    """Helper function to resolve relative and linked paths that might
    give other packages problems.
    """
    if fn_path is None:
        return None
    return os.path.realpath(os.path.expanduser(fn_path))


def get_megalodon_fn(out_dir, out_type):
    return os.path.join(out_dir, OUTPUT_FNS[out_type])


def add_fn_suffix(fn, suffix):
    if suffix is not None:
        base_fn, fn_ext = os.path.splitext(fn)
        fn = base_fn + '.' + suffix + fn_ext
    return fn


def mkdir(out_dir, overwrite):
    if os.path.exists(out_dir):
        if not overwrite:
            raise MegaError(
                '--output-directory exists and --overwrite is not set.')
        if os.path.isfile(out_dir) or os.path.islink(out_dir):
            os.remove(out_dir)
        else:
            shutil.rmtree(out_dir)
    os.mkdir(out_dir)


def prep_out_fn(out_fn, overwrite):
    if os.path.exists(out_fn):
        if overwrite:
            os.remove(out_fn)
        else:
            raise NotImplementedError(
                'ERROR: Output filename exists and --overwrite not set.')
    try:
        open(out_fn, 'w').close()
        os.remove(out_fn)
    except Exception as e:
        sys.stderr.write(
            '*' * 60 + '\nERROR: Attempt to write to output filename ' +
            'location failed with the following error.\n' + '*' * 60 + '\n\n')
        raise e


############################
# Calibration File Loading #
############################

def get_var_calibration_fn(
        guppy_config=None, var_calib_fn=None, disable_var_calib=False):
    if disable_var_calib:
        return None
    if var_calib_fn is not None:
        var_calib_fn = resolve_path(var_calib_fn)
        if not os.path.exists(var_calib_fn):
            raise MegaError(
                'Sequence variants calibration file not found: {}'.format(
                    var_calib_fn))
        return var_calib_fn
    if guppy_config is not None:
        guppy_calib_fn = resolve_path(pkg_resources.resource_filename(
            'megalodon', os.path.join(
                MODEL_DATA_DIR_NAME, guppy_config, VAR_CALIBRATION_FN)))
        if not os.path.exists(guppy_calib_fn):
            raise MegaError(
                'No default sequence variant calibration file found for ' +
                'guppy config: {}'.format(guppy_config))
        return guppy_calib_fn
    raise MegaError('No valid sequence variant calibration specified.')


def get_mod_calibration_fn(
        guppy_config=None, mod_calib_fn=None, disable_mod_calib=False):
    if disable_mod_calib:
        return None
    if mod_calib_fn is not None:
        mod_calib_fn = resolve_path(mod_calib_fn)
        if not os.path.exists(mod_calib_fn):
            raise MegaError(
                'Modified base calibration file not found: {}'.format(
                    mod_calib_fn))
        return mod_calib_fn
    if guppy_config is not None:
        guppy_calib_fn = resolve_path(pkg_resources.resource_filename(
            'megalodon', os.path.join(
                MODEL_DATA_DIR_NAME, guppy_config, MOD_CALIBRATION_FN)))
        if not os.path.exists(guppy_calib_fn):
            raise MegaError(
                'No default modified base calibration file found for guppy ' +
                'config: {}'.format(guppy_config))
        return guppy_calib_fn
    raise MegaError('No valid modified base calibration specified.')


def get_supported_configs_message():
    configs = os.listdir(resolve_path(pkg_resources.resource_filename(
        'megalodon', MODEL_DATA_DIR_NAME)))
    if len(configs) == 0:
        return ('No guppy config calibration files found. Check that ' +
                'megalodon installation is valid.\n')
    out_msg = ('Megalodon support for guppy configs (basecalling and ' +
               'mapping supported for flip-flop configs):\n' +
               'Variant Support    Modbase Support    Config\n')
    for config in sorted(configs):
        config_files = os.listdir(resolve_path(
            pkg_resources.resource_filename('megalodon', os.path.join(
                MODEL_DATA_DIR_NAME, config))))
        out_msg += 'X' + ' ' * 18 if VAR_CALIBRATION_FN in config_files else \
                   ' ' * 19
        out_msg += 'X' + ' ' * 18 if MOD_CALIBRATION_FN in config_files else \
                   ' ' * 19
        out_msg += config + '\n'
    return out_msg


########################
# Stat Aggregation ABC #
########################

class AbstractAggregationClass(ABC):
    @abstractmethod
    def iter_uniq(self):
        return

    @abstractmethod
    def num_uniq(self):
        return


#####################
# Signal Extraction #
#####################

def med_mad(data, factor=None, axis=None, keepdims=False):
    """Compute the Median Absolute Deviation, i.e., the median
    of the absolute deviations from the median, and the median

    Args:
        data (np.ndarray): Data to be scaled
        factor (float): Factor to scale MAD by. Default (None) is to be
            consistent with the standard deviation of a normal distribution
            (i.e. mad( N(0, sigma^2) ) = sigma).
        axis: For multidimensional arrays, which axis to calculate over
        keepdims: If True, axis is kept as dimension of length 1

    Returns:
        A tuple containing the median and MAD of the data
    """
    if factor is None:
        factor = MED_NORM_FACTOR
    dmed = np.median(data, axis=axis, keepdims=True)
    dmad = factor * np.median(abs(data - dmed), axis=axis, keepdims=True)
    if axis is None:
        dmed = dmed.flatten()[0]
        dmad = dmad.flatten()[0]
    elif not keepdims:
        dmed = dmed.squeeze(axis)
        dmad = dmad.squeeze(axis)
    return dmed, dmad


#####################
# File-type Parsers #
#####################

def parse_read_ids(read_ids_fn):
    """ Robust read IDs file parser.

    This function will parse either a file with one read ID per line or a TSV
    file with the "read_id" field in the header.

    Args:
        read_ids_fn (str): Filename/path for read IDs parser.

    Returns:
        Set of read ID strings if successfully parsed or `None` otherwise.
    """
    if read_ids_fn is None:
        return None

    with open(read_ids_fn) as read_ids_fp:
        header = read_ids_fp.readline().strip()
        if len(header) == 0:
            LOGGER.warning(
                ('Read IDs file, {}, first line is empty. Read IDs file ' +
                 'will be ignored.').format(read_ids_fn))
            return None
        header_fields = header.split('\t')
        if len(header_fields) > 1:
            try:
                read_id_field_num = next(
                    i for i, f in enumerate(header_fields) if f == 'read_id')
            except StopIteration:
                LOGGER.warning(
                    ('Read IDs file, {}, contains multiple header fields, ' +
                     'but none are "read_id". Read IDs file will be ' +
                     'ignored.').format(read_ids_fn))
                return None

            # parse read id TSV file
            been_warned = False
            valid_read_ids = set()
            for line_num, line in enumerate(read_ids_fp):
                try:
                    valid_read_ids.add(
                        line.strip().split('\t')[read_id_field_num])
                except Exception:
                    if not been_warned:
                        LOGGER.warning(
                            ('Read IDs TSV file, {}, contains lines with ' +
                             'fewer fields than header. First offending ' +
                             'line: {}').format(read_ids_fn, line_num))
                    been_warned = True
        else:
            valid_read_ids = set(line.strip() for line in read_ids_fp)
            # if first line is not single field header assume this is a
            # read id as well.
            if header_fields[0] != 'read_id':
                valid_read_ids.add(header_fields[0])
    return valid_read_ids


def str_strand_to_int(strand_str):
    """ Convert string stand representation to integer +/-1 as used in
    minimap2/mappy
    """
    if strand_str == '+':
        return 1
    elif strand_str == '-':
        return -1
    return None


def int_strand_to_str(strand_str):
    """ Convert string stand representation to integer +/-1 as used in
    minimap2/mappy
    """
    if strand_str == 1:
        return '+'
    elif strand_str == -1:
        return '-'
    return '.'


def parse_bed_scores(bed_fn):
    bed_scores = defaultdict(dict)
    with open(bed_fn) as bed_fp:
        for line in bed_fp:
            chrm, pos, _, _, score, strand = line.split()[:6]
            bed_scores[(chrm, str_strand_to_int(strand))][
                int(pos)] = float(score)
    return dict(bed_scores)


def parse_bed_scores_np(bed_fn, ref_names_and_lens):
    bed_scores = dict(
        ((chrm, strand), np.zeros(chrm_len, dtype=np.float32))
        for chrm, chrm_len in zip(*ref_names_and_lens)
        for strand in (-1, 1))
    with open(bed_fn) as bed_fp:
        for line in bed_fp:
            chrm, pos, _, _, score, strand = line.split()[:6]
            bed_scores[(chrm, str_strand_to_int(strand))][
                int(pos)] = float(score)
    return bed_scores


def parse_beds(bed_fns, ignore_strand=False, show_prog_bar=True):
    """ Parse bed files.

    Args:
        bed_fns (Iterable): Iterable containing bed paths
        ignore_strand (bool): Set strand values to None
        show_prog_bar (bool): Show twdm progress bar

    Returns:
        Dictionary with keys (chromosome, strand) and values with set of
        0-based coordiantes.
    """
    sites = defaultdict(set)
    for bed_fn in bed_fns:
        with open(bed_fn) as bed_fp:
            bed_iter = (tqdm(bed_fp, desc=bed_fn, smoothing=0)
                        if show_prog_bar else bed_fp)
            for line in bed_iter:
                chrm, start, _, _, _, strand = line.split()[:6]
                start = int(start)
                store_strand = None if ignore_strand else \
                    str_strand_to_int(strand)
                sites[(chrm, store_strand)].add(start)

    # convert to standard dict
    sites = dict(sites)

    return sites


def parse_bed_methyls(
        bed_fns, strand_offset=None, show_prog_bar=True, valid_pos=None,
        limit=None):
    """ Parse bedmethyl files and return two dictionaries containing
    total and methylated coverage. Both dictionaries have top level keys
    (chromosome, strand) and second level keys with 0-based position.

    Args:
        bed_fns (Iterable): Bed methyl file paths
        strand_offset (bool): Set to aggregate negative strand along with
            positive strand values. Positive indicates negative strand sites
            have higher coordinate values.
        show_prog_bar (bool): Show twdm progress bar
        valid_pos (dict): Filter to valid positions, as returned from
            mh.parse_beds
        limit (int): limit the total number of sites to parse
    """
    cov = defaultdict(lambda: defaultdict(int))
    meth_cov = defaultdict(lambda: defaultdict(int))
    n_sites = 0
    for bed_fn in bed_fns:
        with open(bed_fn) as bed_fp:
            bed_iter = (tqdm(bed_fp, desc=bed_fn, smoothing=0)
                        if show_prog_bar else bed_fp)
            for line in bed_iter:
                (chrm, start, _, _, _, strand, _, _, _, num_reads,
                 pct_meth) = line.split()
                start = int(start)
                # convert to 1/-1 strand storage (matching mappy)
                store_strand = str_strand_to_int(strand)
                if strand_offset is not None:
                    # store both strand counts under None
                    store_strand = None
                    # apply offset to reverse strand positions
                    if strand == '-':
                        start -= strand_offset
                # skip any positions not found in valid_pos
                if valid_pos is not None and (
                        (chrm, store_strand) not in valid_pos or
                        start not in valid_pos[(chrm, store_strand)]):
                    continue
                num_reads = int(num_reads)
                if num_reads <= 0:
                    continue
                meth_reads = int(float(pct_meth) * num_reads / 100.0)
                cov[(chrm, store_strand)][start] += num_reads
                meth_cov[(chrm, store_strand)][start] += meth_reads
                n_sites += 1
                if limit is not None and n_sites >= limit:
                    break
        if limit is not None and n_sites >= limit:
            break

    # convert to standard dicts
    cov = dict((k, dict(v)) for k, v in cov.items())
    meth_cov = dict((k, dict(v)) for k, v in meth_cov.items())

    return cov, meth_cov


def text_to_bool(val):
    """ Convert text value to boolean.
    """
    lower_val = str(val).lower()
    if lower_val in TRUE_TEXT_VALUES:
        return True
    elif lower_val in FALSE_TEXT_VALUES:
        return False
    raise MegaError('Invalid boolean string encountered: "{}".'.format(val))


def parse_ground_truth_file(gt_data_fn, include_strand=True):
    """ Parse a ground truth data file. CSV with chrm, pos, is_mod values.
    As generated by create_mod_ground_truth.py.

    Args:
        gt_data_fn (str): Filename to be read and parsed.
        include_strand (boolean): Include strand values in parsed position
            values.

    Returns:
        Two sets of position values. First set is ground truth `True` sites and
        second are `False` sites. Values are (chrm, strand, pos) if
        include_strand is True, else values are (chrm, pos). Strand is encoded
        as +/-1 to match minimap2/mappy strands.
    """
    gt_mod_pos = set()
    gt_can_pos = set()
    with open(gt_data_fn) as fp:
        for line in fp:
            chrm, strand, pos, is_mod = line.strip().split(',')
            pos_key = (chrm, str_strand_to_int(strand), int(pos)) \
                if include_strand else (chrm, int(pos))
            if text_to_bool(is_mod):
                gt_mod_pos.add(pos_key)
            else:
                gt_can_pos.add(pos_key)
    return gt_mod_pos, gt_can_pos


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
