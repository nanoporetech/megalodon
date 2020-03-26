import os
import sys
import shutil
import pkg_resources
import multiprocessing as mp
from abc import ABC, abstractmethod
from collections import namedtuple, OrderedDict
from multiprocessing.queues import Queue as mpQueue

import numpy as np


# TODO move these values into model specific files as they may change from
# model to model. Need to create script to automate HET_FACTOR optimization
# process first
# Also potentially use sepatate insertion and deletion factor
# determine necessity from larger validation run
DEFAULT_SNV_HET_FACTOR = 2.1
DEFAULT_INDEL_HET_FACTOR = 1.6

DEFAULT_EDGE_BUFFER = 2
CONTEXT_MAX_DIST = 5
DEFAULT_MAX_INDEL_SIZE = 50
DEFUALT_MAX_VAR_CNTXTS = 16
DEFAULT_SNV_CONTEXT = 5
DEFAULT_INDEL_CONTEXT = 10
DEFAULT_VAR_CONTEXT_BASES = [DEFAULT_SNV_CONTEXT, DEFAULT_INDEL_CONTEXT]
DEFAULT_MOD_CONTEXT = 5
DEFAULT_CONTEXT_MIN_ALT_PROB = 0.05

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

_MAX_QUEUE_SIZE = 10000

# allow 64GB for memory mapped sqlite file access
MEMORY_MAP_LIMIT = 64000000000
SQLITE_TIMEOUT = 5

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
BC_OUT_FMTS = ('fastq', 'fasta')
BC_MODS_NAME = 'mod_basecalls'
MAP_NAME = 'mappings'
MAP_SUMM_NAME = 'mappings_summary'
MAP_OUT_FMTS = ('bam', 'cram', 'sam')
PR_VAR_NAME = 'per_read_variants'
PR_VAR_TXT_NAME = 'per_read_variants_text'
WHATSHAP_MAP_NAME = 'whatshap_mappings'
VAR_NAME = 'variants'
PR_MOD_NAME = 'per_read_mods'
PR_MOD_TXT_NAME = 'per_read_mods_text'
MOD_NAME = 'mods'
SIG_MAP_NAME = 'signal_mappings'
PR_REF_NAME = 'per_read_refs'
OUTPUT_FNS = {
    BC_NAME: 'basecalls',
    BC_MODS_NAME: 'basecalls.modified_base_scores.hdf5',
    MAP_NAME: 'mappings',
    MAP_SUMM_NAME: 'mappings.summary.txt',
    PR_VAR_NAME: 'per_read_variant_calls.db',
    PR_VAR_TXT_NAME: 'per_read_variant_calls.txt',
    VAR_NAME: 'variants.vcf',
    WHATSHAP_MAP_NAME: 'whatshap_mappings',
    PR_MOD_NAME: 'per_read_modified_base_calls.db',
    PR_MOD_TXT_NAME: 'per_read_modified_base_calls.txt',
    MOD_NAME: 'modified_bases',
    SIG_MAP_NAME: 'signal_mappings.hdf5',
    PR_REF_NAME: 'per_read_references.fasta'
}
LOG_FILENAME = 'log.txt'
# outputs to be selected with command line --outputs argument
OUTPUT_DESCS = OrderedDict([
    (BC_NAME, 'Called bases (FASTA/Q)'),
    (BC_MODS_NAME, 'Basecall-anchored modified base scores (HDF5)'),
    (MAP_NAME, 'Mapped reads (BAM/CRAM/SAM)'),
    (PR_VAR_NAME, 'Per-read, per-site sequence variant scores database'),
    (VAR_NAME, 'Sample-level aggregated sequence variant calls (VCF)'),
    (WHATSHAP_MAP_NAME,
     'Sequence variant annotated mappings for use with whatshap'),
    (PR_MOD_NAME, 'Per-read, per-site modified base scores database'),
    (MOD_NAME, 'Sample-level aggregated modified base calls (modVCF)'),
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

ALIGN_OUTPUTS = set((MAP_NAME, PR_REF_NAME, SIG_MAP_NAME, PR_VAR_NAME,
                     VAR_NAME, WHATSHAP_MAP_NAME, PR_MOD_NAME, MOD_NAME))
GETTER_PROC = namedtuple('getter_proc', ('queue', 'proc', 'conn'))

REF_OUT_INFO = namedtuple('ref_out_info', (
    'pct_idnt', 'pct_cov', 'min_len', 'max_len', 'alphabet',
    'collapse_alphabet', 'annotate_mods', 'annotate_vars', 'mod_thresh',
    'output_sig_maps', 'output_pr_refs'))
REF_OUT_INFO.__new__.__defaults__ = (None, None, None, None, False, False)

# directory names define model preset string
# currently only one model trained
MODEL_DATA_DIR_NAME = 'model_data'
VAR_CALIBRATION_FN = 'megalodon_variant_calibration.npz'
MOD_CALIBRATION_FN = 'megalodon_mod_calibration.npz'


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


#######################
# Filename Extraction #
#######################

def resolve_path(fn_path):
    """Helper function to resolve relative and linked paths that might
    give other packages problems.
    """
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

    return


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
                'megalodon installation is valid.')
    out_msg = ('Megalodon support for guppy configs (basecalling and ' +
               'mapping supported for flip-flop configs):\n' +
               'Variant Support    Modbase Support    Config\n')
    for config in configs:
        config_files = os.listdir(resolve_path(
            pkg_resources.resource_filename('megalodon', os.path.join(
                MODEL_DATA_DIR_NAME, config))))
        out_msg += 'X' + ' ' * 18 if VAR_CALIBRATION_FN in config_files else \
                   ' ' * 19
        out_msg += 'X' + ' ' * 18 if MOD_CALIBRATION_FN in config_files else \
                   ' ' * 19
        out_msg += config + '\n'
    return out_msg


###########################
# Multi-processing Helper #
###########################

class CountingMPQueue(mpQueue):
    """ Minimal version of multiprocessing queue maintaining a queue size
    counter
    """
    def __init__(self, **kwargs):
        super().__init__(ctx=mp.get_context(), **kwargs)
        self.size = mp.Value('i', 0)
        self.maxsize = None
        if 'maxsize' in kwargs:
            self.maxsize = kwargs['maxsize']

    def put(self, *args, **kwargs):
        super().put(*args, **kwargs)
        with self.size.get_lock():
            self.size.value += 1

    def get(self, *args, **kwargs):
        rval = super().get(*args, **kwargs)
        with self.size.get_lock():
            self.size.value -= 1
        return rval

    def qsize(self):
        qsize = max(0, self.size.value)
        if self.maxsize is not None:
            return min(self.maxsize, qsize)
        return qsize

    def empty(self):
        return self.qsize() <= 0


def create_getter_q(getter_func, args, max_size=_MAX_QUEUE_SIZE):
    if max_size is None:
        q = CountingMPQueue()
    else:
        q = CountingMPQueue(maxsize=max_size)
    main_conn, conn = mp.Pipe()
    p = mp.Process(target=getter_func, daemon=True, args=(q, conn, *args))
    p.start()
    return GETTER_PROC(q, p, main_conn)


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

    :param data: A :class:`ndarray` object
    :param factor: Factor to scale MAD by. Default (None) is to be consistent
    with the standard deviation of a normal distribution
    (i.e. mad( N(0, sigma^2) ) = sigma).
    :param axis: For multidimensional arrays, which axis to calculate over
    :param keepdims: If True, axis is kept as dimension of length 1

    :returns: a tuple containing the median and MAD of the data
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


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
