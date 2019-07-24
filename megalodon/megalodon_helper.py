import os
import pkg_resources
import multiprocessing as mp
from collections import namedtuple
from abc import ABC, abstractmethod

import pysam
import numpy as np


# TODO move these values into model specific files as they may change from
# model to model. Need to create script to automate HET_FACTOR optimization
# process first
# Also potentially use sepatate insertion and deletion factor
# determine necessity from larger validation run
DEFAULT_SNV_HET_FACTOR = 2.1
DEFAULT_INDEL_HET_FACTOR = 1.6

CONTEXT_MAX_DIST = 5
DEFAULT_SNV_CONTEXT = 15
DEFAULT_INDEL_CONTEXT = 30
DEFAULT_MOD_CONTEXT = 15

MED_NORM_FACTOR = 1.4826

ALPHABET = 'ACGT'
NP_ALPHABET = np.char.array(list(ALPHABET))
COMP_BASES = dict(zip(map(ord, 'ACGT'), map(ord, 'TGCA')))

_MAX_QUEUE_SIZE = 1000

# VCF spec text
MIN_GL_VALUE = -1000
MAX_PL_VALUE = 255
VCF_VERSION_MI = 'fileformat=VCFv{}'
FILE_DATE_MI = 'fileDate={}'
SOURCE_MI = 'source=ont-megalodon.v{}'
REF_MI = "reference={}"
CONTIG_MI = "contig=<ID={},length={}>"

# outputs specification
BC_NAME = 'basecalls'
BC_OUT_FMTS = ('fasta',)
BC_MODS_NAME = 'mod_basecalls'
MAP_NAME = 'mappings'
MAP_SUMM_NAME = 'mappings_summary'
MAP_OUT_FMTS = ('bam', 'cram', 'sam')
PR_SNP_NAME = 'per_read_snps'
PR_SNP_TXT_NAME = 'per_read_snps_text'
WHATSHAP_MAP_NAME = 'whatshap_mappings'
SNP_NAME = 'snps'
PR_MOD_NAME = 'per_read_mods'
PR_MOD_TXT_NAME = 'per_read_mods_text'
# TOOD add wig/bedgraph modified base output
MOD_NAME = 'mods'
PR_REF_NAME = 'per_read_ref'
OUTPUT_FNS = {
    BC_NAME:'basecalls',
    BC_MODS_NAME:'basecalls.modified_base_scores.hdf5',
    MAP_NAME:'mappings',
    MAP_SUMM_NAME:'mappings.summary.txt',
    PR_SNP_NAME:'per_read_snp_calls.db',
    PR_SNP_TXT_NAME:'per_read_snp_calls.txt',
    SNP_NAME:'variants.vcf',
    WHATSHAP_MAP_NAME:'whatshap_mappings',
    PR_MOD_NAME:'per_read_modified_base_calls.db',
    PR_MOD_TXT_NAME:'per_read_modified_base_calls.txt',
    MOD_NAME:'modified_bases',
    PR_REF_NAME:'per_read_references.fasta'
}
LOG_FILENAME = 'log.txt'
# outputs to be selected with command line --outputs argument
OUTPUT_DESCS = {
    BC_NAME:'Called bases (FASTA)',
    BC_MODS_NAME:'Basecall-anchored modified base scores (HDF5)',
    MAP_NAME:'Mapped reads (BAM/CRAM/SAM)',
    PR_SNP_NAME:'Per-read, per-site SNP scores database',
    SNP_NAME:'Sample-level aggregated SNP calls (VCF)',
    WHATSHAP_MAP_NAME:'SNP annotated mappings for use with whatshap',
    PR_MOD_NAME:'Per-read, per-site modified base scores database',
    MOD_NAME:'Sample-level aggregated modified base calls (modVCF)'
}

# output formats for modified bases and file extensions
MOD_BEDMETHYL_NAME = 'bedmethyl'
MOD_VCF_NAME = 'modvcf'
MOD_WIG_NAME = 'wiggle'
MOD_OUTPUT_FMTS = {
    MOD_BEDMETHYL_NAME:'bedMethyl',
    MOD_VCF_NAME:'modVcf',
    MOD_WIG_NAME:'wig'
}
MOD_OUTPUT_EXTNS = {
    MOD_BEDMETHYL_NAME:'bed',
    MOD_VCF_NAME:'vcf',
    MOD_WIG_NAME:'wig'
}

ALIGN_OUTPUTS = set((MAP_NAME, PR_REF_NAME, PR_SNP_NAME, SNP_NAME,
                     WHATSHAP_MAP_NAME, PR_MOD_NAME, MOD_NAME))

PR_REF_FILTERS = namedtuple(
    'pr_ref_filters', ('pct_idnt', 'pct_cov', 'min_len', 'max_len'))

# directory names define model preset string
# currently only one model trained
MODEL_PRESETS = ['R941.min.high_acc.5mC_6mA_bio_cntxt']
MODEL_PRESET_DESC = (
    'R9.4.1, MinION/GridION, High Accuracy, 5mC(Z) and 6mA(Y) ' +
    'in biological context model')
DEFAULT_MODEL_PRESET = MODEL_PRESETS[0]
MODEL_DATA_DIR_NAME =  'model_data'
MODEL_FN = 'model.checkpoint'
SNP_CALIBRATION_FN = 'megalodon_snp_calibration.npz'
MOD_CALIBRATION_FN = 'megalodon_mod_calibration.npz'


class MegaError(Exception):
    """ Custom megalodon error for more graceful error handling
    """
    pass

############################
##### Helper Functions #####
############################

def nstate_to_nbase(nstate):
    return int(np.sqrt(0.25 + (0.5 * nstate)) - 0.5)

def comp(seq):
    return seq.translate(COMP_BASES)

def revcomp(seq):
    return seq.translate(COMP_BASES)[::-1]

def seq_to_int(seq, alphabet=mh.ALPHABET):
    np_seq = np.array([alphabet.find(b) for b in seq], dtype=np.uintp)
    if np_seq.shape[0] > 0 and np_seq.max() >= len(alphabet):
        raise mh.MegaError('Invalid character in sequence')
    return np_seq

def int_to_seq(np_seq, alphabet=mh.ALPHABET):
    if np_seq.max() >= len(alphabet):
        raise mh.MegaError('Invalid character in sequence')
    return ''.join(mh.ALPHABET[b] for b in np_seq)


###############################
##### Filename Extraction #####
###############################

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


####################################
##### Calibration File Loading #####
####################################

def get_snp_calibration_fn(
        snp_calib_fn=None, disable_snp_calib=False, preset_str=None):
    if disable_snp_calib:
        return None
    elif snp_calib_fn is not None:
        return resolve_path(snp_calib_fn)
    elif preset_str is not None:
        if preset_str not in MODEL_PRESETS:
            raise MegaError('Invalid model preset: {}'.format(preset_str))
        resolve_path(pkg_resources.resource_filename(
            'megalodon', os.path.join(
                MODEL_DATA_DIR_NAME, preset_str, SNP_CALIBRATION_FN)))
    # else return default snp calibration file
    return resolve_path(pkg_resources.resource_filename(
        'megalodon', os.path.join(
            MODEL_DATA_DIR_NAME, DEFAULT_MODEL_PRESET, SNP_CALIBRATION_FN)))

def get_mod_calibration_fn(
        mod_calib_fn=None, disable_mod_calib=False, preset_str=None):
    if disable_mod_calib:
        return None
    elif mod_calib_fn is not None:
        return resolve_path(mod_calib_fn)
    elif preset_str is not None:
        if preset_str not in MODEL_PRESETS:
            raise MegaError('Invalid model preset: {}'.format(preset_str))
        resolve_path(pkg_resources.resource_filename(
            'megalodon', os.path.join(
                MODEL_DATA_DIR_NAME, preset_str, SNP_CALIBRATION_FN)))
    # else return default snp calibration file
    return resolve_path(pkg_resources.resource_filename(
        'megalodon', os.path.join(
            MODEL_DATA_DIR_NAME, DEFAULT_MODEL_PRESET, MOD_CALIBRATION_FN)))

def get_model_fn(model_fn=None, preset_str=None):
    if model_fn is not None:
        return resolve_path(model_fn)
    elif preset_str is not None:
        if preset_str not in MODEL_PRESETS:
            raise MegaError('Invalid model preset: {}'.format(preset_str))
        resolve_path(pkg_resources.resource_filename(
            'megalodon', os.path.join(
                MODEL_DATA_DIR_NAME, preset_str, MODEL_FN)))
    # else return default snp calibration file
    return resolve_path(pkg_resources.resource_filename(
        'megalodon', os.path.join(
            MODEL_DATA_DIR_NAME, DEFAULT_MODEL_PRESET, MODEL_FN)))


###################################
##### Multi-processing Helper #####
###################################

def create_getter_q(getter_func, args, max_size=_MAX_QUEUE_SIZE):
    if max_size is None:
        q = mp.Queue()
    else:
        q = mp.Queue(maxsize=max_size)
    main_conn, conn = mp.Pipe()
    p = mp.Process(target=getter_func, daemon=True, args=(q, conn, *args))
    p.start()
    return q, p, main_conn


################################
##### Stat Aggregation ABC #####
################################

class AbstractAggregationClass(ABC):
    @abstractmethod
    def iter_uniq(self):
        return

    @abstractmethod
    def num_uniq(self):
        return


#############################
##### Signal Extraction #####
#############################

def med_mad(data, factor=None, axis=None, keepdims=False):
    """Compute the Median Absolute Deviation, i.e., the median
    of the absolute deviations from the median, and the median

    :param data: A :class:`ndarray` object
    :param factor: Factor to scale MAD by. Default (None) is to be consistent
    with the standard deviation of a normal distribution
    (i.e. mad( N(0,\sigma^2) ) = \sigma).
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
