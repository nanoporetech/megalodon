import os
import pkg_resources
import multiprocessing as mp
from collections import namedtuple
from abc import ABC, abstractmethod

import numpy as np


# TODO move these values into model specific files as they may change from
# model to model. Need to create script to automate HET_FACTOR optimization
# process first
# Also potentially use sepatate insertion and deletion factor
# determine necessity from larger validation run
DEFAULT_SNV_HET_FACTOR = 0.85
DEFAULT_INDEL_HET_FACTOR = 0.78

MED_NORM_FACTOR = 1.4826

ALPHABET = 'ACGT'
COMP_BASES = dict(zip(map(ord, 'ACGT'), map(ord, 'TGCA')))

_MAX_QUEUE_SIZE = 1000

# VCF spec text
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
MAP_OUT_FMTS = ('bam', 'cram', 'sam')
PR_SNP_NAME = 'per_read_snps'
SNP_NAME = 'snps'
PR_MOD_NAME = 'per_read_mods'
# TOOD add wig/bedgraph modified base output
MOD_NAME = 'mods'
OUTPUT_FNS = {
    BC_NAME:'basecalls',
    BC_MODS_NAME:'basecalls.modified_base_scores.hdf5',
    MAP_NAME:['mappings', 'mappings.summary.txt'],
    PR_SNP_NAME:['per_read_snp_calls.db',
                 'per_read_snp_calls.txt'],
    SNP_NAME:'variants.vcf',
    PR_MOD_NAME:['per_read_modified_base_calls.db',
                 'per_read_modified_base_calls.txt'],
    MOD_NAME:'modified_bases.mvcf'
}
OUTPUT_DESCS = [
    (BC_NAME, 'Called bases (FASTA)'),
    (BC_MODS_NAME, 'Basecall-anchored modified base scores (HDF5)'),
    (MAP_NAME, 'Mapped reads (BAM/CRAM/SAM)'),
    (PR_SNP_NAME, 'Per-read, per-site SNP scores database'),
    (SNP_NAME, 'Sample-level aggregated SNP calls (VCF)'),
    (PR_MOD_NAME, 'Per-read, per-site modified base scores database'),
    (MOD_NAME, 'Sample-level aggregated modified base calls (modVCF)')
]
LOG_FILENAME = 'log.txt'
# special output type, not included in standard --outputs (since it is
# used only in special circumstances)
PR_REF_NAME = 'per_read_ref'
PR_REF_FN = 'per_read_references.fasta'

ALIGN_OUTPUTS = set((MAP_NAME, PR_REF_NAME, PR_SNP_NAME, SNP_NAME,
                     PR_MOD_NAME, MOD_NAME))

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

def nstate_to_nbase(nstate):
    return int(np.sqrt(0.25 + (0.5 * nstate)) - 0.5)

def comp(seq):
    return seq.translate(COMP_BASES)

def revcomp(seq):
    return seq.translate(COMP_BASES)[::-1]

def resolve_path(fn_path):
    """Helper function to resolve relative and linked paths that might
    give other packages problems.
    """
    return os.path.realpath(os.path.expanduser(fn_path))


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

def get_model_fn(model_fn, preset_str=None):
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

def create_getter_q(getter_func, args):
    q = mp.Queue(maxsize=_MAX_QUEUE_SIZE)
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
