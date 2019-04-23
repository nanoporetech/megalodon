import multiprocessing as mp
from abc import ABC, abstractmethod

import numpy as np


DEFAULT_SNV_HET_FACTOR = 0.58
DEFAULT_INDEL_HET_FACTOR = 1.3

MED_NORM_FACTOR = 1.4826

# VCF spec text
MAX_PL_VALUE = 255
VCF_VERSION_MI = 'fileformat=VCFv{}'
FILE_DATE_MI = 'fileDate={}'
SOURCE_MI = 'source=ont-megalodon.v{}'
REF_MI = "reference={}"

ALPHABET = 'ACGT'
BC_NAME = 'basecalls'
BC_OUT_FMTS = ('fasta',)
BC_MODS_NAME = 'mod_basecalls'
MAP_NAME = 'mappings'
MAP_OUT_FMTS = ('bam', 'cram', 'sam')
PR_SNP_NAME = 'per_read_snps'
SNP_NAME = 'snps'
PR_MOD_NAME = 'per_read_mods'
MOD_NAME = 'mods'
ALIGN_OUTPUTS = set((MAP_NAME, PR_SNP_NAME, SNP_NAME, PR_MOD_NAME, MOD_NAME))
OUTPUT_FNS = {
    BC_NAME:'basecalls',
    BC_MODS_NAME:'basecalls.modified_base_scores.hdf5',
    MAP_NAME:['mappings', 'mappings.summary.txt'],
    PR_SNP_NAME:['per_read_snp_calls.db',
                 'per_read_snp_calls.txt'],
    SNP_NAME:'snps.vcf',
    PR_MOD_NAME:['per_read_modified_base_calls.db',
                 'per_read_modified_base_calls.txt'],
    MOD_NAME:'mods.mvcf'
}
LOG_FILENAME = 'log.txt'
COMP_BASES = dict(zip(map(ord, 'ACGT'), map(ord, 'TGCA')))

_MAX_QUEUE_SIZE = 1000


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
