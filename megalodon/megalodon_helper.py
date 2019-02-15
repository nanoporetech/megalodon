import numpy as np

ALPHABET = 'ACGT'
BC_NAME = 'basecalls'
BC_OUT_FMTS = ('fasta',)
MAP_NAME = 'mappings'
MAP_OUT_FMTS = ('bam', 'cram', 'sam')
PR_SNP_NAME = 'per_read_snps'
PR_MOD_NAME = 'per_read_mods'
ALIGN_OUTPUTS = set((MAP_NAME, PR_SNP_NAME, PR_MOD_NAME))
OUTPUT_FNS = {
    BC_NAME:'basecalls',
    MAP_NAME:['mappings', 'mappings.summary.txt'],
    PR_SNP_NAME:'per_read_snp_calls.txt',
    PR_MOD_NAME:'per_read_modified_base_calls.txt',
}
COMP_BASES = dict(zip(map(ord, 'ACGT'), map(ord, 'TGCA')))


class MegaError(Exception):
    """ Custom megalodon error for more graceful error handling
    """
    pass

def nstate_to_nbase(nstate):
    return int(np.sqrt(0.25 + (0.5 * nstate)) - 0.5)

def revcomp(seq):
    return seq.translate(COMP_BASES)[::-1]


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
