from collections import namedtuple

import pysam
from tqdm import tqdm

from megalodon import logging, megalodon_helper as mh
from ._extras_parsers import get_parser_modified_bases_create_motif_bed


LOGGER = logging.get_logger()

MOTIF_INFO = namedtuple('MOTIF_INFO', (
    'bases_before', 'bases_after', 'raw_motif', 'motif', 'rc_motif'))
BED_TMPLT = '{chrom}\t{pos}\t{end}\t.\t.\t{strand}\n'


def parse_motifs(raw_motifs):
    motifs = []
    for raw_motif, bases_before in raw_motifs:
        bases_before = int(bases_before)
        bases_after = len(raw_motif) - bases_before - 1
        motif = mh.compile_motif_pat(raw_motif)
        rc_motif = mh.compile_rev_comp_motif_pat(raw_motif)
        motifs.append(MOTIF_INFO(
            bases_before=bases_before, bases_after=bases_after,
            raw_motif=raw_motif, motif=motif, rc_motif=rc_motif))

    return motifs


def _main(args):
    logging.init_logger()

    # parse motifs
    motifs = parse_motifs(args.motif)
    # open indexed FASTA reference
    ref = pysam.FastaFile(args.reference)

    with open(args.out_filename, 'w') as fp:
        # sort using RefName
        for chrm in tqdm(sorted([mh.RefName(chrm) for chrm in ref.references]),
                         desc='Contigs', smoothing=0, dynamic_ncols=True):
            chrm_seq = ref.fetch(chrm)
            chrm_sites = []
            for motif in motifs:
                for m in motif.motif.finditer(chrm_seq):
                    pos = m.start() + motif.bases_before
                    chrm_sites.append((pos, '+'))
                for m in motif.rc_motif.finditer(chrm_seq):
                    pos = m.start() + motif.bases_after
                    chrm_sites.append((pos, '-'))
            fp.write(''.join(BED_TMPLT.format(
                chrom=chrm, pos=pos, end=pos + 1, strand=strand)
                for pos, strand in sorted(chrm_sites)))


if __name__ == '__main__':
    _main(get_parser_modified_bases_create_motif_bed().parse_args())
