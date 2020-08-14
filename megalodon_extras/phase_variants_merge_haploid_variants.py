import sys

import pysam
import numpy as np
from tqdm import tqdm

from megalodon import logging, megalodon_helper as mh, variants
from ._extras_parsers import get_parser_phase_variants_merge_haploid_variants


LOGGER = logging.get_logger()

_HEADER_LINES = (
    '##fileformat=VCFv4.1',
    '##source=megalodon_haploid_merge',
    '{}',
    '##phasing=megalodon_haploid_merge',
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
    '##FORMAT=<ID=GL,Number=G,Type=Float,Description="Log10 likelihoods ' +
    'for genotypes">',
    '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, ' +
    'Phred-scaled likelihoods for genotypes">',
    '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	' +
    'FORMAT	SAMPLE')
HEADER = '\n'.join(_HEADER_LINES)

CONTIG_HEADER_LINE = "##contig=<ID={},length={}>"

RECORD_LINE = ('{chrm}\t{pos}\t{rid}\t{ref}\t{alts}\t{qual}\t.\tDP={dp:d}\t' +
               'GT:GQ:DP:GL:PL\t{gt}:{gq:.0f}:{dp:d}:{gl}:{pl}\n')


def are_same_var(s_v, h1_v, h2_v):
    if h1_v is None or h2_v is None:
        return False
    return (s_v.chrom == h1_v.chrom == h2_v.chrom and
            s_v.pos == h1_v.pos == h2_v.pos and
            s_v.ref == h1_v.ref == h2_v.ref and
            len(s_v.alts) == len(h1_v.alts) == len(h2_v.alts) and
            all((s_vi == h1_vi == h2_vi)
                for s_vi, h1_vi, h2_vi in zip(s_v.alts, h1_v.alts, h2_v.alts)))


def parse_qual(qual):
    if qual is None:
        return 0
    return int(qual)


def compute_diploid_stats(gl1, gl2):
    probs, gts = [], []
    for a2 in range(len(gl1)):
        for a1 in range(a2 + 1):
            if a1 == a2:
                gts.append('{}|{}'.format(a1, a2))
                probs.append(np.exp(gl1[a1] + gl2[a1]))
                continue
            p12 = np.exp(gl1[a1] + gl2[a2])
            p21 = np.exp(gl1[a2] + gl2[a1])
            if p12 > p21:
                gts.append('{}|{}'.format(a1, a2))
                probs.append(p12)
            else:
                gts.append('{}|{}'.format(a2, a1))
                probs.append(p21)
    with np.errstate(divide='ignore'):
        gl = np.maximum(mh.MIN_GL_VALUE, np.log10(probs))
    raw_pl = -10 * gl
    # "normalized" PL values stored as decsribed by VCF format
    pl = np.minimum(raw_pl - raw_pl.min(), mh.MAX_PL_VALUE)
    s_pl = np.sort(pl)

    gq = np.around(s_pl[1])
    gt = gts[np.argmax(probs)]
    try:
        qual = int(np.minimum(np.around(raw_pl[0]), mh.MAX_PL_VALUE))
    except ValueError:
        qual = mh.MAX_PL_VALUE

    return gt, gq, gl, pl, qual


def get_most_likely_homo(var, pls):
    gt_i = 0
    min_homo_pl = pls[0]
    min_homo_pl_gt = '0|0'
    for a2 in range(len(var.alts) + 1):
        for a1 in range(a2 + 1):
            if a1 == a2 and pls[gt_i] < min_homo_pl:
                min_homo_pl = pls[gt_i]
                min_homo_pl_gt = '{0}|{0}'.format(a1)
            gt_i += 1
    return min_homo_pl_gt


def write_var(curr_s_rec, curr_h1_rec, curr_h2_rec, out_vars, contig):
    s0_attrs = next(iter(curr_s_rec.samples.values()))
    gt0, gq0, gl0, pl0, dp = (
        s0_attrs['GT'], s0_attrs['GQ'], s0_attrs['GL'], s0_attrs['PL'],
        int(s0_attrs['DP']))
    pos, rid, ref, alts = (curr_s_rec.pos, curr_s_rec.id, curr_s_rec.ref,
                           ','.join(curr_s_rec.alts))
    if are_same_var(curr_s_rec, curr_h1_rec, curr_h2_rec):
        s1_attrs = next(iter(curr_h1_rec.samples.values()))
        s2_attrs = next(iter(curr_h2_rec.samples.values()))
        # don't let homozygous calls change a heterozygous call
        # TODO explore settings where this restriction could
        # be relaxed
        if len(set(gt0)) == 1:
            gt = '{}|{}'.format(*gt0)
            gq, gl, pl = gq0, gl0, pl0
            qual = parse_qual(curr_s_rec.qual)
        else:
            gt, gq, gl, pl, qual = compute_diploid_stats(
                s1_attrs['GL'], s2_attrs['GL'])
    else:
        # convert genotype to most likely homozygous call. since whatshap
        # could not phase this site, it is likely not heterozygous
        qual = parse_qual(curr_s_rec.qual)
        gt = get_most_likely_homo(curr_s_rec, pl0)
        gq, gl, pl = gq0, gl0, pl0

    if rid is None:
        rid = '.'
    qual = '.' if qual == 0 else '{:d}'.format(qual)
    gl_fmt = ','.join('{:.2f}' for _ in range(len(gl))).format(*gl)
    pl_fmt = ','.join('{:.0f}' for _ in range(len(pl))).format(*pl)
    out_vars.write(RECORD_LINE.format(
        chrm=contig, pos=pos, rid=rid, ref=ref, alts=alts, qual=qual, dp=dp,
        gt=gt, gq=gq, gl=gl_fmt, pl=pl_fmt))


def iter_contig_vars(
        source_var_iter, h1_var_iter, h2_var_iter, contig, bar,
        force_invalid_vars):
    """ Iterate variants from source and both haplotypes.
    """
    def next_or_none(var_iter):
        """ Extract next variant or return None if iterator is exhausted
        """
        try:
            return next(var_iter)
        except StopIteration:
            return None

    def get_uniq_pos(var):
        return var.pos, var.ref, var.alts

    def init_curr_vars(next_rec):
        # Multiple variants may occur at the same position. To ensure these
        # variants are grouped correctly no matter the sorted order based on
        # alleles, record all variants at each position before returning
        return dict([(get_uniq_pos(next_rec), next_rec)]) \
            if next_rec is not None else {(-1, ): None}

    def get_pos_vars(curr_recs, next_rec, var_iter, source_pos_vars=None):
        """ Extract all variants with the same position from variant iterator
        """
        curr_pos = list(curr_recs.keys())[0][0]
        if source_pos_vars is not None:
            s_pos = list(source_pos_vars.keys())[0][0]
            # source variant not found in haplotype VCF
            if s_pos < curr_pos:
                return curr_recs, next_rec
            # check that haplotype variants occur in source VCF
            while next_rec is not None and curr_pos < s_pos:
                if not force_invalid_vars:
                    bar.close()
                    LOGGER.error((
                        'Variant found in haplotype file which is missing ' +
                        'from source file: {}:{}.\nSet ' +
                        '--force-invalid-variant-processing to force ' +
                        'processing. Results when invalid variants are ' +
                        'encountered is not defined.').format(
                            contig, curr_pos))
                    sys.exit(1)
                bar.write(
                    'WARNING: Variant found in haplotype file which is ' +
                    'missing from source file: {}:{}.'.format(
                        contig, curr_pos))
                curr_recs = init_curr_vars(next_rec)
                next_rec = next_or_none(var_iter)
                curr_pos = list(curr_recs.keys())[0][0]
        while next_rec is not None and next_rec.pos == curr_pos:
            curr_recs[get_uniq_pos(next_rec)] = next_rec
            next_rec = next_or_none(var_iter)
        return curr_recs, next_rec

    # initialize next contig variants
    next_s_rec = next_or_none(source_var_iter)
    if next_s_rec is None:
        return
    curr_s_recs = init_curr_vars(next_s_rec)
    next_s_rec = next_or_none(source_var_iter)

    # initialize haplotype contig variants
    next_h1_rec = next_or_none(h1_var_iter)
    next_h2_rec = next_or_none(h2_var_iter)
    # One of the files has no variants on this contig
    if next_h1_rec is None or next_h2_rec is None:
        while next_s_rec is not None:
            curr_s_recs = init_curr_vars(next_s_rec)
            next_s_rec = next_or_none(source_var_iter)
            curr_s_recs, next_s_rec = get_pos_vars(
                curr_s_recs, next_s_rec, source_var_iter)
            for pos in curr_s_recs:
                yield curr_s_recs[pos], None, None
        for pos in curr_s_recs:
            yield curr_s_recs[pos], None, None
        return
    curr_h1_recs = init_curr_vars(next_h1_rec)
    curr_h2_recs = init_curr_vars(next_h2_rec)
    next_h1_rec = next_or_none(h1_var_iter)
    next_h2_rec = next_or_none(h2_var_iter)

    while next_s_rec is not None:
        # get all variants at the current source VCF position
        curr_s_recs, next_s_rec = get_pos_vars(
            curr_s_recs, next_s_rec, source_var_iter)
        curr_h1_recs, next_h1_rec = get_pos_vars(
            curr_h1_recs, next_h1_rec, h1_var_iter, curr_s_recs)
        curr_h2_recs, next_h2_rec = get_pos_vars(
            curr_h2_recs, next_h2_rec, h2_var_iter, curr_s_recs)
        for pos in curr_s_recs:
            try:
                # variant may not be in the haplotype files
                h1_pos_var = curr_h1_recs[pos]
                h2_pos_var = curr_h2_recs[pos]
            except KeyError:
                h1_pos_var = h2_pos_var = None
            yield curr_s_recs[pos], h1_pos_var, h2_pos_var

        # start new position vars with next and get new next variant
        s_prev_pos = list(curr_s_recs.keys())[0][0]
        curr_s_recs = init_curr_vars(next_s_rec)
        next_s_rec = next_or_none(source_var_iter)
        # if the current haplotype variants matched the source variant get a
        # new haplotype variant, else these variants have not been processed
        if list(curr_h1_recs.keys())[0][0] == s_prev_pos:
            curr_h1_recs = init_curr_vars(next_h1_rec)
            next_h1_rec = next_or_none(h1_var_iter)
        if list(curr_h2_recs.keys())[0][0] == s_prev_pos:
            curr_h2_recs = init_curr_vars(next_h2_rec)
            next_h2_rec = next_or_none(h2_var_iter)

    # process final position
    for pos in curr_s_recs:
        try:
            # variant may not be in the haplotype files
            h1_pos_var = curr_h1_recs[pos]
            h2_pos_var = curr_h2_recs[pos]
        except KeyError:
            h1_pos_var = h2_pos_var = None
        # if get_pos_vars exhausts contig variants then curr_s_recs will be
        # empty
        if curr_s_recs[pos] is not None:
            yield curr_s_recs[pos], h1_pos_var, h2_pos_var


def get_contig_iter(vars_idx, contig):
    try:
        return iter(vars_idx.fetch(contig))
    except ValueError:
        return iter([])


def _main(args):
    logging.init_logger()
    LOGGER.info('Opening VCF files.')
    source_vars = pysam.VariantFile(args.diploid_called_variants)
    h1_vars = pysam.VariantFile(args.haplotype1_variants)
    h2_vars = pysam.VariantFile(args.haplotype2_variants)
    try:
        contigs0 = list(source_vars.header.contigs.keys())
        source_vars.fetch(contigs0[0])
        h1_vars.fetch(next(iter(h1_vars.header.contigs.keys())))
        h2_vars.fetch(next(iter(h2_vars.header.contigs.keys())))
    except ValueError:
        raise mh.MegaError(
            'Variant files must be indexed. Use bgzip and tabix.')

    LOGGER.info('Processing variants.')
    out_vars = open(args.out_vcf, 'w')
    out_vars.write(HEADER.format('\n'.join(
        (CONTIG_HEADER_LINE.format(ctg.name, ctg.length)
         for ctg in source_vars.header.contigs.values()))))
    bar = tqdm(total=len(contigs0), smoothing=0, unit=' contigs',
               dynamic_ncols=True, desc='Variant Processing', mininterval=0)
    for contig in contigs0:
        for curr_s_rec, curr_h1_rec, curr_h2_rec in tqdm(
                iter_contig_vars(get_contig_iter(source_vars, contig),
                                 get_contig_iter(h1_vars, contig),
                                 get_contig_iter(h2_vars, contig),
                                 contig, bar,
                                 args.force_invalid_variant_processing),
                smoothing=0, unit=' variants', dynamic_ncols=True, leave=False,
                desc='{} Variants'.format(contig)):
            write_var(curr_s_rec, curr_h1_rec, curr_h2_rec, out_vars, contig)
        bar.update(1)

    out_vars.close()

    variants.index_variants(args.out_vcf)


if __name__ == '__main__':
    _main(get_parser_phase_variants_merge_haploid_variants().parse_args())
