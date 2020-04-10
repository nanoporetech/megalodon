import sys
import argparse

import pysam
import numpy as np

from megalodon import variants, megalodon_helper as mh


HEADER = """##fileformat=VCFv4.1
##source=megalodon_haploid_merge
{}
##phasing=megalodon_haploid_merge
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Log10 likelihoods for genotypes">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
"""

CONTIG_HEADER_LINE = "##contig=<ID={},length={}>"

RECORD_LINE = ('{chrm}\t{pos}\t{rid}\t{ref}\t{alts}\t{qual}\t.\tDP={dp:d}\t' +
               'GT:GQ:DP:GL:PL\t{gt}:{gq:.0f}:{dp:d}:{gl}:{pl}\n')


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'diploid_called_variants',
        help='Phased variants from which the diploid calls are derived.')
    parser.add_argument(
        'haplotype1_variants',
        help='Variant file for haplotype 1.')
    parser.add_argument(
        'haplotype2_variants',
        help='Variant file for haplotype 2.')
    parser.add_argument(
        '--out-vcf', default='merged_haploid_variants.vcf',
        help='Output name for VCF. Default: %(default)s')

    return parser


def are_same_var(v0, v1, v2):
    if v1 is None or v2 is None:
        return False
    return (v0.chrom == v1.chrom == v2.chrom and
            v0.pos == v1.pos == v2.pos and
            v0.ref == v1.ref == v2.ref and
            len(v0.alts) == len(v1.alts) == len(v2.alts) and
            all((v0i == v1i == v2i)
                for v0i, v1i, v2i in zip(v0.alts, v1.alts, v2.alts)))


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
    alleles = [var.ref, ] + list(var.alts)
    gt_i = 0
    min_homo_pl = pls[0]
    min_homo_pl_gt = '{0}|{0}'.format(alleles[0])
    for a2 in range(len(alleles)):
        for a1 in range(a2 + 1):
            if a1 == a2 and pls[gt_i] < min_homo_pl:
                min_homo_pl = pls[gt_i]
                min_homo_pl_gt = '{0}|{0}'.format(alleles[a1])
            gt_i += 1
    return min_homo_pl_gt


def write_var(curr_v0_rec, curr_v1_rec, curr_v2_rec, out_vars, contig):
    s0_attrs = next(iter(curr_v0_rec.samples.values()))
    gt0, gq0, gl0, pl0, dp = (
        s0_attrs['GT'], s0_attrs['GQ'], s0_attrs['GL'], s0_attrs['PL'],
        int(s0_attrs['DP']))
    pos, rid, ref, alts = (curr_v0_rec.pos, curr_v0_rec.id, curr_v0_rec.ref,
                           ','.join(curr_v0_rec.alts))
    if are_same_var(curr_v0_rec, curr_v1_rec, curr_v2_rec):
        s1_attrs = next(iter(curr_v1_rec.samples.values()))
        s2_attrs = next(iter(curr_v2_rec.samples.values()))
        # don't let homozygous calls change a heterozygous call
        # TODO explore settings where this restriction could
        # be relaxed
        if len(set(gt0)) == 1:
            gt = '{}|{}'.format(*gt0)
            gq, gl, pl = gq0, gl0, pl0
            qual = parse_qual(curr_v0_rec.qual)
        else:
            gt, gq, gl, pl, qual = compute_diploid_stats(
                s1_attrs['GL'], s2_attrs['GL'])
    else:
        # convert genotype to most likely homozygous call. since whatshap
        # could not phase this site, it is likely not heterozygous
        qual = parse_qual(curr_v0_rec.qual)
        gt = get_most_likely_homo(curr_v0_rec, pl0)
        gq, gl, pl = gq0, gl0, pl0

    if rid is None:
        rid = '.'
    qual = '.' if qual == 0 else '{:d}'.format(qual)
    gl_fmt = ','.join('{:.2f}' for _ in range(len(gl))).format(*gl)
    pl_fmt = ','.join('{:.0f}' for _ in range(len(pl))).format(*pl)
    out_vars.write(RECORD_LINE.format(
        chrm=contig, pos=pos, rid=rid, ref=ref, alts=alts, qual=qual, dp=dp,
        gt=gt, gq=gq, gl=gl_fmt, pl=pl_fmt))

    return


def iter_contig_vars(vars0_contig_iter, vars1_contig_iter, vars2_contig_iter):
    """ Iterate variants from source and both haplotypes.
    """
    def next_or_none(vars_iter):
        """ Extract next variant or return None if iterator is exhausted
        """
        try:
            return next(vars_iter)
        except StopIteration:
            return None

    def get_uniq_pos(var):
        return var.pos, var.ref, var.alts

    def init_curr_vars(next_rec):
        return dict([(get_uniq_pos(next_rec), next_rec)]) \
            if next_rec is not None else {(-1, ): None}

    def get_pos_vars(curr_recs, next_rec, vars_iter, v0_vars=None):
        """ Extract all variants with the same position from variant iterator
        """
        curr_pos = list(curr_recs.keys())[0][0]
        if v0_vars is not None:
            v0_pos = list(v0_vars.keys())[0][0]
            if v0_pos < curr_pos:
                return curr_recs, next_rec
            # check that haplotype variants occur in source VCF
            while curr_pos < v0_pos:
                sys.stderr.write(
                    'WARNING: Variant found in haplotype file which is ' +
                    'missing from source file.\n')
                curr_recs = init_curr_vars(next_rec)
                next_rec = next_or_none(vars_iter)
                curr_pos = list(curr_recs.keys())[0][0]
        while next_rec is not None and next_rec.pos == curr_pos:
            curr_recs[get_uniq_pos(next_rec)] = next_rec
            next_rec = next_or_none(vars_iter)
        return curr_recs, next_rec

    # initialize first contig variants
    first_v0_rec = next_or_none(vars0_contig_iter)
    # Multiple variants may occur at the same position. To ensure these
    # variants are grouped correctly no matter the sorted order based on
    # alleles, record all variants at each position before returning
    curr_v0_recs = dict([(get_uniq_pos(first_v0_rec), first_v0_rec), ])
    next_v0_rec = next_or_none(vars0_contig_iter)

    first_v1_rec = next_or_none(vars1_contig_iter)
    first_v2_rec = next_or_none(vars2_contig_iter)
    # One of the files has no variants on this contig
    if first_v1_rec is None or first_v2_rec is None:
        while next_v0_rec is not None:
            for pos in curr_v0_recs:
                yield curr_v0_recs[pos], None, None
            curr_v0_recs, next_v0_rec = get_pos_vars(
                curr_v0_recs, next_v0_rec, vars0_contig_iter)
        for pos in curr_v0_recs:
            yield curr_v0_recs[pos], None, None
        return

    # initialize haplotype current and next variants
    curr_v1_recs = init_curr_vars(first_v1_rec)
    curr_v2_recs = init_curr_vars(first_v2_rec)
    next_v1_rec = next_or_none(vars1_contig_iter)
    next_v2_rec = next_or_none(vars2_contig_iter)

    while next_v0_rec is not None:
        # get all variants at the current source VCF position
        curr_v0_recs, next_v0_rec = get_pos_vars(
            curr_v0_recs, next_v0_rec, vars0_contig_iter)
        curr_v1_recs, next_v1_rec = get_pos_vars(
            curr_v1_recs, next_v1_rec, vars1_contig_iter, curr_v0_recs)
        curr_v2_recs, next_v2_rec = get_pos_vars(
            curr_v2_recs, next_v2_rec, vars2_contig_iter, curr_v0_recs)
        for pos in curr_v0_recs:
            try:
                # variant may not be in the haplotype files
                v1_pos_var = curr_v1_recs[pos]
                v2_pos_var = curr_v2_recs[pos]
            except KeyError:
                v1_pos_var = v2_pos_var = None
            yield curr_v0_recs[pos], v1_pos_var, v2_pos_var
        curr_v0_recs = init_curr_vars(next_v0_rec)
        curr_v1_recs = init_curr_vars(next_v1_rec)
        curr_v2_recs = init_curr_vars(next_v2_rec)

    for pos in curr_v0_recs:
        try:
            # variant may not be in the haplotype files
            v1_pos_var = curr_v1_recs[pos]
            v2_pos_var = curr_v2_recs[pos]
        except KeyError:
            v1_pos_var = v2_pos_var = None
        yield curr_v0_recs[pos], v1_pos_var, v2_pos_var


def get_contig_iter(vars_idx, contig):
    try:
        return iter(vars_idx.fetch(contig))
    except ValueError:
        return iter([])


def main():
    args = get_parser().parse_args()

    sys.stderr.write('Openning VCF files.')
    vars0_idx = pysam.VariantFile(args.diploid_called_variants)
    vars1_idx = pysam.VariantFile(args.haplotype1_variants)
    vars2_idx = pysam.VariantFile(args.haplotype2_variants)
    try:
        contigs0 = list(vars0_idx.header.contigs.keys())
        vars0_idx.fetch(contigs0[0])
        vars1_idx.fetch(next(iter(vars1_idx.header.contigs.keys())))
        vars2_idx.fetch(next(iter(vars2_idx.header.contigs.keys())))
    except ValueError:
        raise mh.MegaError(
            'Variant files must be indexed. Use bgzip and tabix.')

    sys.stderr.write('Processing variants.')
    out_vars = open(args.out_vcf, 'w')
    out_vars.write(HEADER.format('\n'.join((CONTIG_HEADER_LINE.format(
        ctg.name, ctg.length) for ctg in vars0_idx.header.contigs.values()))))
    for contig in contigs0:
        for curr_v0_rec, curr_v1_rec, curr_v2_rec in iter_contig_vars(
                get_contig_iter(vars0_idx, contig),
                get_contig_iter(vars1_idx, contig),
                get_contig_iter(vars2_idx, contig)):
            write_var(curr_v0_rec, curr_v1_rec, curr_v2_rec, out_vars, contig)

    out_vars.close()

    variants.index_variants(args.out_vcf)


if __name__ == '__main__':
    main()
