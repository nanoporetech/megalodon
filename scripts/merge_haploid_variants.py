import argparse

import pysam
import numpy as np

from megalodon import snps, megalodon_helper as mh


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
        help='Variant file for haplotype 1.')
    parser.add_argument(
        '--out-vcf', default='merged_haploid_variants.vcf',
        help='Output name for VCF. Default: %(default)s')

    return parser

def are_same_var(v0, v1, v2):
    return (v0.chrom == v1.chrom == v2.chrom and
            v0.pos == v1.pos == v2.pos and
            v0.ref == v1.ref == v2.ref and
            len(v0.alts) == len(v1.alts) == len(v2.alts) and
            all((v1i == v1i == v2i)
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

def write_var(
        curr_v0_rec, curr_v1_rec, curr_v2_rec, out_vars,
        vars0_contig_iter, vars1_contig_iter, vars2_contig_iter, contig):
    s0_attrs = next(iter(curr_v0_rec.samples.values()))
    gt0, gq0, gl0, pl0, dp = (
        s0_attrs['GT'], s0_attrs['GQ'], s0_attrs['GL'], s0_attrs['PL'],
        int(s0_attrs['DP']))
    pos, rid, ref, alts = (curr_v0_rec.pos, curr_v0_rec.id, curr_v0_rec.ref,
                           ','.join(curr_v0_rec.alts))
    if are_same_var(curr_v0_rec, curr_v1_rec, curr_v2_rec):
        s1_attrs = next(iter(curr_v1_rec.samples.values()))
        s2_attrs = next(iter(curr_v2_rec.samples.values()))
        # don't let haploid calls change a homozygous call
        # TODO explore settings where this restriction could
        # be relaxed
        if len(set(gt0)) == 1:
            gt = '{}|{}'.format(*gt0)
            gq, gl, pl = gq0, gl0, pl0
            qual = parse_qual(curr_v0_rec.qual)
        else:
            gt, gq, gl, pl, qual = compute_diploid_stats(
                s1_attrs['GL'], s2_attrs['GL'])
        curr_v0_rec = next(vars0_contig_iter)
        curr_v1_rec = next(vars1_contig_iter)
        curr_v2_rec = next(vars2_contig_iter)
    else:
        # write un-phased variant back out.
        qual = parse_qual(curr_v0_rec.qual)
        gt = '{}|{}'.format(*gt0)
        gq, gl, pl = gq0, gl0, pl0
        if curr_v0_rec.pos == curr_v1_rec.pos:
            curr_v1_rec = next(vars1_contig_iter)
        if curr_v0_rec.pos == curr_v2_rec.pos:
            curr_v2_rec = next(vars2_contig_iter)
        curr_v0_rec = next(vars0_contig_iter)

    qual = '.' if qual == 0 else '{:d}'.format(qual)
    gl_fmt = ','.join('{:.2f}' for _ in range(len(gl))).format(*gl)
    pl_fmt = ','.join('{:.0f}' for _ in range(len(pl))).format(*pl)
    out_vars.write(RECORD_LINE.format(
        chrm=contig, pos=pos, rid=rid, ref=ref, alts=alts, qual=qual, dp=dp,
        gt=gt, gq=gq, gl=gl_fmt, pl=pl_fmt))

    return curr_v0_rec, curr_v1_rec, curr_v2_rec

def main():
    args = get_parser().parse_args()

    vars0_idx = pysam.VariantFile(args.diploid_called_variants)
    vars1_idx = pysam.VariantFile(args.haplotype1_variants)
    vars2_idx = pysam.VariantFile(args.haplotype2_variants)
    try:
        contigs0 = list(vars0_idx.header.contigs.keys())
        vars0_idx.fetch(next(iter(contigs0)), 0, 0)
        contigs1 = list(vars1_idx.header.contigs.keys())
        vars1_idx.fetch(next(iter(contigs1)), 0, 0)
        contigs2 = list(vars2_idx.header.contigs.keys())
        vars2_idx.fetch(next(iter(contigs2)), 0, 0)
    except ValueError:
        raise mh.MegaError(
            'Variants file must be indexed. Use bgzip and tabix.')

    out_vars = open(args.out_vcf, 'w')
    out_vars.write(HEADER.format('\n'.join((CONTIG_HEADER_LINE.format(
        ctg.name, ctg.length) for ctg in vars0_idx.header.contigs.values()))))
    for contig in set(contigs0).intersection(contigs1).intersection(contigs2):
        vars0_contig_iter = iter(vars0_idx.fetch(contig))
        vars1_contig_iter = iter(vars1_idx.fetch(contig))
        vars2_contig_iter = iter(vars2_idx.fetch(contig))
        try:
            curr_v0_rec = next(vars0_contig_iter)
            curr_v1_rec = next(vars1_contig_iter)
            curr_v2_rec = next(vars2_contig_iter)
        except StopIteration:
            continue
        while True:
            try:
                if curr_v1_rec.pos < curr_v0_rec.pos:
                    # variant in haplotype 1 does not exist in phased variants
                    # this should never happen
                    curr_v1_rec = next(vars1_contig_iter)
                    continue
                elif curr_v2_rec.pos < curr_v0_rec.pos:
                    # variant in haplotype 2 does not exist in phased variants
                    # this should never happen
                    curr_v2_rec = next(vars2_contig_iter)
                    continue

                curr_v0_rec, curr_v1_rec, curr_v2_rec = write_var(
                    curr_v0_rec, curr_v1_rec, curr_v2_rec, out_vars,
                    vars0_contig_iter, vars1_contig_iter, vars2_contig_iter,
                    contig)
            except StopIteration:
                break

    out_vars.close()

    index_var_fn = snps.index_variants(args.out_vcf)

    return


if __name__ == '__main__':
    main()
