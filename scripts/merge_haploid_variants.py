import argparse

import pysam
import numpy as np

from megalodon import snps


HEADER = """##fileformat=VCFv4.1
##source=megalodon_haploid_merge
{}
##phasing=megalodon_haploid_merge
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
"""

CONTIG_HEADER_LINE = "##contig=<ID={},length={}>"

RECORD_LINE = "{}\t{}\t{}\t{}\t{}\t{}\t.\t.\tGT\t{}\n"


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'phased_variants',
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
    return qual

def main():
    args = get_parser().parse_args()

    vars0_idx = pysam.VariantFile(args.phased_variants)
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
                gt0 = next(iter(curr_v0_rec.samples.values()))['GT']
                if are_same_var(curr_v0_rec, curr_v1_rec, curr_v2_rec):
                    s1_attrs = next(iter(curr_v1_rec.samples.values()))
                    s2_attrs = next(iter(curr_v2_rec.samples.values()))
                    # don't let haploid calls change a homozygous call
                    # TODO explore settings where this restriction could
                    # be relaxed
                    if len(set(gt0)) == 1:
                        gt = '{}|{}'.format(*gt0)
                        qual = parse_qual(curr_v0_rec.qual)
                    else:
                        gt1 = s1_attrs['GT'][0]
                        gt2 = s2_attrs['GT'][0]
                        gt = '{}|{}'.format(gt1, gt2)
                        if gt1 != 0 and gt2 == 0:
                            qual = parse_qual(curr_v1_rec.qual)
                        elif gt1 == 0 and gt2 != 0:
                            qual = parse_qual(curr_v2_rec.qual)
                        else:
                            qual = max(parse_qual(curr_v1_rec.qual),
                                       parse_qual(curr_v2_rec.qual))
                    if qual == 0: qual = '.'
                    out_vars.write(RECORD_LINE.format(
                        contig, curr_v1_rec.pos, curr_v1_rec.id,
                        curr_v1_rec.ref, ','.join(curr_v1_rec.alts), qual, gt))
                    curr_v0_rec = next(vars0_contig_iter)
                    curr_v1_rec = next(vars1_contig_iter)
                    curr_v2_rec = next(vars2_contig_iter)
                elif curr_v1_rec.pos < curr_v0_rec.pos:
                    # variant in haplotype 1 does not exist in phased variants
                    # this should never happen
                    curr_v1_rec = next(vars1_contig_iter)
                elif curr_v2_rec.pos < curr_v0_rec.pos:
                    # variant in haplotype 2 does not exist in phased variants
                    # this should never happen
                    curr_v2_rec = next(vars2_contig_iter)
                else:
                    # write phased variant back out.
                    if curr_v0_rec.qual is None or curr_v0_rec.qual == 0:
                        qual = '.'
                    else:
                        qual = int(curr_v0_rec.qual)
                    out_vars.write(RECORD_LINE.format(
                        contig, curr_v0_rec.pos, curr_v0_rec.id,
                        curr_v0_rec.ref, ','.join(curr_v0_rec.alts),
                        qual, '{}|{}'.format(*gt0)))
                    if curr_v0_rec.pos == curr_v1_rec.pos:
                        curr_v1_rec = next(vars1_contig_iter)
                    if curr_v0_rec.pos == curr_v2_rec.pos:
                        curr_v2_rec = next(vars2_contig_iter)
                    curr_v0_rec = next(vars0_contig_iter)
            except StopIteration:
                break

    out_vars.close()

    index_var_fn = snps.index_variants(args.out_vcf)

    return


if __name__ == '__main__':
    main()
