import argparse

import pysam
import numpy as np


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
        'vcf1', help='First sorted and indexed haploid variant file.')
    parser.add_argument(
        'vcf2', help='Second sorted and indexed haploid variant file.')
    parser.add_argument(
        '--out-vcf', default='merged_haploid_variants.vcf',
        help='Output name for VCF. Default: %(default)s')

    return parser

def are_same_var(v1, v2):
    return (v1.chrom == v2.chrom and
            v1.pos == v2.pos and
            v1.ref == v2.ref and
            all((v1i == v2i) for v1i, v2i in zip(v1.alts, v2.alts)))

def main():
    args = get_parser().parse_args()

    vars1_idx = pysam.VariantFile(args.vcf1)
    vars2_idx = pysam.VariantFile(args.vcf2)
    try:
        contigs1 = list(vars1_idx.header.contigs.keys())
        vars1_idx.fetch(next(iter(contigs1)), 0, 0)
        contigs2 = list(vars2_idx.header.contigs.keys())
        vars2_idx.fetch(next(iter(contigs2)), 0, 0)
    except ValueError:
        raise mh.MegaError(
            'Variants file must be indexed. Use bgzip and tabix.')

    out_vars = open(args.out_vcf, 'w')
    out_vars.write(HEADER.format('\n'.join((CONTIG_HEADER_LINE.format(
        ctg.name, ctg.length) for ctg in vars1_idx.header.contigs.values()))))
    for contig in set(contigs1).intersection(contigs2):
        vars1_contig_iter = iter(vars1_idx.fetch(contig))
        vars2_contig_iter = iter(vars2_idx.fetch(contig))
        try:
            curr_v1_rec = next(vars1_contig_iter)
            curr_v2_rec = next(vars2_contig_iter)
        except StopIteration:
            continue
        while True:
            if are_same_var(curr_v1_rec, curr_v2_rec):
                gt = '{}|{}'.format(
                    next(iter(curr_v1_rec.samples.values()))['GT'][0],
                    next(iter(curr_v2_rec.samples.values()))['GT'][0])
                qual = max(0, int(np.around(np.mean(
                    (curr_v1_rec.qual, curr_v2_rec.qual)))))
                if qual == 0: qual = '.'
                out_vars.write(RECORD_LINE.format(
                    contig, curr_v1_rec.pos, curr_v1_rec.id, curr_v1_rec.ref,
                    ','.join(curr_v1_rec.alts), qual, gt))
                try:
                    curr_v1_rec = next(vars1_contig_iter)
                    curr_v2_rec = next(vars2_contig_iter)
                except StopIteration:
                    break
            elif curr_v1_rec.pos < curr_v2_rec.pos:
                try:
                    curr_v1_rec = next(vars1_contig_iter)
                except StopIteration:
                    break
            elif curr_v1_rec.pos > curr_v2_rec.pos:
                try:
                    curr_v2_rec = next(vars2_contig_iter)
                except StopIteration:
                    break
            else:
                # not equal but same pos
                try:
                    curr_v1_rec = next(vars1_contig_iter)
                    curr_v2_rec = next(vars2_contig_iter)
                except StopIteration:
                    break

    out_vars.close()

    return


if __name__ == '__main__':
    main()
