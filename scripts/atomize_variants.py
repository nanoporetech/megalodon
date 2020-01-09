import sys
import argparse

from megalodon import mapping, megalodon_helper as mh, variants

HEADER = """##fileformat=VCFv4.1
##source=megalodon_atomized
{}
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
"""
CONTIG_HEADER_LINE = "##contig=<ID={},length={}>"
RECORD_LINE = ('{chrm}\t{pos}\t{rid}\t{ref}\t{alts}\t.\t.\t.\t.\t.\n')

def get_parser():
    parser = argparse.ArgumentParser(
        description='Atomize variants so this does not have to be ' +
        'completed during read processing')
    parser.add_argument('in_vcf', help='Proposed varitants (VCF)')
    parser.add_argument(
        'reference',
        help='Reference FASTA or minimap2 index file corresponding to VCF.')
    parser.add_argument(
        '--out-vcf', default='atomized_variants.megalodon.vcf',
        help='Output VCF file. Default: %(default)s')
    parser.add_argument(
        '--max-indel-size', type=int, default=50,
        help='Maximum difference in number of reference and alternate ' +
        'bases. Default: %(default)d')
    return parser

def main():
    args = get_parser().parse_args()
    sys.stderr.write('Loading reference\n')
    aligner = mapping.alignerPlus(
        str(args.reference), preset=str('map-ont'), best_n=1)
    aligner.add_ref_lens()
    sys.stderr.write('Loading variants\n')
    var_data = variants.VarData(
        args.in_vcf, args.max_indel_size, keep_var_fp_open=True,
        aligner=aligner)
    contigs = var_data.variants_idx.header.contigs.values()
    sys.stderr.write('Atomizing variants\n')
    with open(args.out_vcf, 'w') as out_vars:
        out_vars.write(HEADER.format('\n'.join((
            CONTIG_HEADER_LINE.format(ctg.name, ctg.length)
            for ctg in contigs))))
        for ctg in contigs:
            chrm_seq = aligner.seq(ctg.name)
            if len(chrm_seq) != ctg.length:
                sys.stderr.write((
                    'WARNING: Mismatched contig lengths ({}) between ' +
                    'reference ({}) and input VCF ({})\n').format(
                        ctg.name, len(chrm_seq), ctg.length))
            map_pos = mapping.MAP_POS(
                chrm=ctg.name, strand=None, start=0, end=len(chrm_seq),
                q_trim_start=None, q_trim_end=None)
            for var in var_data.fetch_read_variants(
                    map_pos, mh.seq_to_int(chrm_seq)):
                out_vars.write(RECORD_LINE.format(
                    chrm=ctg.name, pos=var.ref_start + 1, rid=var.id,
                    ref=var.ref, alts=','.join(var.alts)))

    sys.stderr.write('Indexing output variant file\n')
    variants.index_variants(args.out_vcf)

    return

if __name__ == '__main__':
    main()
