from tqdm import tqdm

from megalodon import variants
from ._extras_parsers import get_parser_phase_variants_whatshap_filter


def is_complex_variant(ref, alts):
    # single base swaps aren't complex
    if any(len(allele) > 1 for allele in alts + [ref]):
        for alt in alts:
            simp_ref, simp_alt, _, _ = variants.simplify_var_seq(ref, alt)
            # if an allele simplifies to a SNV continue
            if len(simp_ref) == 0 and len(simp_alt) == 0:
                continue
            # if simplified sequence does not leave either allele empty
            # then this is a complex variant which cannot be processed by
            # whatshap
            if len(simp_ref) > 0 and len(simp_alt) > 0:
                return True
    return False


def get_qual(vcf_line):
    qual = vcf_line.split()[5]
    try:
        qual = int(qual)
    except ValueError:
        qual = 0
    return qual


def get_pos_ref_alts(vcf_line):
    chrm, pos, _, ref, alts = vcf_line.split()[:5]
    return chrm, int(pos), ref, alts.split(',')


def _main(args):
    out_fp = open(args.out_vcf, 'w')
    filt_fp = None if args.filtered_records is None else open(
        args.filtered_records, 'w')

    with open(args.in_vcf) as fp:
        prev_line = prev_chrm = prev_end = None
        for line in tqdm(fp, desc='Filtering VCF', unit=' lines', smoothing=0):
            if line.startswith('#'):
                out_fp.write(line)
                continue
            chrm, start, ref, alts = get_pos_ref_alts(line)
            # skip complex variants
            if is_complex_variant(ref, alts):
                if filt_fp is not None:
                    filt_fp.write('COMLEX_VARIANT: ' + line)
                continue

            if prev_chrm == chrm and prev_end > start:
                if get_qual(line) > get_qual(prev_line):
                    if filt_fp is not None:
                        filt_fp.write('OVERLAPPING_VARIANT: ' + prev_line)
                    prev_line = line
                else:
                    if filt_fp is not None:
                        filt_fp.write('OVERLAPPING_VARIANT: ' + line)
                continue

            if prev_line is not None:
                out_fp.write(line)
            prev_line = line
            prev_chrm = chrm
            prev_end = start + len(ref)

    if prev_line is not None:
        out_fp.write(line)
    if filt_fp is None:
        filt_fp.close()

    return


if __name__ == '__main__':
    _main(get_parser_phase_variants_whatshap_filter().parse_args())
