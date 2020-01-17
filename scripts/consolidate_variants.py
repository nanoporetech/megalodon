import re
import sys
import argparse
from collections import OrderedDict

import numpy as np
from tqdm import tqdm

from megalodon import megalodon_helper as mh


GT_PAT = re.compile('^(?P<a1>.)(?:(?P<sep>[/|])(?P<a2>.))?$')
INFO_PAT = re.compile('([^;=]+)(?:=([^;]+))?')


def get_parser():
    parser = argparse.ArgumentParser(
        description= 'Consolidate variants including filtering out reference ' +
        'variants and calling overlapping variants.')
    parser.add_argument(
        'variants',
        help='Megalodon called variant file. Must contain GL sample field.')
    parser.add_argument(
        '--output-filename', default='megalodon.consolidated_variants.vcf',
        help='Output filename. Default: %(default)s')
    parser.add_argument(
        '--max-likelihood-ratio', type=float, default=1,
        help='Maximum likelihood ratio ([ref prob] / [max alt prob]) to include
        variant in output. Allows output of uncertain reference calls. ' +
        'Default: 1; Include only sites called as alternative.')
    parser.add_argument(
        '--min-depth', type=int,
        help='Minimum depth to include a variant. Default: No depth filter')
    parser.add_argument(
        '--trim-variants', action='store_true',
        help='Trim extra padding sequence included by megalodon (e.g. around ' +
        'repeat-region indels). Default: Output as found in input variants.')

    ssv_grp = parser.add_argument_group('Strand-specific Variant Arguments')
    ssv_grp.add_argument(
        '--reverse-strand-variants',
        help='Variants file produced only from reads mapping to the reverse ' +
        'strand. If provided, this assumes that the main variants file ' +
        'contains variants only supported by reads from the forward strand. ' +
        'This is used to identify systematic basecalling error variants. ' +
        'Errors made on both strands indicate potential putative variants ' +
        'and are thus excluded. Homopolymer variants occuring on both ' +
        'strands are included by default. Exclude these variants as well ' +
        'by setting --exclude-both-strand-homopolymers .')
    ssv_grp.add_argument(
        '--homopolymer-min-length', type=int, default=4,
        help='Minimum length to consider a variant as a homopolymer. ' +
        'Default: %(default)d')
    ssv_grp.add_argument(
        '--exclude-both-strand-homopolymers', action='store_true',
        help='By default homopolymer variants are included even if they ' +
        'occur on both strands. Set this flag to treat homopolymer variants ' +
        'as other variants.')

    return parser


class Variant(object):
    def __init__(self, line, do_trim_var, max_lr, min_depth):
        self.do_trim_var = do_trim_var
        self.max_lr = max_lr
        self.min_depth = min_depth
        (self.chrm, self.pos, self.id, self.ref, self.alts, self.qual,
         self.filt, self.raw_info, self.raw_fmt,
         self.raw_sample) = line.strip().split('\t')
        self.start = int(self.pos)
        self.end = self.start + len(self.ref)
        self.info = OrderedDict((
            m.groups() for m in INFO_PAT.finditer(self.raw_info)))
        self.sample = OrderedDict(zip(self.raw_fmt.split(':'),
                                      self.raw_sample.split(':')))
        self.alts = self.alts.split(',')
        self.gt_probs = 10 ** np.array([
            float(x) for x in self.sample['GL'].split(',')])
        return

    @property
    def ref_prob(self):
        return self.gt_probs[0]

    @property
    def depth(self):
        return int(self.sample['DP'])

    @property
    def hp_len(self):
        return np.max([np.max(np.diff(np.concatenate([
            [-1,],
            np.where(np.array([
                b1 != b2 for b1, b2 in
                zip(allele_seq[:-1], allele_seq[1:])]))[0],
            [len(allele_seq) - 1,]])))
                       for allele_seq in [self.ref,] + self.alts])

    @property
    def do_output(self):
        if self.min_depth is not None and self.depth < self.min_depth:
            return False
        try:
            lr = self.gt_probs[0] / max(self.gt_probs[1:])
        except ZeroDivisionError:
            return False
        return lr <= self.max_lr

    def overlaps(self, variant2):
        return not (self.start >= variant2.end or variant2.start >= self.end)

    @staticmethod
    def generate_gts(ploidy, num_alleles):
        if ploidy == 1:
            return [[a,] for a in range(num_alleles)]
        assert ploidy == 2, ('Cannot process variants with ploidy ' +
                             'greater than 2.')
        gts = []
        for a2 in range(num_alleles):
            for a1 in range(a2 + 1):
                gts.append((a1, a2))
        return gts

    def set_strand(self, strand):
        self.info[mh.STRAND_FIELD_NAME] = strand
        return

    def select_var(self):
        """ Select most likely variant allele
        """
        gt = GT_PAT.search(self.sample['GT'])
        if gt is None: raise mh.MegaError('Invalid genotype: {}'.format(
                self.sample['GT']))
        gt = gt.groupdict()
        ploidy = 1 if gt['a2'] is None else 2
        prev_gts = self.generate_gts(ploidy, len(self.alts) + 1)

        # pick most likely alt allele
        alt_idx = np.argmax(self.gt_probs[1:])
        prev_alt_gt = prev_gts[alt_idx + 1]
        prev_alt_alleles = sorted([a for a in set(prev_alt_gt) if a != 0])
        self.alts = [self.alts[a - 1] for a in prev_alt_alleles]
        allele_conv = dict(
            (a, new_a + 1) for new_a, a in enumerate(prev_alt_alleles))
        allele_conv[0] = 0
        alt_gt = (
            '{}'.format(allele_conv[prev_alt_alleles[0]]) if ploidy == 1 else
            '{}{}{}'.format(allele_conv[prev_alt_gt[0]], gt['sep'],
                            allele_conv[prev_alt_gt[1]]))
        ref_gt = '0' if ploidy == 1 else '0{}0'.format(gt['sep'])

        new_prob_sum = self.gt_probs[0] + self.gt_probs[alt_idx + 1]
        self.gt_probs = np.array([self.gt_probs[0] / new_prob_sum,
                                  self.gt_probs[alt_idx + 1] / new_prob_sum])
        with np.errstate(divide='ignore'):
            gl = np.maximum(mh.MIN_GL_VALUE, np.log10(self.gt_probs))
        raw_pl = -10 * gl
        pl = np.abs(np.minimum(raw_pl - raw_pl.min(), mh.MAX_PL_VALUE))
        s_pl = np.sort(pl)

        # set output values with single alt allele
        self.sample['GT'] = ref_gt if self.ref_prob > 0.5 else alt_gt
        self.sample['GL'] = ','.join(
            '{:.2f}' for _ in range(gl.shape[0])).format(*gl)
        self.sample['PL'] = ','.join(
            '{:.0f}' for _ in range(gl.shape[0])).format(*np.around(pl))
        self.sample['GQ'] = '{:.0f}'.format(np.around(s_pl[1]))
        self.qual = str(int(np.around(np.minimum(raw_pl[0], mh.MAX_PL_VALUE))))
        return

    def trim_variant(self):
        # trim padding positions from end of variant
        while (min(len(self.ref), min(len(alt) for alt in self.alts)) > 1 and
               all(self.ref[-1] == alt[-1] for alt in self.alts)):
            self.ref = self.ref[:-1]
            for i in range(len(self.alts)):
                self.alts[i] = self.alts[i][:-1]
        # trim padding positions from beginning of variant
        while (min(len(self.ref), min(len(alt) for alt in self.alts)) > 1 and
               all(self.ref[0] == alt[0] for alt in self.alts)):
            self.ref = self.ref[1:]
            for i in range(len(self.alts)):
                self.alts[i] = self.alts[i][1:]
            self.start += 1
        return

    def write_variant(self, out_fp):
        if self.do_trim_var: self.trim_variant()
        out_fp.write('\t'.join(map(str, (
            self.chrm, self.start, self.id, self.ref, ','.join(self.alts),
            self.qual, self.filt,
            ';'.join((self.info[k] if v is None else '{}={}'.format(k,v)
                      for k, v in self.info.items())), self.raw_fmt,
            ':'.join((v for v in self.sample.values()))))) + '\n')


def iter_valid_variants(vars_fn, do_trim_vars, max_lr, min_depth):
    header = ''
    with open(vars_fn) as fp:
        for line in fp:
            if line.startswith('#'):
                header += line
            else:
                prev_var = Variant(line, do_trim_vars, max_lr, min_depth)
                break
        # yeild header first
        yield header

        for line in tqdm(fp, desc=vars_fn, smoothing=0):
            variant = Variant(line, do_trim_vars, max_lr, min_depth)
            if variant.overlaps(prev_var):
                # for overlapping variants select the more likely alternative
                if variant.ref_prob < prev_var.ref_prob:
                    prev_var = variant
            else:
                if prev_var.do_output:
                    prev_var.select_var()
                    yield prev_var
                prev_var = variant

        if prev_var.do_output:
            prev_var.select_var()
            yield prev_var

    return


def main():
    args = get_parser().parse_args()
    out_fp = open(args.output_filename, 'w')
    vars_iter = iter_valid_variants(
        args.variants, args.trim_variants, args.max_likelihood_ratio,
        args.min_depth)
    rev_vars = next_rev_var = None
    if args.reverse_strand_variants:
        rev_vars = iter_valid_variants(
            args.reverse_strand_variants, args.trim_variants,
            args.max_likelihood_ratio, args.min_depth)
        _ = next(rev_vars)
        next_rev_var = next(rev_vars)

    header = next(vars_iter)
    if rev_vars is not None:
        header += ('##INFO=<ID={},Number=1,Type=String,Description' +
                   '="Variant Evidence Strand">\n').format(mh.STRAND_FIELD_NAME)
    out_fp.write(header)
    for var in vars_iter:
        # output non-overlapping variants from reverse strand
        while (rev_vars is not None and next_rev_var is not None and
               next_rev_var.end < var.start):
            next_rev_var.set_strand('-1')
            next_rev_var.write_variant(out_fp)
            next_rev_var = next(rev_vars)
        if next_rev_var is None or var.end < next_rev_var.start:
            # non-overlapping
            if rev_vars is not None: var.set_strand('1')
            var.write_variant(out_fp)
        else:
            # overlapping variants from strands
            if (not args.exclude_both_strand_homopolymers and
                var.hp_len >= args.homopolymer_min_length):
                # note don't set strand here as evidence is from both strands
                var.write_variant(out_fp)
            if rev_vars is not None:
                try:
                    next_rev_var = next(rev_vars)
                except StopIteration:
                    next_rev_var = None

    return


if __name__ == '__main__':
    main()
