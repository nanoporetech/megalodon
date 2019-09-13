import sys
import argparse
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument(
    'variants',
    help='Consolidate variants including filtering out reference variants, ' +
    'calling overlapping variants and removing extra unambiguous padding ' +
    'included by megalodon.')

class Variant(object):
    def __init__(self, line):
        (self.chrm, self.pos, self.id, self.ref, self.alts, self.qual,
         self.filt, self.info, self.raw_fmt,
         self.raw_sample) = line.strip().split('\t')
        self.start = int(self.pos)
        self.end = self.start + len(self.ref)
        self.sample = dict(zip(self.raw_fmt.split(':'),
                               self.raw_sample.split(':')))
        self.score = None
        return

    def get_score(self):
        if self.score is None:
            self.score = int(self.sample['GQ'])
        return self.score

    def overlaps(self, variant2):
        return not (self.start >= variant2.end or variant2.start >= self.end)

    def is_alt(self):
        return self.sample['GT'][0] != '0' or self.sample['GT'][-1] != '0'

    def write_variant(self, out_fp):
        split_alts = self.alts.split(',')
        # trim padding positions from end of variant
        while (min(len(self.ref), min(len(alt) for alt in split_alts)) > 1 and
               all(self.ref[-1] == alt[-1] for alt in split_alts)):
            self.ref = self.ref[:-1]
            for i in range(len(split_alts)):
                split_alts[i] = split_alts[i][:-1]
        # trim padding positions from beginning of variant
        while (min(len(self.ref), min(len(alt) for alt in split_alts)) > 1 and
               all(self.ref[0] == alt[0] for alt in split_alts)):
            self.ref = self.ref[1:]
            for i in range(len(split_alts)):
                split_alts[i] = split_alts[i][1:]
            self.start += 1
        out_fp.write('\t'.join(map(str, (
            self.chrm, self.start, self.id, self.ref, ','.join(split_alts),
            self.qual, self.filt, self.info, self.raw_fmt, self.raw_sample))) +
                     '\n')

def main():
    args = parser.parse_args()
    out_fp = sys.stdout
    prev_var = None
    with open(args.variants) as fp:
        for line in tqdm(fp):
            if line.startswith('#'):
                out_fp.write(line)
                continue
            variant = Variant(line)
            if prev_var is None:
                prev_var = variant
                continue

            if variant.overlaps(prev_var):
                if variant.get_score() > prev_var.get_score():
                    prev_var = variant
            else:
                if prev_var.is_alt():
                    prev_var.write_variant(out_fp)
                prev_var = variant

        if prev_var.is_alt():
            prev_var.write_variant(out_fp)

    return

if __name__ == '__main__':
    main()
