import heapq

import numpy as np
from tqdm import tqdm

from megalodon import logging, megalodon_helper as mh, mods
from ._extras_parsers import get_parser_merge_aggregated_modified_bases


STRANDS = ('+', '-', '.')


def write_unsorted_merge(in_fns, out_fp, bar):
    cov, mod_cov = mh.parse_bed_methyls(in_fns)
    for chrm in sorted(
            mh.RefName(chrm) for chrm in set(chrm for chrm, _ in cov)):
        # convert back to string after sorting
        chrm = str(chrm)
        s_poss = []
        if (chrm, 1) in cov:
            s_poss.extend([(pos, 1) for pos in cov[(chrm, 1)]])
        if (chrm, -1) in cov:
            s_poss.extend([(pos, -1) for pos in cov[(chrm, -1)]])
        for pos, strand in sorted(s_poss):
            pcov = cov[(chrm, strand)][pos]
            out_fp.write(mods.BEDMETHYL_TMPLT.format(
                chrom=chrm, pos=pos, end=pos + 1,
                strand=mh.int_strand_to_str(strand), cov=pcov,
                score=min(int(pcov), 1000),
                perc=np.around(
                    mod_cov[(chrm, strand)][pos] / pcov * 100, 1)) + '\n')
            bar.update()


def iter_bed_methyl_recs(bed_fn, batch_size=10000):
    chrms, poss, strands, cov, pct_mods = [], [], [], [], []
    with open(bed_fn) as bed_fp:
        for line in bed_fp:
            (chrm, pos, _, _, _, strand, _, _, _, num_reads,
             pct_mod) = line.split()
            chrms.append(chrm)
            poss.append(pos)
            strands.append(strand)
            cov.append(num_reads)
            pct_mods.append(pct_mod)
            if len(chrms) >= batch_size:
                cov = np.array(cov, dtype=int)
                yield from zip(
                    chrms, map(int, poss), strands,
                    np.around(np.array(pct_mods, dtype=float) * cov / 100.0),
                    cov)
                chrms, poss, strands, cov, pct_mods = [], [], [], [], []
    if len(chrms) >= batch_size:
        cov = np.array(cov, dtype=int)
        yield from zip(
            chrms, map(int, poss), strands,
            np.around(np.array(pct_mods, dtype=float) * cov / 100.0),
            cov)

        """
            num_reads = int(num_reads)
            meth_reads = int(np.around(float(pct_meth) * num_reads / 100.0))
            yield chrm, int(pos), strand, meth_reads, num_reads
        """


def sorted_merge(rec_iters):
    # hold next record to process for each iterator
    curr_recs = [next(rec_iter, None) for rec_iter in rec_iters]
    # determine current chromosome and position
    curr_chrm = str(sorted(
        mh.RefName(rec[0]) for rec in curr_recs if rec is not None)[0])
    curr_pos = min(rec[1] for rec in curr_recs if rec[0] == curr_chrm)
    while any(rec is not None for rec in curr_recs):
        # collect all counts at the current position
        curr_pos_cnts = dict((strand, [0, 0]) for strand in STRANDS)
        for iter_i, (curr_rec, rec_iter) in enumerate(
                zip(curr_recs, rec_iters)):
            if curr_rec is None:
                continue
            chrm, pos, strand, mod_cov, cov = curr_rec
            while chrm == curr_chrm and pos == curr_pos:
                curr_pos_cnts[strand][0] += mod_cov
                curr_pos_cnts[strand][1] += cov
                curr_recs[iter_i] = next(rec_iter, None)
                if curr_recs[iter_i] is None:
                    break
                chrm, pos, strand, mod_cov, cov = curr_recs[iter_i]
        # return counts for this position
        for strand in STRANDS:
            if curr_pos_cnts[strand][1] > 0:
                yield (curr_chrm, curr_pos, strand, curr_pos_cnts[strand][0],
                       curr_pos_cnts[strand][1])
        # get next chrm/position
        if all(rec[0] != curr_chrm for rec in curr_recs if rec is not None):
            # determine next chromosome to process
            next_chrms = set(rec[0] for rec in curr_recs if rec is not None)
            if len(next_chrms) == 0:
                # no records left to process
                break
            elif len(next_chrms) == 1:
                curr_chrm = next(iter(next_chrms))
            else:
                # if there are multiple next chromosomes find next in
                # version sorted order
                curr_chrm = str(sorted(
                    mh.RefName(chrm) for chrm in next_chrms)[0])
        curr_pos = min(rec[1] for rec in curr_recs
                       if rec is not None and rec[0] == curr_chrm)

def write_batch(out_fp, chrms, poss, strands, mod_covs, covs):
    covs = np.array(covs, dtype=int)
    out_fp.write('\n'.join(
        mods.BEDMETHYL_TMPLT.format(
            chrom=chrm, pos=pos, end=pos + 1, strand=strand,
            cov=cov, score=score, perc=pct_mod)
        for chrm, pos, strand, cov, score, pct_mod in
        zip(chrms, poss, strands, covs, np.minimum(covs, 1000),
            np.around(np.array(mod_covs, dtype=int) * 100 / covs, 1))) + '\n')


def write_sorted_merge(in_fns, out_fp, bar, batch_size=50000):
    chrms, poss, strands, mod_covs, covs = [], [], [], [], []
    for chrm, pos, strand, mod_cov, cov in sorted_merge([
            iter_bed_methyl_recs(in_fn) for in_fn in in_fns]):
        chrms.append(chrm)
        poss.append(pos)
        strands.append(strand)
        mod_covs.append(mod_cov)
        covs.append(cov)
        bar.update()
        if len(chrms) >= batch_size:
            write_batch(out_fp, chrms, poss, strands, mod_covs, covs)
            chrms, poss, strands, mod_covs, covs = [], [], [], [], []
    if len(chrms) >= 0:
        write_batch(out_fp, chrms, poss, strands, mod_covs, covs)


def _main(args):
    logging.init_logger()
    with open(args.output_bed_methyl_file, 'w') as out_fp, \
         tqdm(desc='Records Written', smoothing=0) as bar:
        if args.sorted_inputs:
            write_sorted_merge(args.bed_methyl_files, out_fp, bar)
        else:
            write_unsorted_merge(args.bed_methyl_files, out_fp, bar)


if __name__ == '__main__':
    _main(get_parser_merge_aggregated_modified_bases().parse_args())
