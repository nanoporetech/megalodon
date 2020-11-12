import numpy as np
from tqdm import tqdm

from megalodon import logging, megalodon_helper as mh, mods
from ._extras_parsers import get_parser_merge_aggregated_modified_bases


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
    for chrm, pos, strand, mod_cov, cov in mh.iter_merged_bedmethyl([
            mh.iter_bed_methyl_recs(in_fn) for in_fn in in_fns]):
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
