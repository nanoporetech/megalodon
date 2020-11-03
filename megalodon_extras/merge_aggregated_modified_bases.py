import heapq

import numpy as np
from tqdm import tqdm

from megalodon import logging, megalodon_helper as mh, mods
from ._extras_parsers import get_parser_merge_aggregated_modified_bases


def _main(args):
    logging.init_logger()
    out_fp = open(args.output_bed_methyl_file, 'w')

    bar = tqdm(desc='Records Written', smoothing=0)
    if args.sorted_inputs:
        curr_chrm = curr_pos = curr_strand = None
        mod_cov = cov = 0
        for chrm, pos, strand, mod_covi, covi in heapq.merge(*[
                mh.iter_bed_methyl_recs(in_fn)
                for in_fn in args.bed_methyl_files]):
            if curr_strand != strand or curr_pos != pos or curr_chrm != chrm:
                if curr_chrm is not None and cov > 0:
                    out_fp.write(mods.BEDMETHYL_TMPLT.format(
                        chrom=curr_chrm, pos=curr_pos, end=curr_pos + 1,
                        strand=mh.int_strand_to_str(curr_strand), cov=cov,
                        score=min(cov, 1000),
                        perc=np.around(mod_cov * 100 / cov, 1)) + '\n')
                    bar.update()
                curr_chrm, curr_pos, curr_strand = chrm, pos, strand
                mod_cov = cov = 0
            mod_cov += mod_covi
            cov += covi
        if curr_chrm is not None and cov > 0:
            out_fp.write(mods.BEDMETHYL_TMPLT.format(
                chrom=curr_chrm, pos=curr_pos, end=curr_pos + 1,
                strand=mh.int_strand_to_str(curr_strand), cov=cov,
                score=min(int(cov), 1000),
                perc=np.around(mod_cov / cov * 100, 1)) + '\n')
            bar.update()
    else:
        cov, mod_cov = mh.parse_bed_methyls(args.bed_methyl_files)
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
    bar.close()
    out_fp.close()


if __name__ == '__main__':
    _main(get_parser_merge_aggregated_modified_bases().parse_args())
