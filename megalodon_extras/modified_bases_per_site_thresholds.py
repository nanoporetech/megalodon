import sys
from collections import defaultdict

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from megalodon import logging, mods, megalodon_helper as mh
from ._extras_parsers import get_parser_modified_bases_per_site_thresholds


BED_TMPLT = '{chrom}\t{pos}\t{end}\t.\t{score}\t{strand}\n'
COV_PCTLS = [1.0, 5.0, 10.0, 15.0, 20.0, 25.0]
MAX_MOD_SCORE = 100

LOGGER = logging.get_logger()


def write_per_site_threshs(
        cov, meth_cov, blk_lst_fn, thresh_fn, mod_db_fn, strand_offset,
        gt_cov_min, np_cov_min, target_mod_base):
    def get_gt_cov(chrm, strand, pos):
        cov_pos = pos if strand_offset is None else (
            pos - strand_offset if strand == -1 else pos)
        cov_strand = strand if strand_offset is None else None
        pos_cov = 0
        try:
            pos_cov = cov[(chrm, cov_strand)][cov_pos]
            if pos_cov < gt_cov_min:
                return pos_cov, 0
            return pos_cov, meth_cov[(chrm, cov_strand)][cov_pos]
        except KeyError:
            return 0, 0

    def get_pos_thresh(pos_cov, pos_meth_cov, pos_mod_data):
        if pos_meth_cov == 0:
            return '{:.4f}'.format(-MAX_MOD_SCORE)
        elif pos_meth_cov == pos_cov:
            return '{:.4f}'.format(MAX_MOD_SCORE)

        # determine fractional threshold
        # collate stats per-read
        pos_stats = defaultdict(dict)
        for read_dbid, mod_dbid, lp in pos_mod_data:
            pos_stats[read_dbid][mod_dbid] = lp
        pos_llrs = []
        with np.errstate(divide='ignore'):
            for read_pos_lps in pos_stats.values():
                mt_lps = list(read_pos_lps.values())
                can_lp = np.log1p(-np.exp(mt_lps).sum())
                for mod_dbid, mod_lp in read_pos_lps.items():
                    if mod_dbid == target_mod_dbid:
                        pos_llrs.append(can_lp - mod_lp)
        gt_meth_pct = 100.0 * pos_meth_cov / pos_cov
        return '{:.4f}'.format(np.percentile(pos_llrs, gt_meth_pct))

    blk_lst_fp = open(blk_lst_fn, 'w')
    thresh_fp = open(thresh_fn, 'w')
    mods_db = mods.ModsDb(mod_db_fn)
    target_mod_dbid = mods_db.get_mod_base_dbid(target_mod_base)
    bar = tqdm(total=mods_db.get_num_uniq_stats(), smoothing=0,
               dynamic_ncols=True, unit='stats')
    for (chrm, strand, pos), pos_mod_data in mods_db.iter_pos_scores(
            convert_pos=True):
        bar.update(len(pos_mod_data))
        str_strand = mh.int_strand_to_str(strand)

        # extract ground truth coverage
        pos_cov, pos_meth_cov = get_gt_cov(chrm, strand, pos)
        if pos_cov < gt_cov_min:
            blk_lst_fp.write(
                BED_TMPLT.format(
                    chrom=chrm, pos=pos, end=pos + 1, strand=str_strand,
                    score='GT_COV:{}'.format(pos_cov)))
            continue

        # if nanopore coverage is not deep enough write to blacklist
        target_mod_cov = len(set(
            read_id for read_id, mod_dbid, _ in pos_mod_data
            if mod_dbid == target_mod_dbid))
        if target_mod_cov < np_cov_min:
            # don't blacklist sites covered by other called mods
            if target_mod_cov > 0:
                blk_lst_fp.write(
                    BED_TMPLT.format(
                        chrom=chrm, pos=pos, end=pos + 1, strand=str_strand,
                        score='NP_COV:{}'.format(target_mod_cov)))
            continue

        pos_thresh = get_pos_thresh(pos_cov, pos_meth_cov, pos_mod_data)
        thresh_fp.write(BED_TMPLT.format(
            chrom=chrm, pos=pos, end=pos + 1, strand=str_strand,
            score=pos_thresh))

    bar.close()
    blk_lst_fp.close()
    thresh_fp.close()


def summarize_ground_truth(cov, gt_cov_min, gt_pdf):
    LOGGER.info('Computing ground truth coverage statistics')
    all_gt_cov = np.array([c for cs_cov in cov.values()
                           for c in cs_cov.values()], dtype=int)
    gt_pctls = np.percentile(all_gt_cov, COV_PCTLS)
    thresh_sites_pct = 100 * sum(
        all_gt_cov >= gt_cov_min) / all_gt_cov.shape[0]
    LOGGER.info((
        'Ground truth coverage threshold includes {:.2f}% of ground ' +
        'truth sites (with at least one read).').format(thresh_sites_pct))
    LOGGER.info(
        'Ground truth coverage percentiles (pctl:cov): ' +
        ''.join('{:.0f}:{:.0f}  '.format(*x)
                for x in zip(COV_PCTLS, gt_pctls)))
    if gt_pdf is not None:
        bin_cnts = np.bincount(all_gt_cov)
        plot_thresh = np.percentile(all_gt_cov, [0.1, 99.0])
        bin_cnts = bin_cnts[int(plot_thresh[0]):int(plot_thresh[-1])]
        pdf_fp = PdfPages(gt_pdf)
        plt.bar(np.arange(int(plot_thresh[0]), int(plot_thresh[-1])), bin_cnts)
        pdf_fp.savefig(bbox_inches='tight')
        pdf_fp.close()


def check_matching_attrs(
        ground_truth_beds, strand_offset, mod_db_fn, target_mod_base,
        limit=10000):
    mods_db = mods.ModsDb(mod_db_fn)
    db_strands = (1, -1) if strand_offset is None else (None, )
    db_chrms = set((chrm, strand)
                   for _, chrm, _ in mods_db.iter_chrms()
                   for strand in db_strands)
    cov, meth_cov = mh.parse_bed_methyls(
        ground_truth_beds, strand_offset, show_prog_bar=False, limit=limit)
    if len(db_chrms.intersection(cov.keys())) == 0:
        LOGGER.error(
            ('Using first {} reads from {}, found zero overlapping ' +
             'contig/chromosome names with the mod database.').format(
                 limit, ground_truth_beds[0]))
        LOGGER.info('Database contigs/chromosomes: {}'.format(', '.join(
            map(str, db_chrms))))
        LOGGER.info('BED methyl contigs/chromosomes: {}'.format(', '.join(
            map(str, list(cov.keys())))))
        raise mh.MegaError('No overlapping contigs found.')
    db_mods = set(mod_base for mod_base, _ in mods_db.get_mod_long_names())
    if target_mod_base not in db_mods:
        raise mh.MegaError('Target modified base not found in mods database.')
    mods_db.check_data_covering_index_exists()
    mods_db.close()


def _main(args):
    logging.init_logger()
    LOGGER.info('Parsing ground truth bed methyl files')
    mod_db_fn = mh.get_megalodon_fn(args.megalodon_results_dir, mh.PR_MOD_NAME)
    check_matching_attrs(
        args.ground_truth_beds, args.strand_offset, mod_db_fn, args.mod_base)
    cov, meth_cov = mh.parse_bed_methyls(
        args.ground_truth_beds, args.strand_offset)
    if not args.skip_ground_truth_summary:
        summarize_ground_truth(
            cov, args.ground_truth_cov_min, args.ground_truth_coverage_pdf)
    if args.ground_truth_cov_only:
        sys.exit(1)

    LOGGER.info('Processing nanopore modified base results')
    write_per_site_threshs(
        cov, meth_cov, args.out_blacklist_sites,
        args.out_per_site_mod_thresholds, mod_db_fn, args.strand_offset,
        args.ground_truth_cov_min, args.nanopore_cov_min, args.mod_base)


if __name__ == '__main__':
    _main(get_parser_modified_bases_per_site_thresholds().parse_args())
