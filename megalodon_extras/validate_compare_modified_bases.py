import sys

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from megalodon import megalodon_helper as mh
from ._extras_parsers import get_parser_validate_compare_modified_bases


CMAP = plt.cm.inferno_r


def _main(args):
    pdf_fp = PdfPages(args.out_pdf)
    out_fp = (sys.stdout if args.out_filename is None else
              open(args.out_filename, 'w'))

    # parse bed methyl files
    samp1_cov, samp1_mod_cov = mh.parse_bed_methyls(
        args.sample1_bed_methyl_files, strand_offset=args.strand_offset)
    samp1_all_cov = np.array([cov for ctg_cov in samp1_cov.values()
                              for cov in ctg_cov.values()])
    samp2_cov, samp2_mod_cov = mh.parse_bed_methyls(
        args.sample2_bed_methyl_files, strand_offset=args.strand_offset)
    samp2_all_cov = np.array([cov for ctg_cov in samp2_cov.values()
                              for cov in ctg_cov.values()])
    out_fp.write(
        '{} coverage median: {:.2f}   mean: {:.2f}  sd: {:.2f}\n'.format(
            args.sample_names[0], np.median(samp1_all_cov),
            np.mean(samp1_all_cov), np.std(samp1_all_cov)))
    out_fp.write(
        '{} coverage median: {:.2f}   mean: {:.2f}  sd: {:.2f}\n'.format(
            args.sample_names[1], np.median(samp2_all_cov),
            np.mean(samp2_all_cov), np.std(samp2_all_cov)))

    # parse valid positions
    valid_pos = None
    if args.valid_positions is not None:
        valid_pos = mh.parse_beds(
            args.valid_positions, ignore_strand=args.strand_offset is not None)

    # compute methylation percentages
    samp1_meth_pct, samp2_meth_pct = [], []
    samp1_valid_cov, samp2_valid_cov = [], []
    for ctg in set(samp1_cov).intersection(samp2_cov):
        if valid_pos is not None:
            if ctg not in valid_pos:
                continue
            valid_ctg_pos = valid_pos[ctg]
        for pos in set(samp1_cov[ctg]).intersection(samp2_cov[ctg]):
            if valid_pos is not None and pos not in valid_ctg_pos:
                continue
            samp1_pos_cov = samp1_cov[ctg][pos]
            samp2_pos_cov = samp2_cov[ctg][pos]
            if min(samp1_pos_cov, samp2_pos_cov) < args.coverage_threshold:
                continue
            samp1_meth_pct.append(
                100 * samp1_mod_cov[ctg][pos] / samp1_pos_cov)
            samp2_meth_pct.append(
                100 * samp2_mod_cov[ctg][pos] / samp2_pos_cov)
            samp1_valid_cov.append(samp1_pos_cov)
            samp2_valid_cov.append(samp2_pos_cov)
    out_fp.write(
        '{} valid coverage median: {:.2f}   mean: {:.2f}  sd: {:.2f}\n'.format(
            args.sample_names[0], np.median(samp1_valid_cov),
            np.mean(samp1_valid_cov), np.std(samp1_valid_cov)))
    out_fp.write(
        '{} valid coverage median: {:.2f}   mean: {:.2f}  sd: {:.2f}\n'.format(
            args.sample_names[1], np.median(samp2_valid_cov),
            np.mean(samp2_valid_cov), np.std(samp2_valid_cov)))
    samp1_meth_pct = np.array(samp1_meth_pct)
    samp2_meth_pct = np.array(samp2_meth_pct)
    out_fp.write('{} sites detected\n'.format(samp1_meth_pct.shape[0]))
    corrcoef = np.corrcoef(samp1_meth_pct, samp2_meth_pct)[0, 1]
    out_fp.write('Correlation coefficient: {:.4f}\n'.format(corrcoef))
    out_fp.write('R^2: {:.4f}\n'.format(corrcoef**2))
    out_fp.write('Squared Error (model: y=x): {:.4f}\n'.format(
        np.mean(np.square(samp1_meth_pct - samp2_meth_pct))))

    hm_ticks = [None, ] * args.heatmap_num_bins
    for mod_pct in (0, 25, 50, 75, 100):
        hm_ticks[int((args.heatmap_num_bins - 1) * mod_pct / 100)] = str(
            mod_pct)
    hist_data = np.histogram2d(
        samp2_meth_pct, samp1_meth_pct, bins=args.heatmap_num_bins,
        range=[[0, 100], [0, 100]])[0][::-1]
    with np.errstate(divide='ignore'):
        log_hist_data = np.log10(hist_data)

    plt.figure(figsize=(6, 5))
    hm = sns.heatmap(
        log_hist_data, vmin=max(log_hist_data.min(), 1),
        vmax=log_hist_data.max(), cmap=CMAP, square=True,
        xticklabels=hm_ticks, yticklabels=hm_ticks[::-1])
    hm_cbar = hm.collections[0].colorbar
    hm_cbar.set_ticks(hm_cbar.get_ticks())
    hm_cbar.set_ticklabels([
        '{:.0f}'.format(tl) for tl in 10**hm_cbar.get_ticks()])
    plt.xlabel('{} Percent Modified'.format(args.sample_names[0]))
    plt.ylabel('{} Percent Modified'.format(args.sample_names[1]))
    plt.title('Log10 Counts  N = {}  r = {:.4f}'.format(
        samp1_meth_pct.shape[0], corrcoef))
    pdf_fp.savefig(bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(6, 5))
    sns.heatmap(hist_data, cmap=CMAP, square=True, robust=True,
                xticklabels=hm_ticks, yticklabels=hm_ticks[::-1])
    plt.xlabel('{} Percent Methylated'.format(args.sample_names[0]))
    plt.ylabel('{} Percent Methylated'.format(args.sample_names[1]))
    plt.title('Raw Counts  N = {}  r = {:.4f}'.format(
        samp1_meth_pct.shape[0], corrcoef))
    pdf_fp.savefig(bbox_inches='tight')
    plt.close()

    pdf_fp.close()
    if out_fp is not sys.stdout:
        out_fp.close()


if __name__ == '__main__':
    _main(get_parser_validate_compare_modified_bases().parse_args())
