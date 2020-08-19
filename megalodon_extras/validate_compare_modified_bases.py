import sys

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from megalodon import logging, megalodon_helper as mh
from ._extras_parsers import get_parser_validate_compare_modified_bases


LOGGER = logging.get_logger()

CMAP = plt.cm.inferno_r
COV_SAMP1_COLOR = '#7A7A7A'
COV_SAMP2_COLOR = '#0084A9'


def get_coverage_plot_data(
        samp1_cov, samp2_cov, trim_lower_pct=1, trim_upper_pctl=99):
    samp1_bin_cov = np.bincount(samp1_cov)
    samp2_bin_cov = np.bincount(samp2_cov)
    # find intervals included requested percentile of total coverage
    samp1_pctl = 100 * np.cumsum(samp1_bin_cov) / samp1_bin_cov.sum()
    samp2_pctl = 100 * np.cumsum(samp2_bin_cov) / samp2_bin_cov.sum()
    cov_start = min(np.searchsorted(samp1_pctl, trim_lower_pct, side='left'),
                    np.searchsorted(samp2_pctl, trim_lower_pct, side='left'))
    cov_end = max(np.searchsorted(samp1_pctl, trim_upper_pctl, side='right'),
                  np.searchsorted(samp2_pctl, trim_upper_pctl, side='right'))
    if samp1_bin_cov.shape[0] < cov_end:
        samp1_tmp = samp1_bin_cov
        samp1_bin_cov = np.zeros(cov_end)
        samp1_bin_cov[:samp1_tmp.shape[0]] = samp1_tmp
    if samp2_bin_cov.shape[0] < cov_end:
        samp2_tmp = samp2_bin_cov
        samp2_bin_cov = np.zeros(cov_end)
        samp2_bin_cov[:samp2_tmp.shape[0]] = samp2_tmp
    return (np.arange(cov_start, cov_end), samp1_bin_cov[cov_start:cov_end],
            samp2_bin_cov[cov_start:cov_end])


def plot_cov(
        samp1_cov, samp2_cov, samp1_valid_cov, samp2_valid_cov, samp_names,
        pdf_fp, width=0.35):
    x_ind, samp1_bin_valid_cov, samp2_bin_valid_cov = get_coverage_plot_data(
        samp1_valid_cov, samp2_valid_cov)
    plt.figure(figsize=(6, 5))
    plt.bar(x_ind, samp1_bin_valid_cov, width, color=COV_SAMP1_COLOR,
            label=samp_names[0])
    plt.bar(x_ind + width, samp2_bin_valid_cov, width, color=COV_SAMP2_COLOR,
            label=samp_names[1])
    plt.xlabel('Read Coverage')
    plt.ylabel('Number of Genomic Sites')
    plt.title('Read Coverage (Overlapping and Coverage Filtered Sites)')
    plt.legend()
    pdf_fp.savefig(bbox_inches='tight')
    plt.close()

    x_ind, samp1_bin_cov, samp2_bin_cov = get_coverage_plot_data(
        samp1_cov, samp2_cov)
    plt.figure(figsize=(6, 5))
    plt.bar(x_ind, samp1_bin_cov, width, color=COV_SAMP1_COLOR,
            label=samp_names[0])
    plt.bar(x_ind + width, samp2_bin_cov, width, color=COV_SAMP2_COLOR,
            label=samp_names[1])
    plt.xlabel('Read Coverage')
    plt.ylabel('Number of Genomic Sites')
    plt.title('Read Coverage (All Sites)')
    plt.legend()
    pdf_fp.savefig(bbox_inches='tight')
    plt.close()


def plot_hm(samp1_mod_pct, samp2_mod_pct, hm_num_bins, samp_names,
            out_fp, pdf_fp):
    corrcoef = np.corrcoef(samp1_mod_pct, samp2_mod_pct)[0, 1]
    out_fp.write('Correlation coefficient: {:.4f}\n'.format(corrcoef))
    out_fp.write('R^2: {:.4f}\n'.format(corrcoef**2))
    out_fp.write('RMSE (model: y=x): {:.4f}\n'.format(
        np.sqrt(np.mean(np.square(samp1_mod_pct - samp2_mod_pct)))))

    hm_ticks = [None, ] * hm_num_bins
    for mod_pct in (0, 25, 50, 75, 100):
        hm_ticks[int((hm_num_bins - 1) * mod_pct / 100)] = str(
            mod_pct)
    hist_data = np.histogram2d(
        samp2_mod_pct, samp1_mod_pct, bins=hm_num_bins,
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
    plt.xlabel('{} Percent Modified'.format(samp_names[0]))
    plt.ylabel('{} Percent Modified'.format(samp_names[1]))
    plt.title('Log10 Counts  N = {}  r = {:.4f}'.format(
        samp1_mod_pct.shape[0], corrcoef))
    pdf_fp.savefig(bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(6, 5))
    sns.heatmap(hist_data, cmap=CMAP, square=True, robust=True,
                xticklabels=hm_ticks, yticklabels=hm_ticks[::-1])
    plt.xlabel('{} Percent Methylated'.format(samp_names[0]))
    plt.ylabel('{} Percent Methylated'.format(samp_names[1]))
    plt.title('Raw Counts  N = {}  r = {:.4f}'.format(
        samp1_mod_pct.shape[0], corrcoef))
    pdf_fp.savefig(bbox_inches='tight')
    plt.close()


def compute_filt_mod_pct(
        samp1_cov, samp1_mod_cov, samp2_cov, samp2_mod_cov, cov_thresh,
        samp_names, out_fp):
    common_ctgs = set(samp1_cov).intersection(samp2_cov)
    if len(common_ctgs) == 0:
        LOGGER.error((
            'No common contigs found between provided samples.\n\t"{}" ' +
            'first 5 contigs: {}\n\t"{}" first 5 contigs: {}').format(
                samp_names[0], ','.join(sorted(list(samp1_cov.keys()))[:5]),
                samp_names[1], ','.join(sorted(list(samp2_cov.keys()))[:5])))
        sys.exit(1)

    # compute methylation percentages
    samp1_mod_pct, samp2_mod_pct = [], []
    samp1_valid_cov, samp2_valid_cov = [], []
    for ctg in common_ctgs:
        for pos in set(samp1_cov[ctg]).intersection(samp2_cov[ctg]):
            samp1_pos_cov = samp1_cov[ctg][pos]
            samp2_pos_cov = samp2_cov[ctg][pos]
            if min(samp1_pos_cov, samp2_pos_cov) < cov_thresh:
                continue
            samp1_mod_pct.append(
                100 * samp1_mod_cov[ctg][pos] / samp1_pos_cov)
            samp2_mod_pct.append(
                100 * samp2_mod_cov[ctg][pos] / samp2_pos_cov)
            samp1_valid_cov.append(samp1_pos_cov)
            samp2_valid_cov.append(samp2_pos_cov)

    if len(samp1_valid_cov) == 0:
        LOGGER.error(
            'No overlapping positions identified between provided samples.')
        sys.exit(1)

    out_fp.write(
        '{} valid coverage median: {:.2f}   mean: {:.2f}  sd: {:.2f}\n'.format(
            samp_names[0], np.median(samp1_valid_cov),
            np.mean(samp1_valid_cov), np.std(samp1_valid_cov)))
    out_fp.write(
        '{} valid coverage median: {:.2f}   mean: {:.2f}  sd: {:.2f}\n'.format(
            samp_names[1], np.median(samp2_valid_cov),
            np.mean(samp2_valid_cov), np.std(samp2_valid_cov)))
    out_fp.write('{} valid corresponding sites detected\n'.format(
        len(samp1_mod_pct)))

    return (np.array(samp1_mod_pct), np.array(samp2_mod_pct),
            np.array(samp1_valid_cov), np.array(samp2_valid_cov))


def parse_inputs(
        samp1_bm_fns, samp2_bm_fns, strand_offset, samp_names, valid_pos_fn,
        out_fp):
    # parse valid positions
    valid_pos = None
    if valid_pos_fn is not None:
        valid_pos = mh.parse_beds(
            valid_pos_fn, ignore_strand=strand_offset is not None)

    # parse bed methyl files
    samp1_cov, samp1_mod_cov = mh.parse_bed_methyls(
        samp1_bm_fns, strand_offset=strand_offset, valid_pos=valid_pos)
    samp1_all_cov = np.array([cov for ctg_cov in samp1_cov.values()
                              for cov in ctg_cov.values()])
    samp2_cov, samp2_mod_cov = mh.parse_bed_methyls(
        samp2_bm_fns, strand_offset=strand_offset, valid_pos=valid_pos)
    samp2_all_cov = np.array([cov for ctg_cov in samp2_cov.values()
                              for cov in ctg_cov.values()])
    out_fp.write(
        '{} coverage median: {:.2f}   mean: {:.2f}  sd: {:.2f}\n'.format(
            samp_names[0], np.median(samp1_all_cov),
            np.mean(samp1_all_cov), np.std(samp1_all_cov)))
    out_fp.write(
        '{} coverage median: {:.2f}   mean: {:.2f}  sd: {:.2f}\n'.format(
            samp_names[1], np.median(samp2_all_cov),
            np.mean(samp2_all_cov), np.std(samp2_all_cov)))

    return (samp1_cov, samp1_mod_cov, samp1_all_cov,
            samp2_cov, samp2_mod_cov, samp2_all_cov)


def _main(args):
    logging.init_logger()
    pdf_fp = PdfPages(args.out_pdf)
    out_fp = (sys.stdout if args.out_filename is None else
              open(args.out_filename, 'w'))

    (samp1_cov, samp1_mod_cov, samp1_all_cov,
     samp2_cov, samp2_mod_cov, samp2_all_cov) = parse_inputs(
         args.sample1_bed_methyl_files, args.sample2_bed_methyl_files,
         args.strand_offset, args.sample_names, args.valid_positions, out_fp)
    (samp1_mod_pct, samp2_mod_pct,
     samp1_valid_cov, samp2_valid_cov) = compute_filt_mod_pct(
         samp1_cov, samp1_mod_cov, samp2_cov, samp2_mod_cov,
         args.coverage_threshold, args.sample_names, out_fp)
    plot_hm(
        samp1_mod_pct, samp2_mod_pct, args.heatmap_num_bins, args.sample_names,
        out_fp, pdf_fp)
    plot_cov(samp1_all_cov, samp2_all_cov, samp1_valid_cov, samp2_valid_cov,
             args.sample_names, pdf_fp)

    pdf_fp.close()
    if out_fp is not sys.stdout:
        out_fp.close()


if __name__ == '__main__':
    _main(get_parser_validate_compare_modified_bases().parse_args())
