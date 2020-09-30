import os
import sys
from collections import defaultdict, namedtuple

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.metrics import (
    roc_curve, auc, precision_recall_curve, average_precision_score)

from megalodon import logging, mapping, megalodon_helper as mh, mods
from ._extras_parsers import get_parser_validate_results


LOGGER = logging.get_logger()

PLOT_MIN_BC_ACC = 80
MOD_BANDWIDTH = 0.2
BC_BANDWIDTH = 0.2
LEN_BANDWIDTH = 10
GRIDSIZE = 1000

BC_LEGEND_LABEL = 'Sample'
DEFAULT_VS_LABEL = 'All Sites'

ACC_METRICS_HEADER = (
    '{: <17}{: <15}{: <15}{: <11}{: <15}{: <15}{: <15}{}\n'.format(
        'Median_Accuracy', 'Mean_Accuracy', 'Mode_Accuracy', 'Num_Reads',
        'Longest_Aligned_Len', 'Median_Aligned_Len', 'Mean_Aligned_Len',
        'Sample_Label'))
ACC_METRICS_TMPLT = (
    '{: <17.4}{: <15.4f}{: <15.1f}{: <11d}{: <15.1f}{: <15.1f}{: <15.1f}{}\n')

MOD_MISSING_MSG = ('{0} not found in "{1}", "{2}", "{3}" sites. ' +
                   'Skipping validation for "{0}" + "{1}" + "{2}".')

VAL_MOD_DATA = namedtuple('VAL_MOD_DATA', (
    'acc', 'parsim_acc', 'aligned_lens', 'mod_data', 'ctrl_data', 'label'))

MOD_VAL_METRICS_HEADER = (
    '{: <12}{: <19}{: <20}{: <9}{: <20}{: <19}{: <10}{}  {}\n'.format(
        'Optimal_F1', 'Optimal_Threshold', 'Mean_Avg_Precision', 'ROC_AUC',
        'Num_Modified_Stats', 'Num_Control_Stats', 'Mod_Base', 'Sample_Label',
        'Valid_Sites_Label'))
MOD_VAL_METRICS_TMPLT = (
    '{: <12.6f}{: <19.4f}{: <20.6f}{: <9.6f}{: <20d}{: <19d}{: <10}{}  {}\n')


def plot_nmap_reads(pdf_fp, samp_labs, nmapped_reads):
    LOGGER.info('Plotting number of mapped reads')


def plot_pr(pdf_fp, pr_data):
    for mod_base, mod_pr_data in pr_data.items():
        LOGGER.info('Plotting {} precision-recall curves'.format(mod_base))
        plt.figure(figsize=(8, 7))
        for lab, prec, recall in mod_pr_data:
            plt.step(recall, prec, label=lab, where='post')
        plt.ylim([-0.05, 1.05])
        plt.xlim([-0.05, 1.05])
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.title(('Modified Base "{}" Precision-Recall Curves').format(
            mod_base))
        plt.legend()
        pdf_fp.savefig(bbox_inches='tight')
        plt.close()


def plot_roc(pdf_fp, roc_data):
    for mod_base, mod_roc_data in roc_data.items():
        LOGGER.info('Plotting {} ROC curves'.format(mod_base))
        plt.figure(figsize=(8, 7))
        for lab, fpr, tpr in mod_roc_data:
            plt.step(fpr, tpr, label=lab)
        plt.ylim([-0.05, 1.05])
        plt.xlim([-0.05, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(('Modified Base "{}" ROC Curves').format(mod_base))
        plt.legend()
        pdf_fp.savefig(bbox_inches='tight')
        plt.close()


def plot_kde(pdf_fp, kde_data):
    for samp_lab, mod_stats, ctrl_stats in kde_data:
        LOGGER.info(
            'Plotting {} modified base statistics densities'.format(
                samp_lab))
        plt.figure(figsize=(8, 5))
        sns.kdeplot(mod_stats, shade=True, bw=MOD_BANDWIDTH,
                    gridsize=GRIDSIZE, label='Yes')
        sns.kdeplot(ctrl_stats, shade=True, bw=MOD_BANDWIDTH,
                    gridsize=GRIDSIZE, label='No')
        plt.legend(prop={'size': 16}, title='Is Modified')
        plt.xlabel('Log Likelihood Ratio\nMore Likely Modified <--> ' +
                   'More Likely Canonical')
        plt.ylabel('Density')
        plt.title(samp_lab)
        pdf_fp.savefig(bbox_inches='tight')
        plt.close()


def compute_mod_sites_stats(
        mod_stats, ctrl_stats, balance_classes, mod_base, samp_lab, vs_lab,
        out_fp):
    if balance_classes:
        # randomly downsample sample with more observations
        if mod_stats.shape[0] > ctrl_stats.shape[0]:
            mod_stats = np.random.choice(
                mod_stats, ctrl_stats.shape[0], replace=False)
        elif mod_stats.shape[0] < ctrl_stats.shape[0]:
            ctrl_stats = np.random.choice(
                ctrl_stats, mod_stats.shape[0], replace=False)

    is_can = np.repeat([0, 1], [mod_stats.shape[0], ctrl_stats.shape[0]])
    all_stats = np.concatenate([mod_stats, ctrl_stats])
    if any(np.isnan(all_stats)):
        LOGGER.warning((
            'Encountered {} NaN modified base scores.').format(
                sum(np.isnan(all_stats))))
        all_stats = all_stats[~np.isnan(all_stats)]
        if all_stats.shape[0] == 0:
            raise mh.MegaError('All modified base scores are NaN')
    inf_idx = np.isinf(all_stats)
    if any(inf_idx):
        all_stats[inf_idx] = np.max(all_stats[~inf_idx])
    neginf_idx = np.isinf(all_stats)
    if any(neginf_idx):
        all_stats[neginf_idx] = np.min(all_stats[~neginf_idx])
    LOGGER.info(
        'Computing PR/ROC for {} from {} at {}'.format(
            mod_base, samp_lab, vs_lab))
    # compute roc and presicion recall
    precision, recall, thresh = precision_recall_curve(is_can, all_stats)
    prec_recall_sum = precision + recall
    valid_idx = np.where(prec_recall_sum > 0)
    all_f1 = (2 * precision[valid_idx] * recall[valid_idx] /
              prec_recall_sum[valid_idx])
    optim_f1_idx = np.argmax(all_f1)
    optim_f1 = all_f1[optim_f1_idx]
    optim_thresh = thresh[optim_f1_idx]
    avg_prcn = average_precision_score(is_can, all_stats)

    fpr, tpr, _ = roc_curve(is_can, all_stats)
    roc_auc = auc(fpr, tpr)

    out_fp.write(
        MOD_VAL_METRICS_TMPLT.format(
            optim_f1, optim_thresh, avg_prcn, roc_auc, mod_stats.shape[0],
            ctrl_stats.shape[0], mod_base, samp_lab, vs_lab))
    pr_data = ('{} at {} mAP={:0.2f}'.format(
        samp_lab, vs_lab, avg_prcn), precision, recall)
    roc_data = ('{} at {} AUC={:0.2f}'.format(
        samp_lab, vs_lab, roc_auc), fpr, tpr)
    kde_data = ('{} from {} at {}'.format(
        mod_base, samp_lab, vs_lab), mod_stats, ctrl_stats)

    return pr_data, roc_data, kde_data


def report_mod_metrics(
        mod_samps_data, ctrl_samps_data, balance_classes, vs_labs,
        out_fp, pdf_fp):
    LOGGER.info('Computing modified base metrics')
    if vs_labs is None:
        vs_labs = [DEFAULT_VS_LABEL, ]
    all_mods_data = [msd.mod_data for msd in mod_samps_data]
    samp_labs = [msd.label for msd in mod_samps_data]
    if ctrl_samps_data is None:
        all_ctrl_data = [msd.ctrl_data for msd in mod_samps_data]
    else:
        # handle case where single control sample is provided for all mod
        # samples
        if len(ctrl_samps_data) == 1:
            all_ctrl_data = [
                ctrl_samps_data[0].mod_data for _ in mod_samps_data]
        else:
            all_ctrl_data = [
                ctrl_samp_data.mod_data for ctrl_samp_data in ctrl_samps_data]

    # extract all modified bases from all samples and all valid sites
    all_mod_bases = set((
        mod_base for samp_data in all_mods_data + all_ctrl_data
        for vs_samp_data in samp_data for mod_base in vs_samp_data))
    out_fp.write(MOD_VAL_METRICS_HEADER)
    all_pr_data, all_roc_data = defaultdict(list), defaultdict(list)
    all_kde_data = []
    # loop over samples
    for mod_samp_data, ctrl_samp_data, samp_lab in zip(
            all_mods_data, all_ctrl_data, samp_labs):
        # loop over valid site sets
        for vs_samp_mod_data, vs_samp_ctrl_data, vs_lab in zip(
                mod_samp_data, ctrl_samp_data, vs_labs):
            # loop over modified bases
            for mod_base in all_mod_bases:
                # check that mod_base exists in both data set
                if mod_base not in vs_samp_mod_data or \
                   vs_samp_mod_data[mod_base].shape[0] == 0:
                    LOGGER.warning(MOD_MISSING_MSG.format(
                        mod_base, samp_lab, vs_lab, 'modified'))
                    continue
                if mod_base not in vs_samp_ctrl_data or \
                   vs_samp_ctrl_data[mod_base].shape[0] == 0:
                    LOGGER.warning(MOD_MISSING_MSG.format(
                        mod_base, samp_lab, vs_lab, 'control'))
                    continue
                try:
                    # compute modified base metrics
                    pr_data, roc_data, kde_data = compute_mod_sites_stats(
                        vs_samp_mod_data[mod_base],
                        vs_samp_ctrl_data[mod_base],
                        balance_classes, mod_base, samp_lab, vs_lab, out_fp)
                    all_pr_data[mod_base].append(pr_data)
                    all_roc_data[mod_base].append(roc_data)
                    all_kde_data.append(kde_data)
                except mh.MegaError as e:
                    LOGGER.warning(str(e))

    plot_pr(pdf_fp, all_pr_data)
    plot_roc(pdf_fp, all_roc_data)
    plot_kde(pdf_fp, all_kde_data)


def plot_acc(pdf_fp, samps_val_data):
    # check that there are accuracies to be plotted, else return
    if all(samp_val_data.acc is None for samp_val_data in samps_val_data):
        return

    LOGGER.info('Plotting mapping accuracy distribution(s)')
    plt.figure(figsize=(8, 5))
    for samp_val_data in samps_val_data:
        if samp_val_data.acc is not None:
            sns.kdeplot(samp_val_data.acc, shade=False, bw=BC_BANDWIDTH,
                        gridsize=GRIDSIZE, label=samp_val_data.label)
    plt.legend(title=BC_LEGEND_LABEL)
    plt.xlabel('Mapping Accuracy')
    plt.ylabel('Density')
    plt.title('Mapping Accuracy')
    plt.xlim(PLOT_MIN_BC_ACC, 100)
    pdf_fp.savefig(bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(8, 5))
    for samp_val_data in samps_val_data:
        if samp_val_data.parsim_acc is not None:
            sns.kdeplot(samp_val_data.parsim_acc, shade=False, bw=BC_BANDWIDTH,
                        gridsize=GRIDSIZE, label=samp_val_data.label)
    plt.legend(title=BC_LEGEND_LABEL)
    plt.xlabel('Mapping Accuracy')
    plt.ylabel('Density')
    plt.title('Mapping Accuracy (Parsimonious: match - ins / ref_len)')
    plt.xlim(PLOT_MIN_BC_ACC, 100)
    pdf_fp.savefig(bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(8, 5))
    for samp_val_data in samps_val_data:
        if samp_val_data.aligned_lens is not None:
            sns.kdeplot(samp_val_data.aligned_lens, shade=False,
                        bw=LEN_BANDWIDTH, gridsize=GRIDSIZE,
                        label=samp_val_data.label)
    plt.legend(title=BC_LEGEND_LABEL)
    plt.xlabel('Aligned Length (Log10 scale)')
    plt.ylabel('Density')
    plt.title('Aligned Length (alignment_length - num_insertions)')
    plt.xscale('log', basex=10)
    pdf_fp.savefig(bbox_inches='tight')
    plt.close()

    samp_labs, nmapped_reads = [], []
    for samp_val_data in samps_val_data:
        if samp_val_data.acc is not None:
            samp_labs.append(samp_val_data.label)
            nmapped_reads.append(len(samp_val_data.acc))
    plt.figure(figsize=(8, 5))
    with sns.axes_style("whitegrid"):
        sns.barplot(x=samp_labs, y=nmapped_reads, hue=samp_labs, dodge=False)
    plt.legend([], [], frameon=False)
    plt.xlabel('Samples')
    plt.ylabel('Number of Mapped Reads')
    plt.title('Number of Mapped Reads')
    pdf_fp.savefig(bbox_inches='tight')
    plt.close()


def report_acc_metrics(res_dir, out_fp, samp_lab):
    try:
        bc_data = mapping.parse_map_summary_file(mh.get_megalodon_fn(
            res_dir, mh.MAP_SUMM_NAME))
        bc_acc = np.array([r_data.pct_identity for r_data in bc_data])
        parsim_acc = np.array([
            100 * (r_data.num_match - r_data.num_ins) /
            (r_data.num_align - r_data.num_ins) for r_data in bc_data])
        aligned_lens = np.array([r_data.num_align - r_data.num_ins
                                 for r_data in bc_data])
        # crude mode by rounding to 1 decimal
        uniq_acc, acc_counts = np.unique(np.around(
            bc_acc, 1), return_counts=True)
        mode_bc_acc = uniq_acc[np.argmax(acc_counts)]
        out_fp.write(ACC_METRICS_TMPLT.format(
            np.median(bc_acc), np.mean(bc_acc), mode_bc_acc, len(bc_data),
            np.max(aligned_lens), np.median(aligned_lens),
            np.mean(aligned_lens), samp_lab))
    except FileNotFoundError:
        bc_acc = parsim_acc = aligned_lens = None
        LOGGER.info('Mappings not found for {}'.format(res_dir))

    return bc_acc, parsim_acc, aligned_lens


def parse_mod_data(
        res_dir, out_fp, valid_sites, include_strand, samp_lab,
        ctrl_sites=None):
    mod_acc, parsim_acc, aligned_lens = report_acc_metrics(
        res_dir, out_fp, samp_lab)

    ctrl_data = None
    mods_db_fn = mh.get_megalodon_fn(res_dir, mh.PR_MOD_NAME)
    if os.path.exists(mods_db_fn):
        if ctrl_sites is not None:
            all_site_stats = mods.extract_stats_at_valid_sites(
                mods_db_fn, valid_sites + ctrl_sites,
                include_strand=include_strand)
            mods_data = all_site_stats[:len(valid_sites)]
            ctrl_data = all_site_stats[len(valid_sites):]
        elif valid_sites is not None:
            mods_data = mods.extract_stats_at_valid_sites(
                mods_db_fn, valid_sites, include_strand=include_strand)
        else:
            mods_data = [mods.extract_all_stats(mods_db_fn), ]
    else:
        mods_data = None

    return VAL_MOD_DATA(
        mod_acc, parsim_acc, aligned_lens, mods_data, ctrl_data, samp_lab)


def parse_valid_sites(valid_sites_fns, gt_data_fn, include_strand):
    if valid_sites_fns is None and gt_data_fn is None:
        return None, None, None

    # if ground truth file provided, parse first
    if gt_data_fn is not None:
        LOGGER.info('Reading ground truth file')
        gt_mod_pos, gt_ctrl_pos = mh.parse_ground_truth_file(
            gt_data_fn, include_strand=include_strand)
        if valid_sites_fns is None:
            # if ground truth provided, but not valid sites return parsed
            # ground truth sites.
            return [gt_mod_pos, ], None, [gt_ctrl_pos, ]

    # parse valid sites files and intersect with ground truth (if provided)
    LOGGER.info('Reading valid sites data')
    valid_sites, vs_labs = [], []
    ctrl_sites = None if gt_data_fn is None else []
    for vs_lab, valid_sites_fn in valid_sites_fns:
        try:
            vs_i_sites = mh.parse_beds([valid_sites_fn, ])
        except FileNotFoundError:
            LOGGER.warning('Could not find valid sites file: {}'.format(
                valid_sites_fn))
            continue

        vs_i_sites = set((
            (chrm, strand, pos) for (chrm, strand), cs_pos in vs_i_sites
            for pos in cs_pos))
        if gt_data_fn is None:
            valid_sites.append(vs_i_sites)
        else:
            ctrl_sites.append(vs_i_sites.intersection(gt_ctrl_pos))
            valid_sites.append(vs_i_sites.intersection(gt_mod_pos))
        vs_labs.append(vs_lab)

    if len(valid_sites) == 0:
        return None, None, None

    return valid_sites, vs_labs, ctrl_sites


def _main(args):
    logging.init_logger(quiet=args.quiet)
    pdf_fp = PdfPages(args.out_pdf)
    out_fp = (sys.stdout if args.out_filename is None else
              open(args.out_filename, 'w'))
    do_report_mod_metrics = (args.control_megalodon_results_dirs is not None or
                             args.ground_truth_data is not None)
    if args.control_megalodon_results_dirs is not None and \
       args.ground_truth_data is not None:
        LOGGER.warning('Cannot provide both control results and ground ' +
                       'truth file. Ignoring ground truth file.')
        args.ground_truth_data = None
    if args.control_megalodon_results_dirs is not None and \
       len(args.control_megalodon_results_dirs) > 1 and \
       len(args.control_megalodon_results_dirs) \
       != len(args.megalodon_results_dirs):
        LOGGER.error(
            'Must provide either one control results directory for all ' +
            'modified results directories or a control directory for each ' +
            'modified base results directory.')
        sys.exit(1)
    valid_sites, vs_labs, ctrl_sites = parse_valid_sites(
        args.valid_sites, args.ground_truth_data, args.strand_specific_sites)

    out_fp.write(ACC_METRICS_HEADER)
    LOGGER.info('Reading Megalodon results data')
    if args.results_labels is None:
        samp_labs = ['Sample {}'.format(samp_i + 1)
                     for samp_i in range(len(args.megalodon_results_dirs))]
    else:
        assert len(args.megalodon_results_dirs) == len(args.results_labels), (
            'Must be a label in --results-labels for each provided ' +
            'megalodon_results_dir')
        samp_labs = args.results_labels
    mod_samps_data = [
        parse_mod_data(
            mega_dir, out_fp, valid_sites, args.strand_specific_sites,
            samp_lab, ctrl_sites)
        for samp_lab, mega_dir in zip(samp_labs, args.megalodon_results_dirs)]
    ctrl_samps_data = None
    # if control is not specified via ground truth file, and control results
    # dir was provided parse control data
    if args.control_megalodon_results_dirs is not None:
        LOGGER.info('Reading Megalodon control data results')
        if len(args.control_megalodon_results_dirs) > 1:
            ctrl_samps_data = [
                parse_mod_data(
                    mega_dir, out_fp, valid_sites,
                    args.strand_specific_sites,
                    '{} Control'.format(samp_lab))
                for samp_lab, mega_dir in
                zip(samp_labs, args.control_megalodon_results_dirs)]
        else:
            # handle case with a single control for all mod dirs
            ctrl_samps_data = [parse_mod_data(
                args.control_megalodon_results_dirs[0], out_fp,
                valid_sites, args.strand_specific_sites, 'Control'), ]
        plot_acc(pdf_fp, mod_samps_data + ctrl_samps_data)
    else:
        plot_acc(pdf_fp, mod_samps_data)

    if do_report_mod_metrics:
        # enter newline between basecall accuracy and mod base results
        out_fp.write('\n')
        report_mod_metrics(
            mod_samps_data, ctrl_samps_data, not args.allow_unbalance_classes,
            vs_labs, out_fp, pdf_fp)

    pdf_fp.close()
    if out_fp is not sys.stdout:
        out_fp.close()


if __name__ == '__main__':
    _main(get_parser_validate_results().parse_args())
