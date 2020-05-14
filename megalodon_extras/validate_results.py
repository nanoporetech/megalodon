import sys
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.metrics import (
    roc_curve, auc, precision_recall_curve, average_precision_score)

from megalodon import mapping, megalodon_helper as mh, mods
from ._extras_parsers import get_parser_validate_results


VERBOSE = False

PLOT_MIN_BC_ACC = 80
MOD_BANDWIDTH = 0.2
MOD_GRIDSIZE = 1000
BC_BANDWIDTH = 0.2
BC_GRIDSIZE = 1000

BC_LEGEND_LABEL = 'Sample'
BC_SAMPLE_NAME = 'Sample'
BC_CONTROL_NAME = 'Control'

DEFAULT_LABEL = 'All Sites'

MOD_MISSING_MSG = ('***** WARNING ***** {0} not found in {1} {2} sites. ' +
                   'Skipping {0} validation.')


def compute_mod_sites_stats(
        vs_mods_data, vs_ctrl_data, balance_classes,
        mod_base, vs_label, out_fp, pdf_fp):
    if balance_classes:
        # randomly downsample sample with more observations
        if vs_mods_data.shape[0] > vs_ctrl_data.shape[0]:
            vs_mods_data = np.random.choice(
                vs_mods_data, vs_ctrl_data.shape[0], replace=False)
        elif vs_mods_data.shape[0] < vs_ctrl_data.shape[0]:
            vs_ctrl_data = np.random.choice(
                vs_ctrl_data, vs_mods_data.shape[0], replace=False)

    is_can = np.repeat([0, 1], [vs_mods_data.shape[0], vs_ctrl_data.shape[0]])
    all_data = np.concatenate([vs_mods_data, vs_ctrl_data])
    if any(np.isnan(all_data)):
        sys.stderr.write((
            '***** WARNING ***** Encountered {} NaN modified base ' +
            'scores.\n').format(sum(np.isnan(all_data))))
        all_data = all_data[~np.isnan(all_data)]
        if all_data.shape[0] == 0:
            raise mh.MegaError('All modified base scores are NaN')
    inf_idx = np.isinf(all_data)
    if any(inf_idx):
        all_data[inf_idx] = np.max(all_data[~inf_idx])
    neginf_idx = np.isinf(all_data)
    if any(neginf_idx):
        all_data[neginf_idx] = np.min(all_data[~neginf_idx])
    if VERBOSE:
        sys.stderr.write(
            'Computing PR/ROC for {} at {}\n'.format(mod_base, vs_label))
    # compute roc and presicion recall
    precision, recall, thresh = precision_recall_curve(is_can, all_data)
    prec_recall_sum = precision + recall
    valid_idx = np.where(prec_recall_sum > 0)
    all_f1 = (2 * precision[valid_idx] * recall[valid_idx] /
              prec_recall_sum[valid_idx])
    optim_f1_idx = np.argmax(all_f1)
    optim_f1 = all_f1[optim_f1_idx]
    optim_thresh = thresh[optim_f1_idx]
    avg_prcn = average_precision_score(is_can, all_data)

    fpr, tpr, _ = roc_curve(is_can, all_data)
    roc_auc = auc(fpr, tpr)

    out_fp.write((
        'Modified base metrics for {} at {}:\t{:.6f} (at {:.4f} )\t' +
        '{:.6f}\t{:.6f}\t{}\t{}\n').format(
            mod_base, vs_label, optim_f1, optim_thresh, avg_prcn, roc_auc,
            vs_mods_data.shape[0], vs_ctrl_data.shape[0]))

    if VERBOSE:
        sys.stderr.write('Plotting {} at {}\n'.format(mod_base, vs_label))
    plt.figure(figsize=(11, 7))
    sns.kdeplot(vs_mods_data, shade=True, bw=MOD_BANDWIDTH,
                gridsize=MOD_GRIDSIZE, label='Yes')
    sns.kdeplot(vs_ctrl_data, shade=True, bw=MOD_BANDWIDTH,
                gridsize=MOD_GRIDSIZE, label='No')
    plt.legend(prop={'size': 16}, title='Is Modified')
    plt.xlabel('Log Likelihood Ratio\nMore Modified <--> Less Modified')
    plt.ylabel('Density')
    plt.title('Modified Base: {}  Sites: {}'.format(mod_base, vs_label))
    pdf_fp.savefig(bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(8, 7))
    plt.step(recall, precision, where='post')
    plt.ylim([-0.05, 1.05])
    plt.xlim([-0.05, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title(('Modified Base: {}  Sites: {}  Precision-Recall ' +
               'curve: AP={:0.2f}').format(mod_base, vs_label, avg_prcn))
    pdf_fp.savefig(bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(8, 7))
    plt.plot(fpr, tpr)
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(('Modified Base: {}  Sites: {}  ROC curve: ' +
               'AUC={:0.2f}').format(mod_base, vs_label, roc_auc))
    pdf_fp.savefig(bbox_inches='tight')
    plt.close()


def report_mod_metrics(
        mods_data, ctrl_data, balance_classes, vs_labels, out_fp, pdf_fp):
    if VERBOSE:
        sys.stderr.write('Computing modified base metrics\n')
    if vs_labels is None:
        vs_labels = [DEFAULT_LABEL, ]
        assert len(mods_data) == 1 and len(ctrl_data) == 1
    all_mod_bases = set((mod_base for vs_mods_data in mods_data
                         for mod_base in vs_mods_data)).union((
                                 mod_base for vs_ctrl_data in ctrl_data
                                 for mod_base in vs_ctrl_data))
    out_fp.write('Modified Base Metrics: Optimal F1 :: Optimal Threshold :: ' +
                 'Average Precision :: ROC AUC :: Num. Modified Sites :: ' +
                 'Num. Canonical Sites\n')
    for vs_mods_data, vs_ctrl_data, vs_label in zip(
            mods_data, ctrl_data, vs_labels):
        for mod_base in all_mod_bases:
            # check that mod_base exists in both data set
            if mod_base not in vs_mods_data or \
               vs_mods_data[mod_base].shape[0] == 0:
                sys.stderr.write(MOD_MISSING_MSG.format(
                    mod_base, vs_label, 'modified'))
                continue
            if mod_base not in vs_ctrl_data or \
               vs_ctrl_data[mod_base].shape[0] == 0:
                sys.stderr.write(MOD_MISSING_MSG.format(
                    mod_base, vs_label, 'control'))
                continue
            try:
                # compute modified base metrics
                compute_mod_sites_stats(
                    vs_mods_data[mod_base], vs_ctrl_data[mod_base],
                    balance_classes, mod_base, vs_label, out_fp, pdf_fp)
            except mh.MegaError as e:
                sys.stderr.write(str(e) + '\n')


def plot_acc(pdf_fp, mod_acc, mod_parsim_acc, ctrl_acc, ctrl_parsim_acc):
    if VERBOSE:
        sys.stderr.write('Plotting mapping accuracy distribution(s)\n')
    plt.figure(figsize=(11, 7))
    sns.kdeplot(mod_acc, shade=True,
                bw=BC_BANDWIDTH, gridsize=BC_GRIDSIZE, label=BC_SAMPLE_NAME)
    if ctrl_acc is not None:
        sns.kdeplot(ctrl_acc, shade=True,
                    bw=BC_BANDWIDTH, gridsize=BC_GRIDSIZE,
                    label=BC_CONTROL_NAME)
    plt.legend(prop={'size': 16}, title=BC_LEGEND_LABEL)
    plt.xlabel('Mapping Accuracy')
    plt.ylabel('Density')
    plt.title('Mapping Accuracy')
    plt.xlim(PLOT_MIN_BC_ACC, 100)
    pdf_fp.savefig(bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(11, 7))
    sns.kdeplot(mod_parsim_acc, shade=True,
                bw=BC_BANDWIDTH, gridsize=BC_GRIDSIZE, label=BC_SAMPLE_NAME)
    if ctrl_parsim_acc is not None:
        sns.kdeplot(ctrl_parsim_acc, shade=True,
                    bw=BC_BANDWIDTH, gridsize=BC_GRIDSIZE,
                    label=BC_CONTROL_NAME)
    plt.legend(prop={'size': 16}, title=BC_LEGEND_LABEL)
    plt.xlabel('Mapping Accuracy')
    plt.ylabel('Density')
    plt.title('Mapping Accuracy (Parsimonious: match - ins / ref_len)')
    plt.xlim(PLOT_MIN_BC_ACC, 100)
    pdf_fp.savefig(bbox_inches='tight')
    plt.close()


def report_acc_metrics(res_dir, out_fp):
    try:
        bc_data = mapping.parse_map_summary_file(mh.get_megalodon_fn(
            res_dir, mh.MAP_SUMM_NAME))
        bc_acc = np.array([r_data.pct_identity for r_data in bc_data])
        parsim_acc = np.array([
            100 * (r_data.num_match - r_data.num_ins) /
            (r_data.num_align - r_data.num_ins) for r_data in bc_data])
        mean_bc_acc = np.mean(bc_acc)
        med_bc_acc = np.median(bc_acc)
        # crude mode by rounding to 1 decimal
        uniq_acc, acc_counts = np.unique(np.around(
            bc_acc, 1), return_counts=True)
        mode_bc_acc = uniq_acc[np.argmax(acc_counts)]
        out_fp.write(
            ('Mapping metrics for {} :\t{:.4f}\t{:.4f}\t' +
             '{:.1f}\t{}\n').format(res_dir, med_bc_acc, mean_bc_acc,
                                    mode_bc_acc, len(bc_data)))
    except FileNotFoundError:
        bc_acc = parsim_acc = None
        if VERBOSE:
            sys.stderr.write(
                'WARNING: Mappings not found for {}\n'.format(res_dir))

    return bc_acc, parsim_acc


def parse_mod_data(
        res_dir, out_fp, valid_sites, include_strand, ctrl_sites=None):
    mod_acc, parsim_acc = report_acc_metrics(res_dir, out_fp)

    ctrl_data = None
    mods_db_fn = mh.get_megalodon_fn(res_dir, mh.PR_MOD_NAME)
    try:
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
    except (FileNotFoundError, mh.MegaError):
        mods_data = None

    return mod_acc, parsim_acc, mods_data, ctrl_data


def parse_valid_sites(valid_sites_fns, gt_data_fn, include_strand):
    if valid_sites_fns is None and gt_data_fn is None:
        return None, None, None

    # if ground truth file provided, parse first
    if gt_data_fn is not None:
        if VERBOSE:
            sys.stderr.write('Reading ground truth file\n')
        gt_mod_pos, gt_ctrl_pos = mh.parse_ground_truth_file(
            gt_data_fn, include_strand=include_strand)
        if valid_sites_fns is None:
            # if ground truth provided, but not valid sites return parsed
            # ground truth sites.
            return [gt_mod_pos, ], None, [gt_ctrl_pos, ]

    # parse valid sites files and intersect with ground truth (if provided)
    if VERBOSE:
        sys.stderr.write('Reading valid sites data\n')
    valid_sites, vs_labels = [], []
    ctrl_sites = None if gt_data_fn is None else []
    for vs_label, valid_sites_fn in valid_sites_fns:
        try:
            vs_i_sites = mh.parse_beds([valid_sites_fn, ])
        except FileNotFoundError:
            sys.stderr.write(
                'Could not find valid sites file: {}'.format(
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
        vs_labels.append(vs_label)

    if len(valid_sites) == 0:
        return None, None, None

    return valid_sites, vs_labels, ctrl_sites


def _main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    pdf_fp = PdfPages(args.out_pdf)
    out_fp = (sys.stdout if args.out_filename is None else
              open(args.out_filename, 'w'))
    valid_sites, vs_labels, ctrl_sites = parse_valid_sites(
        args.valid_sites, args.ground_truth_data, args.strand_specific_sites)

    out_fp.write('Mapping metrics: Median Alignment Accuracy :: ' +
                 'Mean Alignment Accuracy :: Mode Alignment Accuracy :: ' +
                 'Num. of Mapped Reads\n')
    if VERBOSE:
        sys.stderr.write('Reading Megalodon data results\n')
    mod_acc, mod_parsim_acc, mods_data, ctrl_data = parse_mod_data(
        args.megalodon_results_dir, out_fp, valid_sites,
        args.strand_specific_sites, ctrl_sites)
    ctrl_acc = ctrl_parsim_acc = None
    # if control is not specified via ground truth file, and control results
    # dir was provided parse control data
    if ctrl_data is None and args.control_megalodon_results_dir is not None:
        if VERBOSE:
            sys.stderr.write('Reading Megalodon control data results\n')
        ctrl_acc, ctrl_parsim_acc, ctrl_data, _ = parse_mod_data(
            args.control_megalodon_results_dir, out_fp, valid_sites,
            args.strand_specific_sites)

    if mod_acc is not None:
        plot_acc(pdf_fp, mod_acc, mod_parsim_acc, ctrl_acc, ctrl_parsim_acc)

    # if modified and control data is not available then mod metrics cannot
    # be computed
    if mods_data is None or ctrl_data is None:
        pdf_fp.close()
        if out_fp is not sys.stdout:
            out_fp.close()
        sys.exit()

    report_mod_metrics(
        mods_data, ctrl_data, not args.allow_unbalance_classes,
        vs_labels, out_fp, pdf_fp)

    pdf_fp.close()
    if out_fp is not sys.stdout:
        out_fp.close()


if __name__ == '__main__':
    _main(get_parser_validate_results().parse_args())
