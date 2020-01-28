import os
import sys
import argparse
import numpy as np
import pandas as pd

import matplotlib
if sys.platform == 'darwin':
    matplotlib.use("TkAgg")
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from sklearn.metrics import (
    roc_curve, auc, precision_recall_curve, average_precision_score)

from megalodon import megalodon_helper as mh


VERBOSE = False

PLOT_MIN_BC_ACC = 80
MOD_BANDWIDTH = 0.2
MOD_GRIDSIZE = 1000
BC_BANDWIDTH = 0.2
BC_GRIDSIZE = 1000

BC_LEGEND_LABEL = 'Sample'
BC_SAMPLE_NAME = 'Sample'
BC_CONTROL_NAME = 'Control'

#BC_LEGEND_LABEL = 'Model'
#BC_SAMPLE_NAME = 'Categorical\nModified Bases\nFlip-Flop'
#BC_CONTROL_NAME = '"High Accuracy"\nFlip-flop'


def compute_mod_sites_stats(
        m_dat, motif, mod_base, v_name, out_fp, pdf_fp, balance_classes):
    motif_m_dat = m_dat[(m_dat['motif'] == motif) &
                        (m_dat['mod_base'] == mod_base)]
    if motif_m_dat.shape[0] == 0:
        # motif/mod_base can conflict with valid sites
        return
    if balance_classes:
        num_mod = sum(motif_m_dat['is_mod'])
        num_unmod = motif_m_dat.shape[0] - num_mod
        if num_mod > num_unmod:
            bal_idx = np.full(motif_m_dat.shape[0], True)
            bal_idx[np.random.choice(
                np.where(motif_m_dat['is_mod'])[0],
                num_mod - num_unmod, replace=False)] = False
            motif_m_dat = motif_m_dat[bal_idx]
        elif num_mod < num_unmod:
            bal_idx = np.full(motif_m_dat.shape[0], True)
            bal_idx[np.random.choice(
                np.where(~motif_m_dat['is_mod'])[0],
                num_unmod - num_mod, replace=False)] = False
            motif_m_dat = motif_m_dat[bal_idx]
    if VERBOSE: sys.stderr.write(
            'Computing PR/ROC for {} in {} at {}\n'.format(
                mod_base, motif, v_name))
    # compute roc and presicion recall
    precision, recall, thresh = precision_recall_curve(
        motif_m_dat['is_mod'], motif_m_dat['llr'])
    prec_recall_sum = precision + recall
    valid_idx = np.where(prec_recall_sum > 0)
    all_f1 = (2 * precision[valid_idx] * recall[valid_idx] /
              prec_recall_sum[valid_idx])
    optim_f1_idx = np.argmax(all_f1)
    optim_f1 = all_f1[optim_f1_idx]
    optim_thresh = thresh[optim_f1_idx]
    avg_prcn = average_precision_score(
        motif_m_dat['is_mod'], motif_m_dat['llr'])

    fpr, tpr, _ = roc_curve(
        motif_m_dat['is_mod'], motif_m_dat['llr'])
    roc_auc = auc(fpr, tpr)

    out_fp.write((
        'Modified base metrics for {} in {} at {}:\t{:.6f} (at {:.4f} )\t' +
        '{:.6f}\t{:.6f}\t{}\t{}\n').format(
            mod_base, motif, v_name, optim_f1, optim_thresh, avg_prcn, roc_auc,
            sum(motif_m_dat['is_mod']), sum(~motif_m_dat['is_mod'])))

    if VERBOSE: sys.stderr.write('Plotting {} in {} at {}\n'.format(
            mod_base, motif, v_name))
    plt.figure(figsize=(11, 7))
    sns.kdeplot(motif_m_dat[motif_m_dat['is_mod']]['llr'],
                shade=True, bw=MOD_BANDWIDTH, gridsize=MOD_GRIDSIZE,
                label='Yes')
    sns.kdeplot(motif_m_dat[~motif_m_dat['is_mod']]['llr'],
                shade=True, bw=MOD_BANDWIDTH, gridsize=MOD_GRIDSIZE,
                label='No')
    plt.legend(prop={'size':16}, title='Is Modified?')
    plt.xlabel('Log Likelihood Ratio\nLess Modified <--> More Modified')
    plt.ylabel('Density')
    plt.title('Modified Base: {}\tMotif: {}\tSites: {}'.format(
        mod_base, motif, v_name))
    pdf_fp.savefig(bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(8, 7))
    plt.step(recall, precision, where='post')
    plt.ylim([-0.05, 1.05])
    plt.xlim([-0.05, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title(('Modified Base: {}\tMotif: {}\tSites: {}\tPrecision-Recall ' +
               'curve: AP={:0.2f}').format(mod_base, motif, v_name, avg_prcn))
    pdf_fp.savefig(bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(8, 7))
    plt.plot(fpr, tpr)
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(('Modified Base: {}\tMotif: {}\tSites: {}\tROC curve: ' +
               'auc={:0.2f}').format(mod_base, motif, v_name, roc_auc))
    pdf_fp.savefig(bbox_inches='tight')
    plt.close()
    return

def report_mod_metrics(
        m_dat, args, out_fp, pdf_fp, valid_sites, balance_classes):
    # cap -inf log probs to lowest other value
    mod_is_ninf = np.isneginf(m_dat['mod_log_prob'])
    m_dat.loc[mod_is_ninf, 'mod_log_prob'] = np.min(
        m_dat['mod_log_prob'][~mod_is_ninf])
    can_is_ninf = np.isneginf(m_dat['can_log_prob'])
    m_dat.loc[can_is_ninf, 'can_log_prob'] = np.min(
        m_dat['can_log_prob'][~can_is_ninf])
    m_dat['llr'] = m_dat['mod_log_prob'] - m_dat['can_log_prob']
    uniq_grps = m_dat.groupby(['mod_base', 'motif']).size().reset_index()

    if valid_sites is not None:
        m_idx = m_dat.set_index(['chrm', 'pos']).index

    out_fp.write('Modified Base Metrics: Optimal F1 :: Optimal Threshold :: ' +
                 'Average Precision :: ROC AUC :: Num. Modified Sites :: ' +
                 'Num. Canonical Sites\n')
    for v_name, valid_sites_i in valid_sites:
        filt_m_dat = m_dat if valid_sites_i is None else m_dat[
            m_idx.isin(valid_sites_i)]
        for mod_base, motif in zip(uniq_grps.mod_base, uniq_grps.motif):
            compute_mod_sites_stats(
                filt_m_dat, motif, mod_base, v_name, out_fp, pdf_fp,
                balance_classes)

    return

def merge_mods_data(mod_dat, ctrl_dat, gt_dat, mod_chrm_sw):
    if VERBOSE: sys.stderr.write('Merging modified base data\n')
    # merge scores with known mod sites
    if ctrl_dat is not None:
        mod_dat['is_mod'] = np.full(mod_dat.shape[0], True)
        ctrl_dat['is_mod'] = np.full(ctrl_dat.shape[0], False)
        m_dat = mod_dat.append(ctrl_dat)
    elif gt_dat is not None:
        m_dat = pd.merge(mod_dat, gt_dat, on=['chrm', 'pos'], sort=False)
    else:
        m_dat = mod_dat
        m_dat['is_mod'] = np.array([
            chrm.startswith(mod_chrm_sw) for chrm in m_dat['chrm']])

    return m_dat

def plot_acc(mod_acc, ctrl_acc, mod_parsim_acc, ctrl_parsim_acc, pdf_fp):
    if VERBOSE: sys.stderr.write('Plotting mapping accuracy distribution(s)\n')
    plt.figure(figsize=(11, 7))
    sns.kdeplot(mod_acc, shade=True,
                bw=BC_BANDWIDTH, gridsize=BC_GRIDSIZE, label=BC_SAMPLE_NAME)
    if ctrl_acc is not None:
        sns.kdeplot(ctrl_acc, shade=True,
                    bw=BC_BANDWIDTH, gridsize=BC_GRIDSIZE,
                    label=BC_CONTROL_NAME)
    plt.legend(prop={'size':16}, title=BC_LEGEND_LABEL)
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
    plt.legend(prop={'size':16}, title=BC_LEGEND_LABEL)
    plt.xlabel('Mapping Accuracy')
    plt.ylabel('Density')
    plt.title('Mapping Accuracy (Parsimonious: match - ins / ref_len)')
    plt.xlim(PLOT_MIN_BC_ACC, 100)
    pdf_fp.savefig(bbox_inches='tight')
    plt.close()

    return

def report_acc_metrics(res_dir, out_fp):
    try:
        bc_dat = pd.read_csv(mh.get_megalodon_fn(res_dir, mh.MAP_SUMM_NAME),
                             sep='\t')
        bc_acc = bc_dat['pct_identity']
        parsim_acc = 100 * (bc_dat['num_match'] - bc_dat['num_ins']) / \
                     (bc_dat['num_align'] - bc_dat['num_ins'])
        mean_bc_acc = np.mean(bc_acc)
        med_bc_acc = np.median(bc_acc)
        # crude mode by rounding to 1 decimal
        uniq_acc, acc_counts = np.unique(np.around(
            bc_acc, 1), return_counts=True)
        mode_bc_acc = uniq_acc[np.argmax(acc_counts)]
        out_fp.write(
            ('Mapping metrics for {} :\t{:.4f}\t{:.4f}\t' +
             '{:.1f}\t{}\n').format(res_dir, med_bc_acc, mean_bc_acc,
                                    mode_bc_acc, bc_dat.shape[0]))
    except FileNotFoundError:
        bc_acc = parsim_acc = None
        if VERBOSE: sys.stderr.write(
                'WARNING: Mappings not found for {}\n'.format(res_dir))

    return bc_acc, parsim_acc

def parse_mod_data(args, out_fp):
    if VERBOSE: sys.stderr.write('Reading megalodon data\n')
    mod_acc, parsim_acc = report_acc_metrics(args.megalodon_results_dir, out_fp)

    try:
        mod_dat = pd.read_csv(
            mh.get_megalodon_fn(args.megalodon_results_dir,
                                mh.PR_MOD_TXT_NAME), sep='\t')
    except FileNotFoundError:
        mod_dat = None

    return mod_dat, mod_acc, parsim_acc

def parse_control_mods(args, out_fp):
    ctrl_acc = ctrl_parsim_acc = ctrl_dat = gt_dat = mod_chrm_sw = None
    if args.control_megalodon_results_dir is not None:
        if VERBOSE: sys.stderr.write('Reading control mods data\n')
        ctrl_acc, ctrl_parsim_acc = report_acc_metrics(
            args.control_megalodon_results_dir, out_fp)
        try:
            ctrl_dat = pd.read_csv(
                mh.get_megalodon_fn(args.control_megalodon_results_dir,
                                    mh.PR_MOD_TXT_NAME), sep='\t')
        except FileNotFoundError:
            ctrl_dat = None

    elif args.ground_truth_data is not None:
        if VERBOSE: sys.stderr.write('Reading ground truth data\n')
        gt_dat = pd.read_csv(
            args.ground_truth_data, header=None,
            names=['chrm', 'pos', 'is_mod'])
    elif args.mod_chrms_startswith is not None:
        mod_chrm_sw = args.mod_chrms_startswith
    else:
        sys.stderr.write(
            '*' * 20 + '  WARNING: No modified base ground truth provided.\n')

    return ctrl_acc, ctrl_parsim_acc, ctrl_dat, gt_dat, mod_chrm_sw

def parse_valid_sites(valid_sites_arg):
    if valid_sites_arg is None:
        return [('', None),]
    if VERBOSE: sys.stderr.write('Reading valid sites data\n')
    valid_sites = []
    for name, valid_sites_fn in valid_sites_arg:
        try:
            valid_sites.append((name, pd.read_csv(
                valid_sites_fn, header=None, sep='\t',
                names=['chrm', 'pos', 'strand'],
                usecols=[0, 1, 5]).set_index(['chrm', 'pos']).index))
        except FileNotFoundError:
            print('Could not find valid sites file: {}'.format(valid_sites_fn))
            continue
    if len(valid_sites) == 0:
        return [('', None),]
    return valid_sites

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'megalodon_results_dir',
        help='Output directory from megalodon with mappings and per_read_mods ' +
        'in outputs. Must have --write-mods-text set for mods validation.')
    parser.add_argument(
        '--control-megalodon-results-dir',
        help='Megalodon output directory with modified base control sample.')
    parser.add_argument(
        '--ground-truth-data',
        help='Ground truth csv with (chrm, pos, is_mod) values.')
    parser.add_argument(
        '--mod-chrms-startswith',
        help='String prefix for all mapped chromosomes with ground ' +
        'truth modifications. All other sites will be assumed unmodified.')
    parser.add_argument(
        '--valid-sites', nargs=2, action='append',
        help='Name and BED file containing sites over which to restrict ' +
        'modified base results. Multiple sets of valid sites may be provided.')
    parser.add_argument(
        '--out-pdf', default='megalodon_validation.pdf',
        help='Output pdf filename. Default: %(default)s')
    parser.add_argument(
        '--out-filename',
        help='Output filename for text summary. Default: stdout')
    parser.add_argument(
        '--balance-classes', action='store_true',
        help='Balance size of modified and canonical classes for each ' +
        'comparison made.')
    parser.add_argument(
        '--quiet', action='store_true',
        help='Suppress progress information.')

    return parser

def main():
    args = get_parser().parse_args()
    global VERBOSE
    VERBOSE = not args.quiet
    pdf_fp = PdfPages(args.out_pdf)
    out_fp = sys.stdout if args.out_filename is None else \
             open(args.out_filename, 'w')
    valid_sites = parse_valid_sites(args.valid_sites)

    out_fp.write('Mapping metrics: Median Alignment Accuracy :: ' +
                 'Mean Alignment Accuracy :: Mode Alignment Accuracy :: ' +
                 'Num. of Mapped Reads\n')
    mod_dat, mod_acc, mod_parsim_acc = parse_mod_data(args, out_fp)

    ctrl_acc, ctrl_parsim_acc, ctrl_dat, gt_dat, mod_chrm_sw \
        = parse_control_mods(args, out_fp)
    if mod_acc is not None:
        plot_acc(mod_acc, ctrl_acc, mod_parsim_acc, ctrl_parsim_acc, pdf_fp)
    # could just compute mapping metrics
    if all(d is None for d in (ctrl_dat, gt_dat, mod_chrm_sw)):
        pdf_fp.close()
        if out_fp is not sys.stdout: out_fp.close()
        return
    m_dat = merge_mods_data(mod_dat, ctrl_dat, gt_dat, mod_chrm_sw)
    report_mod_metrics(
        m_dat, args, out_fp, pdf_fp, valid_sites, args.balance_classes)
    pdf_fp.close()
    if out_fp is not sys.stdout: out_fp.close()

    return

if __name__ == '__main__':
    main()
