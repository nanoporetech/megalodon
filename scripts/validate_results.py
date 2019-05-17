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

MAP_FN = mh.OUTPUT_FNS[mh.MAP_NAME][1]
MODS_FN = mh.OUTPUT_FNS[mh.PR_MOD_NAME][1]

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
    '--out-pdf', default='megalodon_validation.pdf',
    help='Output pdf filename. Default: %(default)s')
parser.add_argument(
    '--out-filename',
    help='Output filename for text summary. Default: stdout')
parser.add_argument(
    '--quiet', action='store_true',
    help='Suppress progress information.')


def report_mod_metrics(m_dat, args, out_fp, pdf_fp):
    for motif in np.unique(m_dat['motif']):
        motif_m_dat = m_dat[m_dat['motif'] == motif]
        if VERBOSE: sys.stderr.write('Computing PR/ROC for {}\n'.format(motif))
        out_fp.write(('Modified base class distribution for {}:\n\t' +
                      'Number of modified observations:    {}\n\t' +
                      'Number of unmodified observations:  {}\n').format(
                          motif,
                          sum(motif_m_dat['is_mod']),
                          sum(~motif_m_dat['is_mod'])))
        # compute roc and presicion recall
        precision, recall, _ = precision_recall_curve(
            motif_m_dat['is_mod'], -motif_m_dat['score'])
        avg_prcn = average_precision_score(
            motif_m_dat['is_mod'], -motif_m_dat['score'])

        fpr, tpr, _ = roc_curve(motif_m_dat['is_mod'], -motif_m_dat['score'])
        roc_auc = auc(fpr, tpr)

        out_fp.write(('Modified base metrics for {}:\n\t' +
                      'Average precision:  {:.6f}\n\t' +
                      'ROC AUC:            {:.6f}\n').format(
                          motif, avg_prcn, roc_auc))

        if VERBOSE: sys.stderr.write('Plotting {}\n'.format(motif))
        plt.figure(figsize=(11, 7))
        sns.kdeplot(motif_m_dat[motif_m_dat['is_mod']]['score'], shade=True,
                    bw=MOD_BANDWIDTH, gridsize=MOD_GRIDSIZE, label='Yes')
        sns.kdeplot(motif_m_dat[~motif_m_dat['is_mod']]['score'], shade=True,
                    bw=MOD_BANDWIDTH, gridsize=MOD_GRIDSIZE, label='No')
        plt.legend(prop={'size':16}, title='Is Modified?')
        plt.xlabel('Modified Base Score')
        plt.ylabel('Density')
        plt.title('Motif: {}'.format(motif))
        pdf_fp.savefig(bbox_inches='tight')
        plt.close()

        plt.figure(figsize=(8, 7))
        plt.step(recall, precision, where='post')
        plt.ylim([-0.05, 1.05])
        plt.xlim([-0.05, 1.05])
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.title('Motif: {}\tPrecision-Recall curve: AP={:0.2f}'.format(
            motif, avg_prcn))
        pdf_fp.savefig(bbox_inches='tight')
        plt.close()

        plt.figure(figsize=(8, 7))
        plt.plot(fpr, tpr)
        plt.xlim([-0.05, 1.05])
        plt.ylim([-0.05, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Motif: {}\tROC curve: auc={:0.2f}'.format(motif, roc_auc))
        pdf_fp.savefig(bbox_inches='tight')
        plt.close()

    return

def merge_mods_data(mod_dat, ctrl_dat, gt_dat, mod_chrm_sw, out_fp):
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
        bc_dat = pd.read_csv(
            os.path.join(res_dir, MAP_FN), sep='\t')
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
            ('Mapping metrics for {} ({} mapped reads):\n\t' +
             'Mean Pct. Identity:    {:.4f}\n\t' +
             'Median Pct. Identity:  {:.4f}\n\t' +
             'Mode Pct. Identity:    {:.1f}\n').format(
                 res_dir, bc_dat.shape[0],
                 mean_bc_acc, med_bc_acc, mode_bc_acc))
    except FileNotFoundError:
        bc_acc = None
        if VERBOSE: sys.stderr.write(
                '*' * 20 + 'Mappings not found for {}\n'.format(res_dir))

    return bc_acc, parsim_acc

def parse_mod_data(args, out_fp):
    if VERBOSE: sys.stderr.write('Reading megalodon data\n')
    mod_acc, parsim_acc = report_acc_metrics(args.megalodon_results_dir, out_fp)

    try:
        mod_dat = pd.read_csv(
            os.path.join(args.megalodon_results_dir, MODS_FN),
            sep='\t', header=None,
            names=['read_id' ,'chrm', 'strand', 'pos', 'score',
                   'motif', 'mod_base'])
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
                os.path.join(args.control_megalodon_results_dir, MODS_FN),
                sep='\t', header=None,
                names=['read_id' ,'chrm', 'strand', 'pos', 'score',
                       'motif', 'mod_base'])
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

def main():
    args = parser.parse_args()
    global VERBOSE
    VERBOSE = not args.quiet
    pdf_fp = PdfPages(args.out_pdf)
    out_fp = sys.stdout if args.out_filename is None else \
             open(args.out_filename, 'w')

    mod_dat, mod_acc, mod_parsim_acc = parse_mod_data(args, out_fp)

    # TODO add SNP validation (a bit more complicated)

    ctrl_acc, ctrl_parsim_acc, ctrl_dat, gt_dat, mod_chrm_sw \
        = parse_control_mods(args, out_fp)
    plot_acc(mod_acc, ctrl_acc, mod_parsim_acc, ctrl_parsim_acc, pdf_fp)
    # could just compute mapping metrics
    if all(d is None for d in (ctrl_dat, gt_dat, mod_chrm_sw)):
        pdf_fp.close()
        if out_fp is not sys.stdout: out_fp.close()
        return
    m_dat = merge_mods_data(mod_dat, ctrl_dat, gt_dat, mod_chrm_sw, out_fp)
    report_mod_metrics(m_dat, args, out_fp, pdf_fp)
    pdf_fp.close()
    if out_fp is not sys.stdout: out_fp.close()

    return

if __name__ == '__main__':
    main()
