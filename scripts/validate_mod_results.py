import os
import sys
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from sklearn.metrics import (
    roc_curve, auc, precision_recall_curve, average_precision_score)

from megalodon import megalodon_helper as mh

MAP_FN = mh.OUTPUT_FNS[mh.MAP_NAME][1]
MODS_FN = mh.OUTPUT_FNS[mh.PR_MOD_NAME][1]

VERBOSE = False


parser = argparse.ArgumentParser()
parser.add_argument(
    'megalodon_results_dir',
    help='Output directory from basecall_and_more.py script ' +
    '(run with mappings and per_read_mods outputs).')
parser.add_argument(
    '--ground-truth-data',
    help='Ground truth csv with (chrm, pos, is_mod) values.')
parser.add_argument(
    '--mod-chrms-startswith',
    help='String prefix for all mapped chromosomes with ground ' +
    'truth modifications. All other sites will be assumed unmodified.')
parser.add_argument(
    '--out-pdf', default='mod_base_refactor_validation.pdf',
    help='Output pdf filename.')
parser.add_argument(
    '--verbose', action='store_true',
    help='Output progress information on top of performance metrics.')


def main():
    args = parser.parse_args()
    global VERBOSE
    VERBOSE = args.verbose

    if VERBOSE: print('Reading ground truth data')
    gt_dat, mod_chrm_sw = None, None
    if args.ground_truth_data is not None:
        gt_dat = pd.read_csv(
            args.ground_truth_data, header=None,
            names=['chrm', 'pos', 'is_mod'])
    elif args.mod_chrms_startswith is not None:
        mod_chrm_sw = args.mod_chrms_startswith
    else:
        sys.stderr.write('ERROR: Must provide grounf truth.\n')
        sys.exit(1)

    if VERBOSE: print('Reading megalodon data')
    mega_dat = pd.read_csv(
        os.path.join(args.megalodon_results_dir, MODS_FN), sep='\t', header=None,
        names=['read_id' ,'chrm', 'strand', 'pos', 'score',
               'motif', 'mod_base'])
    mega_dat['chrm'] = mega_dat['chrm'].str.replace('humanGRCh38_', '')
    try:
        bc_dat = pd.read_csv(
            os.path.join(args.megalodon_results_dir, MAP_FN), sep='\t')
        mean_bc_acc = np.mean(bc_dat['pct_identity'])
        med_bc_acc = np.median(bc_dat['pct_identity'])
    except FileNotFoundError:
        mean_bc_acc, med_bc_acc = 0, 0
        if VERBOSE: print('Basecalls not found')

    if VERBOSE: print('Merging ground truth information')
    # merge scores with known mod sites
    if gt_dat is not None:
        m_dat = pd.merge(mega_dat, gt_dat, on=['chrm', 'pos'], sort=False)
    else:
        m_dat = mega_dat
        m_dat['is_mod'] = np.array([
            chrm.startswith(mod_chrm_sw) for chrm in m_dat['chrm']])

    if VERBOSE: print('Modified base class distribution:\n' +
                      str(m_dat['is_mod'].value_counts()))
    if VERBOSE: print('Computing PR/ROC')
    # compute roc and presicion recall
    precision, recall, _ = precision_recall_curve(
        m_dat['is_mod'], -m_dat['score'])
    avg_prcn = average_precision_score(m_dat['is_mod'], -m_dat['score'])

    fpr, tpr, _ = roc_curve(m_dat['is_mod'], -m_dat['score'])
    roc_auc = auc(fpr, tpr)

    if VERBOSE: print('Metrics (bc_mean_acc, bc_med_acc, pr_auc, roc_auc)')
    print('{}\t{:.4f}\t{:.4f}\t{:.6f}\t{:.6f}'.format(
        args.megalodon_results_dir, mean_bc_acc, med_bc_acc, avg_prcn, roc_auc))

    if VERBOSE: print('Plotting')
    with PdfPages(args.out_pdf) as pdf:
        bw = 0.2
        gs = 1000
        plt.figure(figsize=(11, 7))
        sns.kdeplot(m_dat[m_dat['is_mod']]['score'], shade=True,
                    bw=bw, gridsize=gs, label='Yes')
        sns.kdeplot(m_dat[~m_dat['is_mod']]['score'], shade=True,
                    bw=bw, gridsize=gs, label='No')
        plt.legend(prop={'size':16}, title='Is Modified?')
        plt.xlabel('Modified Base Score')
        plt.ylabel('Density')
        pdf.savefig(bbox_inches='tight')
        plt.close()

        plt.figure(figsize=(8, 7))
        plt.step(recall, precision, where='post')
        plt.ylim([-0.05, 1.05])
        plt.xlim([-0.05, 1.05])
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.title('Precision-Recall curve: AP={0:0.2f}'.format(avg_prcn))
        pdf.savefig(bbox_inches='tight')
        plt.close()

        plt.figure(figsize=(8, 7))
        plt.plot(fpr, tpr)
        plt.xlim([-0.05, 1.05])
        plt.ylim([-0.05, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('ROC curve: auc={:0.2f}'.format(roc_auc))
        pdf.savefig(bbox_inches='tight')
        plt.close()

    return

if __name__ == '__main__':
    main()
