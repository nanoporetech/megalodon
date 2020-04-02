import os
import sys
import argparse

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.metrics import (
    roc_curve, auc, precision_recall_curve, average_precision_score)

from megalodon import megalodon_helper as mh


VERBOSE = False

MOD_BANDWIDTH = 1.0
MOD_GRIDSIZE = 1000


def compute_val_metrics(
        mod_cov, mod_mod_cov, mod_test_sites,
        ctrl_cov, ctrl_mod_cov, ctrl_test_sites,
        out_fp, pdf_fp, balance_classes, ignore_strand, valid_pos_fn=None):
    samp = 'sample'
    if valid_pos_fn is not None:
        samp = os.path.basename(valid_pos_fn)
        valid_pos = mh.parse_beds(
            [valid_pos_fn, ], ignore_strand=ignore_strand)
        mod_test_sites = dict((ctg, valid_pos[ctg].intersection(ctg_sites))
                              for ctg, ctg_sites in mod_test_sites.items()
                              if ctg in valid_pos)
        ctrl_test_sites = dict((ctg, valid_pos[ctg].intersection(ctg_sites))
                               for ctg, ctg_sites in ctrl_test_sites.items()
                               if ctg in valid_pos)

    mod_pct_meths = np.array([100 * mod_mod_cov[ctg][pos] / mod_cov[ctg][pos]
                              for ctg, ctg_poss in mod_test_sites.items()
                              for pos in ctg_poss])
    ctrl_pct_meths = np.array(
        [100 * ctrl_mod_cov[ctg][pos] / ctrl_cov[ctg][pos]
         for ctg, ctg_poss in ctrl_test_sites.items()
         for pos in ctg_poss])
    if balance_classes:
        if mod_pct_meths.shape[0] > ctrl_pct_meths.shape[0]:
            mod_pct_meths = np.random.choice(
                mod_pct_meths, ctrl_pct_meths.shape[0], replace=False)
        elif mod_pct_meths.shape[0] < ctrl_pct_meths.shape[0]:
            ctrl_pct_meths = np.random.choice(
                ctrl_pct_meths, mod_pct_meths.shape[0], replace=False)
    pct_meths = np.concatenate([mod_pct_meths, ctrl_pct_meths])
    is_mod = np.repeat(
        (1, 0), (mod_pct_meths.shape[0], ctrl_pct_meths.shape[0]))
    if is_mod.shape[0] == 0:
        sys.stderr.write('Skipping "{}". No vaild sites available.\n'.format(
            samp))
        return

    precision, recall, thresh = precision_recall_curve(is_mod, pct_meths)
    prec_recall_sum = precision + recall
    valid_idx = np.where(prec_recall_sum > 0)
    all_f1 = (2 * precision[valid_idx] * recall[valid_idx] /
              prec_recall_sum[valid_idx])
    optim_f1_idx = np.argmax(all_f1)
    optim_f1 = all_f1[optim_f1_idx]
    optim_thresh = thresh[optim_f1_idx]
    avg_prcn = average_precision_score(is_mod, pct_meths)

    fpr, tpr, _ = roc_curve(is_mod, pct_meths)
    roc_auc = auc(fpr, tpr)

    out_fp.write((
        'Modified base metrics for {}:\t{:.6f} (at {:.4f} )\t' +
        '{:.6f}\t{:.6f}\t{}\t{}\n').format(
            samp, optim_f1, optim_thresh, avg_prcn, roc_auc,
            mod_pct_meths.shape[0], ctrl_pct_meths.shape[0]))

    if VERBOSE:
        sys.stderr.write('Plotting {}\n'.format(samp))
    plt.figure(figsize=(11, 7))
    sns.kdeplot(mod_pct_meths, shade=True, bw=MOD_BANDWIDTH,
                gridsize=MOD_GRIDSIZE, label='Yes')
    sns.kdeplot(ctrl_pct_meths, shade=True, bw=MOD_BANDWIDTH,
                gridsize=MOD_GRIDSIZE, label='No')
    plt.legend(prop={'size': 16}, title='Is Modified?')
    plt.xlabel('Percent Methylated')
    plt.ylabel('Density')
    plt.title(samp)
    pdf_fp.savefig(bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(8, 7))
    plt.step(recall, precision, where='post')
    plt.ylim([-0.05, 1.05])
    plt.xlim([-0.05, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title(('{}   Precision-Recall curve: AP={:0.2f}').format(
        samp, avg_prcn))
    pdf_fp.savefig(bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(8, 7))
    plt.plot(fpr, tpr)
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(('{}   ROC curve: auc={:0.2f}').format(samp, roc_auc))
    pdf_fp.savefig(bbox_inches='tight')
    plt.close()


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--modified-bed-methyl-files', nargs='+', required=True,
        help='Bed methyl files from modified sample(s).')
    # TODO add option to provide ground truth file within a sample
    parser.add_argument(
        '--control-bed-methyl-files', nargs='+', required=True,
        help='Bed methyl files from control sample(s).')
    parser.add_argument(
        '--valid-positions', action='append',
        help='BED file containing positions to be considered. Multiple ' +
        'files may be provided')
    parser.add_argument(
        '--coverage-threshold', type=int, default=1,
        help='Only include sites with sufficient coverage. ' +
        'Default: 1 (= All sites)')
    parser.add_argument(
        '--strand-offset', type=int,
        help='Offset to combine stranded results. Positive value indicates ' +
        'reverse strand sites have higher position values. Default treat ' +
        'strands independently.')
    parser.add_argument(
        '--allow-unbalance-classes', action='store_true',
        help='Allow unbalanced classes in modified base metric computation. ' +
        'Default: Balance size of modified and canonical classes for each ' +
        'comparison made.')
    parser.add_argument(
        '--out-pdf', default='megalodon_agg_validation.pdf',
        help='Output pdf filename. Default: %(default)s')
    parser.add_argument(
        '--out-filename',
        help='Output filename for text summary. Default: stdout')

    return parser


def main():
    args = get_parser().parse_args()
    pdf_fp = PdfPages(args.out_pdf)
    out_fp = (sys.stdout if args.out_filename is None else
              open(args.out_filename, 'w'))

    mod_cov, mod_mod_cov = mh.parse_bed_methyls(
        args.modified_bed_methyl_files, strand_offset=args.strand_offset)
    mod_all_cov = np.array([cov for ctg_cov in mod_cov.values()
                            for cov in ctg_cov.values()])
    ctrl_cov, ctrl_mod_cov = mh.parse_bed_methyls(
        args.control_bed_methyl_files, strand_offset=args.strand_offset)
    ctrl_all_cov = np.array([cov for ctg_cov in ctrl_cov.values()
                            for cov in ctg_cov.values()])
    sys.stderr.write('Mod coverage median: {:.2f}   mean: {:.2f}\n'.format(
        np.median(mod_all_cov), np.mean(mod_all_cov)))
    sys.stderr.write('Control coverage median: {:.2f}   mean: {:.2f}\n'.format(
        np.median(ctrl_all_cov), np.mean(ctrl_all_cov)))
    mod_test_sites = {}
    for ctg in mod_cov:
        mod_test_sites[ctg] = set(pos for pos, cov in mod_cov[ctg].items()
                                  if cov >= args.coverage_threshold)
    ctrl_test_sites = {}
    for ctg in ctrl_cov:
        ctrl_test_sites[ctg] = set(pos for pos, cov in ctrl_cov[ctg].items()
                                   if cov >= args.coverage_threshold)

    if args.valid_positions is None:
        compute_val_metrics(
            mod_cov, mod_mod_cov, mod_test_sites,
            ctrl_cov, ctrl_mod_cov, ctrl_test_sites,
            out_fp, pdf_fp, not args.allow_unbalance_classes,
            args.strand_offset is not None)
    else:
        for fn in args.valid_positions:
            compute_val_metrics(
                mod_cov, mod_mod_cov, mod_test_sites,
                ctrl_cov, ctrl_mod_cov, ctrl_test_sites,
                out_fp, pdf_fp, not args.allow_unbalance_classes,
                args.strand_offset is not None, fn)

    pdf_fp.close()
    if out_fp is not sys.stdout:
        out_fp.close()


if __name__ == '__main__':
    main()
