import os
import sys
from collections import defaultdict, namedtuple

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.metrics import (
    roc_curve, auc, precision_recall_curve, average_precision_score)

from megalodon import logging, megalodon_helper as mh
from ._extras_parsers import get_parser_validate_aggregated_modified_bases


LOGGER = logging.get_logger()

MOD_BANDWIDTH = 1.0
MOD_GRIDSIZE = 1000

STRAND_CONV = {'+': 1, '-': -1, '.': None}
IS_MOD_VALS = set(('true', 't', 'on', 'yes', 'y', '1'))

MOD_SAMPLE = namedtuple('MOD_SAMPLE', ('cov', 'mod_cov', 'test_sites'))

MOD_VAL_METRICS_HEADER = (
    '{: <12}{: <19}{: <20}{: <9}{: <20}{: <19}{}\n'.format(
        'Optimal_F1', 'Optimal_Threshold', 'Mean_Avg_Precision', 'ROC_AUC',
        'Num_Modified_Stats', 'Num_Control_Stats', 'Sample'))
MOD_VAL_METRICS_TMPLT = (
    '{: <12.6f}{: <19.4f}{: <20.6f}{: <9.6f}{: <20d}{: <19d}{}\n')


def parse_mod_sample(bm_files, strand_offset, cov_thresh, samp_name):
    cov, mod_cov = mh.parse_bed_methyls(
        bm_files, strand_offset=strand_offset)
    all_cov = np.array([cov for ctg_cov in cov.values()
                        for cov in ctg_cov.values()])
    LOGGER.info('{} coverage median: {:.2f}   mean: {:.2f}'.format(
        samp_name, np.median(all_cov), np.mean(all_cov)))
    test_sites = {}
    for ctg in mod_cov:
        test_sites[ctg] = set(pos for pos, cov in cov[ctg].items()
                              if cov >= cov_thresh)
    return MOD_SAMPLE(cov, mod_cov, test_sites)


def parse_ground_truth_data(gt_csvs):
    def convert_strand(strand_str):
        try:
            return STRAND_CONV[strand_str]
        except KeyError:
            return None

    all_gt_data = {}
    for gt_csv in gt_csvs:
        gt_data = defaultdict(list)
        with open(gt_csv) as gt_fp:
            for line in gt_fp:
                chrm, strand, pos, is_mod = line.strip().split(',')
                gt_data[(chrm, convert_strand(strand))].append((
                    int(pos), is_mod.lower() in IS_MOD_VALS))
        all_gt_data[gt_csv] = dict(gt_data)

    return all_gt_data


def compute_val_metrics(
        mod_samp, ctrl_samp, gt_data, out_fp, pdf_fp, balance_classes,
        ignore_strand, samp_name='sample', valid_pos_fn=None):
    # extract ground truth either from mod and control samples or ground truth
    # data
    if gt_data is None:
        if valid_pos_fn is not None:
            valid_pos = mh.parse_beds(
                [valid_pos_fn, ], ignore_strand=ignore_strand)
            mod_samp = mod_samp._replace(test_sites=dict(
                (ctg, valid_pos[ctg].intersection(ctg_sites))
                for ctg, ctg_sites in mod_samp.test_sites.items()
                if ctg in valid_pos))
            ctrl_samp = ctrl_samp._replace(test_sites=dict(
                (ctg, valid_pos[ctg].intersection(ctg_sites))
                for ctg, ctg_sites in ctrl_samp.test_sites.items()
                if ctg in valid_pos))
        mod_pct_mod = np.array([
            100 * mod_samp.mod_cov[ctg][pos] / mod_samp.cov[ctg][pos]
            for ctg, ctg_poss in mod_samp.test_sites.items()
            for pos in ctg_poss])
        ctrl_pct_mod = np.array(
            [100 * ctrl_samp.mod_cov[ctg][pos] / ctrl_samp.cov[ctg][pos]
             for ctg, ctg_poss in ctrl_samp.test_sites.items()
             for pos in ctg_poss])
    else:
        mod_pct_mod, ctrl_pct_mod = [], []
        for ctg, pos_is_mod in gt_data.items():
            try:
                ctg_cov = mod_samp.cov[ctg]
                ctg_mod_cov = mod_samp.mod_cov[ctg]
            except KeyError:
                continue
            for pos, is_mod in pos_is_mod:
                try:
                    pos_cov = ctg_cov[pos]
                    pos_mod_cov = ctg_mod_cov[pos]
                except KeyError:
                    continue
                if is_mod:
                    mod_pct_mod.append(100 * pos_mod_cov / pos_cov)
                else:
                    ctrl_pct_mod.append(100 * pos_mod_cov / pos_cov)
        mod_pct_mod = np.array(mod_pct_mod)
        ctrl_pct_mod = np.array(ctrl_pct_mod)

    if balance_classes:
        if mod_pct_mod.shape[0] > ctrl_pct_mod.shape[0]:
            mod_pct_mod = np.random.choice(
                mod_pct_mod, ctrl_pct_mod.shape[0], replace=False)
        elif mod_pct_mod.shape[0] < ctrl_pct_mod.shape[0]:
            ctrl_pct_mod = np.random.choice(
                ctrl_pct_mod, mod_pct_mod.shape[0], replace=False)
    all_pct_mod = np.concatenate([mod_pct_mod, ctrl_pct_mod])
    if all_pct_mod.shape[0] == 0:
        LOGGER.info('Skipping "{}". No vaild sites available.'.format(
            samp_name))
        return
    is_mod = np.repeat(
        (1, 0), (mod_pct_mod.shape[0], ctrl_pct_mod.shape[0]))

    precision, recall, thresh = precision_recall_curve(is_mod, all_pct_mod)
    prec_recall_sum = precision + recall
    valid_idx = np.where(prec_recall_sum > 0)
    all_f1 = (2 * precision[valid_idx] * recall[valid_idx] /
              prec_recall_sum[valid_idx])
    optim_f1_idx = np.argmax(all_f1)
    optim_f1 = all_f1[optim_f1_idx]
    optim_thresh = thresh[optim_f1_idx]
    avg_prcn = average_precision_score(is_mod, all_pct_mod)

    fpr, tpr, _ = roc_curve(is_mod, all_pct_mod)
    roc_auc = auc(fpr, tpr)

    out_fp.write(
        MOD_VAL_METRICS_TMPLT.format(
            optim_f1, optim_thresh, avg_prcn, roc_auc, mod_pct_mod.shape[0],
            ctrl_pct_mod.shape[0], samp_name))

    LOGGER.info('Plotting {}'.format(samp_name))
    plt.figure(figsize=(11, 7))
    sns.kdeplot(mod_pct_mod, shade=True, bw=MOD_BANDWIDTH,
                gridsize=MOD_GRIDSIZE, label='Yes')
    sns.kdeplot(ctrl_pct_mod, shade=True, bw=MOD_BANDWIDTH,
                gridsize=MOD_GRIDSIZE, label='No')
    plt.legend(prop={'size': 16}, title='Is Modified?')
    plt.xlabel('Percent Modified')
    plt.ylabel('Density')
    plt.title(samp_name)
    pdf_fp.savefig(bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(8, 7))
    plt.step(recall, precision, where='post')
    plt.ylim([-0.05, 1.05])
    plt.xlim([-0.05, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title(('{}   Precision-Recall curve: AP={:0.2f}').format(
        samp_name, avg_prcn))
    pdf_fp.savefig(bbox_inches='tight')
    plt.close()

    plt.figure(figsize=(8, 7))
    plt.plot(fpr, tpr)
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(('{}   ROC curve: auc={:0.2f}').format(samp_name, roc_auc))
    pdf_fp.savefig(bbox_inches='tight')
    plt.close()


def _main(args):
    logging.init_logger()
    pdf_fp = PdfPages(args.out_pdf)
    out_fp = (sys.stdout if args.out_filename is None else
              open(args.out_filename, 'w'))
    out_fp.write(MOD_VAL_METRICS_HEADER)

    mod_samp = parse_mod_sample(
        args.modified_bed_methyl_files, args.strand_offset,
        args.coverage_threshold, 'Mod')
    ctrl_samp = all_gt_data = None
    if args.control_bed_methyl_files is not None:
        if args.ground_truth_csvs is not None:
            LOGGER.warning(
                'Cannot process both control data and ground truth data. ' +
                'Ignoring ground truth CSV.')
            ctrl_samp = parse_mod_sample(
                args.control_bed_methyl_files, args.strand_offset,
                args.coverage_threshold, 'Control')
    elif args.ground_truth_csvs is not None:
        if args.valid_positions is not None:
            LOGGER.warning(
                'Cannot process both ground truth data and valid sites. ' +
                'Ignoring valid sites.')
            args.valid_positions = None
        all_gt_data = parse_ground_truth_data(args.ground_truth_csvs)
    else:
        LOGGER.error(
            'Must provide either --ground-truth-csvs or ' +
            '--control-bed-methyl-files.')
        sys.exit(1)

    if args.valid_positions is not None:
        for vp_fn in args.valid_positions:
            compute_val_metrics(
                mod_samp, ctrl_samp, None, out_fp, pdf_fp,
                not args.allow_unbalance_classes,
                args.strand_offset is not None,
                os.path.basename(vp_fn), vp_fn)
    elif all_gt_data is not None:
        for gt_fn, gt_data in all_gt_data.items():
            compute_val_metrics(
                mod_samp, ctrl_samp, gt_data,
                out_fp, pdf_fp, not args.allow_unbalance_classes,
                args.strand_offset is not None, os.path.basename(gt_fn))
    else:
        compute_val_metrics(
            mod_samp, ctrl_samp, None, out_fp, pdf_fp,
            not args.allow_unbalance_classes, args.strand_offset is not None)

    pdf_fp.close()
    if out_fp is not sys.stdout:
        out_fp.close()


if __name__ == '__main__':
    _main(get_parser_validate_aggregated_modified_bases().parse_args())
