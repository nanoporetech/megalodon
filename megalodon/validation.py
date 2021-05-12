import sys

import numpy as np
import seaborn as sns
import matplotlib

if True:
    # Agg appears to be the most robust backend when only saving plots.
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import (
    roc_curve,
    auc,
    precision_recall_curve,
    average_precision_score,
)

from megalodon import logging, megalodon_helper as mh

LOGGER = logging.get_logger()

# BANDWIDTH2 supports seaborn<0.11 when bw_adjust was introduced
MOD_BANDWIDTH = 0.9
MOD_BANDWIDTH2 = 0.2
GRIDSIZE = 1000

MOD_VAL_METRICS_HEADER = (
    "{: <12}{: <19}{: <20}{: <9}{: <20}{: <19}{: <10}{}  {}\n".format(
        "Optimal_F1",
        "Optimal_Threshold",
        "Mean_Avg_Precision",
        "ROC_AUC",
        "Num_Modified_Stats",
        "Num_Control_Stats",
        "Mod_Base",
        "Sample_Label",
        "Valid_Sites_Label",
    )
)
MOD_VAL_METRICS_TMPLT = (
    "{: <12.6f}{: <19.4f}{: <20.6f}{: <9.6f}{: <20d}{: <19d}{: <10}{}  {}\n"
)


def plot_pr(pdf_fp, pr_data):
    for mod_base, mod_pr_data in pr_data.items():
        LOGGER.info("Plotting {} precision-recall curves".format(mod_base))
        plt.figure(figsize=(8, 7))
        for lab, prec, recall in mod_pr_data:
            plt.step(recall, prec, label=lab, where="post")
        plt.ylim([-0.05, 1.05])
        plt.xlim([-0.05, 1.05])
        plt.xlabel("Recall")
        plt.ylabel("Precision")
        plt.title(
            ('Modified Base "{}" Precision-Recall Curves').format(mod_base)
        )
        plt.legend()
        pdf_fp.savefig(bbox_inches="tight")
        plt.close()


def plot_roc(pdf_fp, roc_data):
    for mod_base, mod_roc_data in roc_data.items():
        LOGGER.info("Plotting {} ROC curves".format(mod_base))
        plt.figure(figsize=(8, 7))
        for lab, fpr, tpr in mod_roc_data:
            plt.step(fpr, tpr, label=lab)
        plt.ylim([-0.05, 1.05])
        plt.xlim([-0.05, 1.05])
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.title(('Modified Base "{}" ROC Curves').format(mod_base))
        plt.legend()
        pdf_fp.savefig(bbox_inches="tight")
        plt.close()


def plot_kde(pdf_fp, kde_data):
    for samp_lab, mod_stats, ctrl_stats in kde_data:
        LOGGER.info(
            "Plotting {} modified base statistics densities".format(samp_lab)
        )
        plt.figure(figsize=(8, 5))
        try:
            sns.kdeplot(
                mod_stats,
                shade=True,
                bw_adjust=MOD_BANDWIDTH,
                gridsize=GRIDSIZE,
                label="Yes",
            )
            sns.kdeplot(
                ctrl_stats,
                shade=True,
                bw_adjust=MOD_BANDWIDTH,
                gridsize=GRIDSIZE,
                label="No",
            )
        except AttributeError:
            sns.kdeplot(
                mod_stats,
                shade=True,
                bw=MOD_BANDWIDTH2,
                gridsize=GRIDSIZE,
                label="Yes",
            )
            sns.kdeplot(
                ctrl_stats,
                shade=True,
                bw=MOD_BANDWIDTH2,
                gridsize=GRIDSIZE,
                label="No",
            )
        plt.legend(prop={"size": 16}, title="Is Modified")
        plt.xlabel(
            "Log Likelihood Ratio\nMore Likely Modified <--> "
            + "More Likely Canonical"
        )
        plt.ylabel("Density")
        plt.title(samp_lab)
        pdf_fp.savefig(bbox_inches="tight")
        plt.close()


def compute_mod_sites_stats(
    mod_stats, ctrl_stats, balance_classes, mod_base, samp_lab, vs_lab, out_fp
):
    if balance_classes:
        # randomly downsample sample with more observations
        if mod_stats.shape[0] > ctrl_stats.shape[0]:
            mod_stats = np.random.choice(
                mod_stats, ctrl_stats.shape[0], replace=False
            )
        elif mod_stats.shape[0] < ctrl_stats.shape[0]:
            ctrl_stats = np.random.choice(
                ctrl_stats, mod_stats.shape[0], replace=False
            )

    is_can = np.repeat([0, 1], [mod_stats.shape[0], ctrl_stats.shape[0]])
    all_stats = np.concatenate([mod_stats, ctrl_stats])
    if any(np.isnan(all_stats)):
        LOGGER.warning(
            ("Encountered {} NaN modified base scores.").format(
                sum(np.isnan(all_stats))
            )
        )
        all_stats = all_stats[~np.isnan(all_stats)]
        if all_stats.shape[0] == 0:
            raise mh.MegaError("All modified base scores are NaN")
    inf_idx = np.isinf(all_stats)
    if any(inf_idx):
        all_stats[inf_idx] = np.max(all_stats[~inf_idx])
    neginf_idx = np.isinf(all_stats)
    if any(neginf_idx):
        all_stats[neginf_idx] = np.min(all_stats[~neginf_idx])
    LOGGER.info(
        "Computing PR/ROC for {} from {} at {}".format(
            mod_base, samp_lab, vs_lab
        )
    )
    # compute roc and presicion recall
    precision, recall, thresh = precision_recall_curve(is_can, all_stats)
    prec_recall_sum = precision + recall
    valid_idx = np.where(prec_recall_sum > 0)
    all_f1 = (
        2
        * precision[valid_idx]
        * recall[valid_idx]
        / prec_recall_sum[valid_idx]
    )
    optim_f1_idx = np.argmax(all_f1)
    optim_f1 = all_f1[optim_f1_idx]
    optim_thresh = thresh[optim_f1_idx]
    avg_prcn = average_precision_score(is_can, all_stats)

    fpr, tpr, _ = roc_curve(is_can, all_stats)
    roc_auc = auc(fpr, tpr)

    out_fp.write(
        MOD_VAL_METRICS_TMPLT.format(
            optim_f1,
            optim_thresh,
            avg_prcn,
            roc_auc,
            mod_stats.shape[0],
            ctrl_stats.shape[0],
            mod_base,
            samp_lab,
            vs_lab,
        )
    )
    pr_data = (
        "{} at {} mAP={:0.2f}".format(samp_lab, vs_lab, avg_prcn),
        precision,
        recall,
    )
    roc_data = (
        "{} at {} AUC={:0.2f}".format(samp_lab, vs_lab, roc_auc),
        fpr,
        tpr,
    )
    kde_data = (
        "{} from {} at {}".format(mod_base, samp_lab, vs_lab),
        mod_stats,
        ctrl_stats,
    )

    return pr_data, roc_data, kde_data


if __name__ == "__main__":
    sys.stderr.write("This is a module. See commands with `megalodon -h`")
    sys.exit(1)
