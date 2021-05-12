import sys
from collections import defaultdict

import numpy as np
import matplotlib

if True:
    # Agg appears to be the most robust backend when only saving plots.
    matplotlib.use("Agg")
from matplotlib.backends.backend_pdf import PdfPages

from megalodon import logging, megalodon_helper as mh, mods, validation
from ._extras_parsers import get_parser_calibrate_modified_bases


LOGGER = logging.get_logger()


def extract_llrs(llr_fn):
    llrs_data = np.load(llr_fn)
    mod_bases = llrs_data[mods.GT_ALL_MOD_BASE_STR]
    mod_base_llrs = {}
    for mod_base in mod_bases:
        mod_base_llrs[mod_base] = (
            llrs_data[mods.GT_MOD_LLR_STR.format(mod_base)],
            llrs_data[mods.GT_CAN_LLR_STR.format(mod_base)],
        )

    return mod_base_llrs


def _main(args):
    logging.init_logger()

    LOGGER.info("Parsing log-likelihood ratios")
    mod_base_llrs = extract_llrs(args.ground_truth_llrs)

    out_fp = (
        sys.stdout
        if args.out_filename is None
        else open(args.out_filename, "w")
    )
    out_fp.write(validation.MOD_VAL_METRICS_HEADER)
    pdf_fp = None if args.out_pdf is None else PdfPages(args.out_pdf)
    all_pr_data, all_roc_data = defaultdict(list), defaultdict(list)
    all_kde_data = []
    for mod_base, (mod_llrs, can_llrs) in mod_base_llrs.items():
        LOGGER.info(f'Computing "{mod_base}" modified base validation.')
        try:
            pr_data, roc_data, kde_data = validation.compute_mod_sites_stats(
                mod_llrs,
                can_llrs,
                not args.allow_unbalance_classes,
                mod_base,
                "Megalodon Calibration Data",
                "Sample",
                out_fp,
            )
            all_pr_data[mod_base].append(pr_data)
            all_roc_data[mod_base].append(roc_data)
            all_kde_data.append(kde_data)
        except mh.MegaError as e:
            LOGGER.warning(str(e))
    validation.plot_pr(pdf_fp, all_pr_data)
    validation.plot_roc(pdf_fp, all_roc_data)
    validation.plot_kde(pdf_fp, all_kde_data)
    if pdf_fp is not None:
        pdf_fp.close()


if __name__ == "__main__":
    _main(get_parser_calibrate_modified_bases().parse_args())
