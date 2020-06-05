import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from megalodon import calibration, logging, megalodon_helper as mh, mods
from ._extras_parsers import get_parser_calibrate_modified_bases


LOGGER = logging.get_logger()
PROB_COLORS = ("#bcbddc", "#807dba", "#6a51a3")


def plot_calib(
        pdf_fp, mod_base, smooth_ls, s_ref, sm_ref, s_alt, sm_alt,
        mono_prob, prob_alt, prob_threshs, add_prob_thresh):
    f, axarr = plt.subplots(3, sharex=True, figsize=(11, 7))
    axarr[0].plot(smooth_ls, s_ref, color='orange')
    axarr[0].plot(smooth_ls, sm_ref, color='red')
    axarr[0].plot(smooth_ls, s_alt, color='grey')
    axarr[0].plot(smooth_ls, sm_alt, color='blue')
    axarr[0].set_ylabel(
        'Probability Density\nred/orange=canonical\nblue/grey=modified')
    axarr[0].set_title(mod_base + ' Calibration')
    axarr[1].plot(smooth_ls, mono_prob, color='orange')
    axarr[1].plot(
        smooth_ls, 1 / (np.exp(smooth_ls) + 1), color='purple')
    axarr[1].set_ylabel(
        'Emperical Modified\nProbability\norange=calibrated\npurple=raw')
    axarr[2].plot(
        smooth_ls, np.log((1 - prob_alt) / prob_alt), color='red')
    axarr[2].plot(
        smooth_ls, np.log((1 - mono_prob) / mono_prob), color='orange')
    axarr[2].set_ylabel('Calibrated LLR\norage=monotonic')
    axarr[2].set_xlabel('Theoretical LLR (NN Score)')
    if add_prob_thresh:
        # indicate the cutoff points for several common cutoff locations
        thresh_f = np.log((1 - mono_prob) / mono_prob)
        for p, col in zip(prob_threshs, PROB_COLORS):
            llr_x = np.log(p / (1 - p))
            thresh_val = np.argmin(np.abs(thresh_f - llr_x))
            nthresh_val = np.argmin(np.abs(thresh_f + llr_x))
            prop_filt = (sum(sm_ref[nthresh_val:thresh_val]) +
                         sum(sm_alt[nthresh_val:thresh_val])) / (
                             sum(sm_ref) + sum(sm_alt))
            for i in range(2):
                axarr[i].axvline(x=smooth_ls[thresh_val], color=col)
                axarr[i].axvline(x=smooth_ls[nthresh_val], color=col)
            axarr[2].axvline(x=smooth_ls[thresh_val], color=col)
            axarr[2].axvline(
                x=smooth_ls[nthresh_val], color=col, label=(
                    '--mod-binary-threshold={} (filters {:.0f}%)').format(
                        p, 100 * prop_filt))
            axarr[2].legend(fontsize='small')

    pdf_fp.savefig(bbox_inches='tight')
    plt.close()


def extract_llrs(llr_fn):
    llrs_data = np.load(llr_fn)
    mod_bases = llrs_data[mods.GT_ALL_MOD_BASE_STR]
    mod_base_llrs = {}
    for mod_base in mod_bases:
        mod_base_llrs[mod_base] = (
            llrs_data[mods.GT_MOD_LLR_STR.format(mod_base)],
            llrs_data[mods.GT_CAN_LLR_STR.format(mod_base)])

    return mod_base_llrs


def _main(args):
    logging.init_logger()
    mh.prep_out_fn(args.out_filename, args.overwrite)

    LOGGER.info('Parsing log-likelihood ratios')
    mod_base_llrs = extract_llrs(args.ground_truth_llrs)

    pdf_fp = None if args.out_pdf is None else PdfPages(args.out_pdf)
    save_kwargs = {}
    for mod_base, (mod_llrs, can_llrs) in mod_base_llrs.items():
        LOGGER.info(
            'Computing {} modified base calibration.'.format(mod_base))
        mod_calib, mod_llr_range, plot_data = calibration.compute_calibration(
            can_llrs, mod_llrs, args.max_input_llr,
            args.num_calibration_values, args.smooth_bandwidth,
            args.min_density, args.diff_epsilon, args.llr_clip_buffer,
            pdf_fp is not None, num_proc=args.processes)
        save_kwargs[mod_base + calibration.LLR_RANGE_SUFFIX] = mod_llr_range
        save_kwargs[mod_base + calibration.CALIB_TABLE_SUFFIX] = mod_calib
        if pdf_fp is not None:
            plot_calib(pdf_fp, mod_base, *plot_data, args.pdf_prob_thresholds,
                       not args.plot_without_prob_thresholds)
    if pdf_fp is not None:
        pdf_fp.close()

    # save calibration table for reading into mod calibration table
    LOGGER.info('Saving calibrations to file.')
    mod_bases = list(mod_base_llrs.keys())
    np.savez(
        args.out_filename,
        stratify_type=calibration.MOD_BASE_STRAT_TYPE,
        smooth_nvals=args.num_calibration_values,
        mod_bases=mod_bases,
        **save_kwargs)


if __name__ == '__main__':
    _main(get_parser_calibrate_modified_bases().parse_args())
