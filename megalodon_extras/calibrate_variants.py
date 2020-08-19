import os
import sys
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from megalodon import calibration, logging, megalodon_helper as mh
from ._extras_parsers import get_parser_calibrate_variants


LOGGER = logging.get_logger()
INVALID_CALIB_MSG = (
    'Encountered invalid distributions for calibration. Not saving ' +
    'calibration file, but pdf will be plotted in order to identify ' +
    'potential issues.')


def plot_calib(
        pdf_fp, var_type, smooth_ls, s_ref, sm_ref, s_alt, sm_alt,
        mono_prob, prob_alt):
    f, axarr = plt.subplots(3, sharex=True, figsize=(11, 7))
    axarr[0].plot(smooth_ls, s_ref, color='orange')
    axarr[0].plot(smooth_ls, sm_ref, color='red')
    axarr[0].plot(smooth_ls, s_alt, color='grey')
    axarr[0].plot(smooth_ls, sm_alt, color='blue')
    axarr[0].set_ylabel(
        'Probability Density\nred/orange=canonical\nblue/grey=modified')
    axarr[0].set_title(var_type + ' Calibration')
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
    pdf_fp.savefig(bbox_inches='tight')
    plt.close()


def extract_llrs(llr_fn, max_indel_len=None):
    snp_ref_llrs, ins_ref_llrs, del_ref_llrs = (
        defaultdict(list) for _ in range(3))
    with open(llr_fn) as llr_fp:
        for line in llr_fp:
            is_ref_correct, llr, ref_seq, alt_seq = line.split()
            llr = float(llr)
            if is_ref_correct != 'True':
                continue
            if np.isnan(llr):
                continue
            if max_indel_len is not None and \
               np.abs(len(ref_seq) - len(alt_seq)) > max_indel_len:
                continue
            if len(ref_seq) == 1 and len(alt_seq) == 1:
                snp_ref_llrs[(ref_seq, alt_seq)].append(llr)
            else:
                if len(ref_seq) > len(alt_seq):
                    del_ref_llrs[len(ref_seq) - len(alt_seq)].append(llr)
                else:
                    ins_ref_llrs[len(alt_seq) - len(ref_seq)].append(llr)

    if min(len(snp_ref_llrs), len(ins_ref_llrs), len(del_ref_llrs)) == 0:
        raise mh.MegaError(
            'Variant statistics file does not contain sufficient data for ' +
            'calibration.')

    return snp_ref_llrs, ins_ref_llrs, del_ref_llrs


def prep_out(out_fn, overwrite):
    if os.path.exists(out_fn):
        if overwrite:
            os.remove(out_fn)
        else:
            raise NotImplementedError(
                'ERROR: --out-filename exists and --overwrite not set.')
    try:
        open(out_fn, 'w').close()
        os.remove(out_fn)
    except Exception:
        LOGGER.error(
            'Attempt to write to --out-filename location failed with the ' +
            'following error.')
        raise


def _main(args):
    logging.init_logger()
    prep_out(args.out_filename, args.overwrite)

    LOGGER.info('Parsing log-likelihood ratios')
    snp_ref_llrs, ins_ref_llrs, del_ref_llrs = extract_llrs(
        args.ground_truth_llrs)
    # add calibration for a generic varaint (mostly multiple SNPs
    # as single variant; but not an indel)
    generic_var_llrs = [llr for snp_type_llrs in snp_ref_llrs.values()
                        for llr in snp_type_llrs]
    # downsample to same level as other snp types
    snp_ref_llrs[
        (calibration.GENERIC_BASE,
         calibration.GENERIC_BASE)] = np.random.choice(
             generic_var_llrs, int(len(generic_var_llrs) / 12), replace=False)
    max_indel_len = max(ins_ref_llrs)
    assert set(ins_ref_llrs) == set(del_ref_llrs), (
        'Must test same range of lengths for insertions and deletions')
    assert set(ins_ref_llrs) == set(range(1, max_indel_len + 1)), (
        'Must test every length in length range for indels')

    do_save_calib = True
    pdf_fp = None if args.out_pdf is None else PdfPages(args.out_pdf)
    LOGGER.info('Computing stratified single-base SNP calibration.')
    snp_calibs = {}
    for (ref_seq, alt_seq), snp_llrs in sorted(snp_ref_llrs.items()):
        LOGGER.info('Computing ' + ref_seq + ' -> ' + alt_seq +
                    ' SNP calibration.')
        try:
            snp_calib, snp_llr_range, plot_data \
                = calibration.compute_mirrored_calibration(
                    np.array(snp_llrs), args.max_input_llr,
                    args.num_calibration_values, args.smooth_bandwidth,
                    args.min_density, args.diff_epsilon, args.llr_clip_buffer,
                    pdf_fp is not None, num_proc=args.processes)
        except mh.MegaError:
            do_save_calib = False
            LOGGER.error(INVALID_CALIB_MSG)
        snp_calibs[(ref_seq, alt_seq)] = (snp_calib, snp_llr_range)
        if pdf_fp is not None:
            plot_calib(pdf_fp, 'SNP: ' + ref_seq + ' -> ' + alt_seq,
                       *plot_data)
    LOGGER.info('Computing deletion calibration.')
    del_calibs = {}
    for del_len, del_llrs in sorted(del_ref_llrs.items()):
        LOGGER.info('Computing deletion length {} calibration.'.format(
            del_len))
        try:
            del_calib, del_llr_range, plot_data \
                = calibration.compute_mirrored_calibration(
                    np.array(del_llrs), args.max_input_llr,
                    args.num_calibration_values, args.smooth_bandwidth,
                    args.min_density, args.diff_epsilon, args.llr_clip_buffer,
                    pdf_fp is not None, num_proc=args.processes)
        except mh.MegaError:
            do_save_calib = False
            LOGGER.error(INVALID_CALIB_MSG)
        del_calibs[del_len] = (del_calib, del_llr_range)
        if pdf_fp is not None:
            plot_calib(pdf_fp, 'Deletion Length ' + str(del_len), *plot_data)
    LOGGER.info('Computing insertion calibration.')
    ins_calibs = {}
    for ins_len, ins_llrs in sorted(ins_ref_llrs.items()):
        LOGGER.info('Computing insertion length {} calibration.'.format(
            ins_len))
        try:
            ins_calib, ins_llr_range, plot_data \
                = calibration.compute_mirrored_calibration(
                    np.array(ins_llrs), args.max_input_llr,
                    args.num_calibration_values, args.smooth_bandwidth,
                    args.min_density, args.diff_epsilon, args.llr_clip_buffer,
                    pdf_fp is not None, num_proc=args.processes)
        except mh.MegaError:
            do_save_calib = False
            LOGGER.error(INVALID_CALIB_MSG)
        ins_calibs[ins_len] = (ins_calib, ins_llr_range)
        if pdf_fp is not None:
            plot_calib(pdf_fp, 'Insertion Length ' + str(ins_len), *plot_data)

    if pdf_fp is not None:
        pdf_fp.close()

    if not do_save_calib:
        sys.exit(1)
    # save calibration table for reading into variant calibration table
    LOGGER.info('Saving calibrations to file.')
    snp_llr_range_save_data, snp_calib_save_data = {}, {}
    for (ref_seq, alt_seq), (snp_calib, snp_llr_range) in snp_calibs.items():
        snp_calib_save_data[
            calibration.SNP_CALIB_TMPLT.format(ref_seq, alt_seq)] = snp_calib
        snp_llr_range_save_data[
            calibration.SNP_LLR_RNG_TMPLT.format(
                ref_seq, alt_seq)] = snp_llr_range
    del_llr_range_save_data, del_calib_save_data = {}, {}
    for del_len, (del_calib, del_llr_range) in del_calibs.items():
        del_calib_save_data[
            calibration.DEL_CALIB_TMPLT.format(del_len)] = del_calib
        del_llr_range_save_data[
            calibration.DEL_LLR_RNG_TMPLT.format(del_len)] = del_llr_range
    ins_llr_range_save_data, ins_calib_save_data = {}, {}
    for ins_len, (ins_calib, ins_llr_range) in ins_calibs.items():
        ins_calib_save_data[
            calibration.INS_CALIB_TMPLT.format(ins_len)] = ins_calib
        ins_llr_range_save_data[
            calibration.INS_LLR_RNG_TMPLT.format(ins_len)] = ins_llr_range
    np.savez(
        args.out_filename,
        stratify_type=calibration.VAR_CALIB_TYPE,
        smooth_nvals=args.num_calibration_values,
        max_indel_len=max_indel_len,
        **snp_calib_save_data,
        **snp_llr_range_save_data,
        **del_calib_save_data,
        **del_llr_range_save_data,
        **ins_calib_save_data,
        **ins_llr_range_save_data
    )


if __name__ == '__main__':
    _main(get_parser_calibrate_variants().parse_args())
