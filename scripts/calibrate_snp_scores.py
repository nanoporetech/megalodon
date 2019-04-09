import sys
import argparse
from collections import defaultdict

import numpy as np
from tqdm import tqdm


DO_PLOT = False
SMOOTH_BW = 0.8
SMOOTH_MAX = 200
SMOOTH_NVALS = 1001


def extract_and_filter_llrs(args):
    all_llrs = defaultdict(dict)
    homo_alt_calls = set()
    sys.stderr.write('Extracting log-likelihood ratios.\n')
    for line in tqdm(open(args.snps_vcf), smoothing=0):
        if line.startswith('##'): continue
        if line.startswith('#'):
            header = line[1:].split()
            continue
        fields = dict(zip(header, line.split()))
        fields['INFO'] = dict(f.split('=') for f in fields['INFO'].split(':'))
        if int(fields['INFO']['DP']) < args.coverage_threshold: continue

        fields['SAMPLE'] = dict((f, s) for f, s in zip(
            fields['FORMAT'].split(':'), fields['SAMPLE'].split(':')))
        if fields['SAMPLE']['GT'] == '1/1':
            homo_alt_calls.add((fields['CHROM'], int(fields['POS'])))
        llrs = list(map(float, fields['SAMPLE']['LLRS'].split(',')))
        all_llrs[(fields['CHROM'], fields['POS'])][
            (fields['REF'], fields['ALT'])] = llrs

    # create filter to remove all SNPs near homozygous alt calls
    expand_homo_alt_filter = set()
    for chrom, pos in homo_alt_calls:
        for filt_pos in range(pos - args.homozygous_alt_mask_width,
                              pos + args.homozygous_alt_mask_width + 1):
            expand_homo_alt_filter.add((chrom, filt_pos))

    return np.array([
        llr
        for chrm_pos, alts in all_llrs.items()
        for ref_alt, alt_llrs in alts.items()
        for llr in alt_llrs
        if chrm_pos not in expand_homo_alt_filter])

def compute_calibration(filt_llrs, args):
    def guassian(x):
        return (np.exp(-x ** 2 / (2 * SMOOTH_BW ** 2)) /
                (SMOOTH_BW * np.sqrt(2 * np.pi)))


    sys.stderr.write('Computing emperical density.\n')
    smooth_ls = np.linspace(-args.max_input_llr, args.max_input_llr,
                            args.num_calibration_values, endpoint=True)

    smooth_vals = np.zeros(smooth_ls.shape[0])
    for llr in tqdm(filt_llrs, smoothing=0):
        smooth_vals += guassian(smooth_ls - llr)
    smooth_vals /= filt_llrs.shape[0]

    sys.stderr.write('Computing emperical likelihood.\n')
    peak_site = np.argmax(smooth_vals)
    # force monotonic increasing before peak and monotonic decreasing after
    mono_smooth_vals = np.concatenate([
        np.mean(np.stack([
            np.maximum.accumulate(smooth_vals[:peak_site]),
            np.minimum.accumulate(smooth_vals[:peak_site][::-1])[::-1]]),
                axis=0),
        np.mean(np.stack([
            np.minimum.accumulate(smooth_vals[peak_site:]),
            np.maximum.accumulate(smooth_vals[peak_site:][::-1])[::-1]]),
                axis=0)])
    prob_alt = (mono_smooth_vals[::-1] /
                (mono_smooth_vals + mono_smooth_vals[::-1]))
    # need to compute this for non-symmetric distributions
    prob_mp = args.num_calibration_values // 2
    # force monotonic decreasing with reverse maximum before p=0.5 and
    # forward minimum after p=0.5
    mono_prob = np.concatenate([
        np.maximum.accumulate(prob_alt[:prob_mp][::-1])[::-1],
        np.minimum.accumulate(prob_alt[prob_mp:])])

    if DO_PLOT:
        import matplotlib
        if sys.platform == 'darwin':
            matplotlib.use("TkAgg")
        import matplotlib.pyplot as plt
        #plt.ion()

        f, axarr = plt.subplots(3, sharex=True)
        axarr[0].plot(smooth_ls, mono_smooth_vals, color='red')
        axarr[0].plot(smooth_ls, mono_smooth_vals[::-1])
        axarr[1].plot(smooth_ls, mono_prob, color='orange')
        axarr[1].plot(smooth_ls, 1 / (np.exp(smooth_ls) + 1), color='purple')
        axarr[2].plot(smooth_ls, np.log((1 - prob_alt) / prob_alt), color='red')
        axarr[2].plot(smooth_ls, np.log((1 - mono_prob) / mono_prob),
                      color='orange')
        plt.show()

    return np.log((1 - mono_prob) / mono_prob)


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--snps-vcf', default='megalodon_results/snps.vcf',
        help='VCF file produced from megalodon (with --write-vcf-llr). ' +
        'Default: %(default)s')
    parser.add_argument(
        '--coverage-threshold', type=int, default=10,
        help='Minimum coverage to include a site. Default: %(default)d')
    parser.add_argument(
        '--homozygous-alt-mask-width', type=int, default=10,
        help='Window to mask around homozygous alt calls. Default: %(default)d')
    parser.add_argument(
        '--max-input-llr', type=int, default=SMOOTH_MAX,
        help='Maximum log-likelihood ratio to compute calibration. ' +
        'Default: %(default)d')
    parser.add_argument(
        '--num-calibration-values', type=int, default=SMOOTH_NVALS,
        help='Number of discrete calibration values to compute. ' +
        'Default: %(default)d')
    parser.add_argument(
        '--out-filename', default='megalodon_snp_calibration.npz',
        help='Filename to output calibration values. Default: %(default)s')

    return parser

def main():
    args = get_parser().parse_args()

    filt_llrs = extract_and_filter_llrs(args)

    mono_llr = compute_calibration(filt_llrs, args)

    # save valibration table for reading into SNP table
    np.savez(
        args.out_filename,
        stratify_type='none',
        smooth_max=args.max_input_llr,
        smooth_nvals=args.num_calibration_values,
        global_calibration_table=mono_llr)

    return

if __name__ == '__main__':
    main()
