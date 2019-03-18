import sys
import argparse
from collections import defaultdict

import numpy as np
from tqdm import tqdm
from scipy import stats

DO_PLOT = False
SMOOTH_BW = 0.01

def extract_and_filter_llhrs(args):
    all_llhrs = defaultdict(dict)
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
        llhrs = list(map(float, fields['SAMPLE']['LLHRS'].split(',')))
        all_llhrs[(fields['CHROM'], fields['POS'])][
            (fields['REF'], fields['ALT'])] = llhrs

    # create filter to remove all SNPs near homozygous alt calls
    expand_homo_alt_filter = set()
    for chrom, pos in homo_alt_calls:
        for filt_pos in range(pos - args.homozygous_alt_mask_width,
                              pos + args.homozygous_alt_mask_width + 1):
            expand_homo_alt_filter.add((chrom, filt_pos))

    return np.array([
        llhr
        for chrm_pos, alts in all_llhrs.items()
        for ref_alt, alt_llhrs in alts.items()
        for llhr in alt_llhrs
        if chrm_pos not in expand_homo_alt_filter])

def compute_calibration(filt_llhrs, args):
    sys.stderr.write('Computing emperical density.\n')
    kern = stats.gaussian_kde(filt_llhrs, bw_method=SMOOTH_BW)
    smooth_ls = np.linspace(-args.max_input_llhr, args.max_input_llhr,
                            args.num_calibration_values, endpoint=True)
    smooth_vals = kern(smooth_ls)

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


SMOOTH_MAX, SMOOTH_NVALS = 200, 1001
def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--snps-vcf', default='megalodon_results/snps.vcf',
        help='VCF file produced from megalodon (with VCF_OUTPUT_LLHRS=True). ' +
        'Default: %(default)s')
    parser.add_argument(
        '--coverage-threshold', type=int, default=10,
        help='Minimum coverage to include a site. Default: %(default)d')
    parser.add_argument(
        '--homozygous-alt-mask-width', type=int, default=10,
        help='Window to mask around homozygous alt calls. Default: %(default)d')
    parser.add_argument(
        '--max-input-llhr', type=int, default=200,
        help='Maximum log-likelihood ratio to compute calibration. ' +
        'Default: %(default)d')
    parser.add_argument(
        '--num-calibration-values', type=int, default=1001,
        help='Number of discrete calibration values to compute. ' +
        'Default: %(default)d')
    parser.add_argument(
        '--out-filename', default='megalodon_snp_calibration.npz',
        help='Filename to output calibration values. Default: %(default)s')

    return parser

def main():
    args = get_parser().parse_args()

    filt_llhrs = extract_and_filter_llhrs(args)

    mono_llhr = compute_calibration(filt_llhrs, args)

    # save valibration table for reading into SNP table
    np.savez(
        args.out_filename,
        stratify_type='none',
        smooth_max=args.max_input_llhr,
        smooth_nvals=args.num_calibration_values,
        global_calibration_table=mono_llhr)

    return

if __name__ == '__main__':
    main()
