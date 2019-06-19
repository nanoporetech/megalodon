import os
import sys
import argparse
import numpy as np
import pandas as pd

from megalodon import megalodon_helper as mh


VERBOSE = False


def output_mods_data(mod_dat, ctrl_dat, gt_dat, mod_chrm_sw, out_fn):
    if VERBOSE: sys.stderr.write('Merging modified base data\n')
    # merge data with known mod sites
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

    m_dat['llr'] = m_dat['can_log_prob'] - m_dat['mod_log_prob']
    with open(out_fn, 'w') as fp:
        for _, pos_dat in m_dat.iterrows():
            fp.write('{}\t{}\t{}\n'.format(
                pos_dat.is_mod, pos_dat.llr, pos_dat.mod_base))

    return

def parse_mod_data(args):
    if VERBOSE: sys.stderr.write('Reading megalodon data\n')
    try:
        mod_dat = pd.read_csv(
            mh.get_megalodon_fn(args.megalodon_results_dir,
                                mh.PR_MOD_TXT_NAME), sep='\t')
    except FileNotFoundError:
        sys.write('ERROR: Must provide a valid Megalodon result directory.')
        sys.exit(1)

    return mod_dat

def parse_control_mods(args):
    ctrl_dat = gt_dat = mod_chrm_sw = None
    if args.control_megalodon_results_dir is not None:
        if VERBOSE: sys.stderr.write('Reading control mods data\n')
        try:
            ctrl_dat = pd.read_csv(
                mh.get_megalodon_fn(args.megalodon_results_dir,
                                    mh.PR_MOD_TXT_NAME), sep='\t')
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
        sys.stderr.write('ERROR: Must provide a control data type.\n')
        sys.exit(1)

    return ctrl_dat, gt_dat, mod_chrm_sw


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
    '--out-filename', default='mod_calibration_statistics.txt',
    help='Output filename for text summary. Default: %(default)s')
parser.add_argument(
    '--quiet', action='store_true',
    help='Suppress progress information.')


def main():
    args = parser.parse_args()
    global VERBOSE
    VERBOSE = not args.quiet

    mod_dat = parse_mod_data(args)
    ctrl_dat, gt_dat, mod_chrm_sw = parse_control_mods(args)
    output_mods_data(mod_dat, ctrl_dat, gt_dat, mod_chrm_sw, args.out_filename)

    return

if __name__ == '__main__':
    main()
