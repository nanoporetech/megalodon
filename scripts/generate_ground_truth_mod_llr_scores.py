import sys
import argparse
import numpy as np
from collections import defaultdict

from megalodon import megalodon_helper as mh, mods


VERBOSE = False


def output_mods_data(all_mod_llrs, all_can_llrs, out_fn):
    if VERBOSE:
        sys.stderr.write('Merging modified base data\n')
    all_mod_bases = list(set(all_mod_llrs.keys()).intersection(
        all_can_llrs.keys()))
    if len(set(all_mod_llrs.keys()).difference(all_mod_bases)) > 0:
        sys.stderr.write(
            'WARNING: Modified base(s) found in modified dataset which were ' +
            'not found in canonical dataset: {}'.format(','.join(
                set(all_mod_llrs.keys()).difference(all_mod_bases))))
    if len(set(all_can_llrs.keys()).difference(all_mod_bases)) > 0:
        sys.stderr.write(
            'WARNING: Modified base(s) found in modified dataset which were ' +
            'not found in canonical dataset: {}'.format(','.join(
                set(all_mod_llrs.keys()).difference(all_mod_bases))))
    mod_base_stats = {mods.GT_ALL_MOD_BASE_STR: all_mod_bases}
    for mod_base in all_mod_bases:
        mod_base_stats[mods.GT_MOD_LLR_STR.format(
            mod_base)] = all_mod_llrs[mod_base]
        mod_base_stats[mods.GT_CAN_LLR_STR.format(
            mod_base)] = all_can_llrs[mod_base]
    np.savez(out_fn, **mod_base_stats)


def get_samp_stats(mods_db_fn):
    all_llrs = defaultdict(list)
    mods_db = mods.ModsDb(mods_db_fn)
    for mods_pos_llrs, pos in mods_db.iter_pos_scores_from_covering_index():
        for mod_base, mod_pos_llrs in mods_pos_llrs.items():
            all_llrs[mod_base].append(mod_pos_llrs)
    all_llrs = dict((mod_base, np.concatenate(mod_llrs))
                    for mod_base, mod_llrs in all_llrs.items())
    return all_llrs


def get_samp_stats_w_ground_truth(mods_db_fn, gt_mod_pos, gt_can_pos):
    all_mod_llrs = defaultdict(list)
    all_can_llrs = defaultdict(list)
    mods_db = mods.ModsDb(mods_db_fn)
    for mods_pos_llrs, (
            chrm, _, pos) in mods_db.iter_pos_scores_from_covering_index(
                return_pos=True):
        if (chrm, pos) in gt_mod_pos:
            for mod_base, mod_pos_llrs in mods_pos_llrs.items():
                all_mod_llrs[mod_base].append(mod_pos_llrs)
        elif (chrm, pos) in gt_can_pos:
            for mod_base, mod_pos_llrs in mods_pos_llrs.items():
                all_can_llrs[mod_base].append(mod_pos_llrs)
    all_mod_llrs = dict((mod_base, np.concatenate(mod_llrs))
                        for mod_base, mod_llrs in all_mod_llrs.items())
    all_can_llrs = dict((mod_base, np.concatenate(mod_llrs))
                        for mod_base, mod_llrs in all_can_llrs.items())
    return all_mod_llrs, all_can_llrs


def parse_ground_truth_data(gt_data_fn):
    gt_mod_pos = set()
    gt_can_pos = set()
    with open(gt_data_fn) as fp:
        for line in fp:
            chrm, pos, is_mod = line.strip().split(',')
            if bool(is_mod):
                gt_mod_pos.add((chrm, int(pos)))
            else:
                gt_can_pos.add((chrm, int(pos)))
    return gt_mod_pos, gt_can_pos


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
    '--out-filename', default='mod_calibration_statistics.npz',
    help='Output filename for text summary. Default: %(default)s')
parser.add_argument(
    '--quiet', action='store_true',
    help='Suppress progress information.')


def main():
    args = parser.parse_args()
    global VERBOSE
    VERBOSE = not args.quiet

    db_fn = mh.get_megalodon_fn(args.megalodon_results_dir,
                                mh.PR_MOD_TXT_NAME)
    if args.ground_truth_data is not None:
        gt_mod_pos, gt_can_pos = parse_ground_truth_data(
            args.ground_truth_data)
        all_mod_llrs, all_can_llrs = get_samp_stats_w_ground_truth(
            db_fn, gt_mod_pos, gt_can_pos)
    else:
        all_mod_llrs = get_samp_stats(db_fn)
        all_can_llrs = get_samp_stats(mh.get_megalodon_fn(
            args.control_megalodon_results_dir, mh.PR_MOD_TXT_NAME))
    output_mods_data(all_mod_llrs, all_can_llrs, args.out_filename)


if __name__ == '__main__':
    main()
