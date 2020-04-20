import sys
import argparse
import numpy as np
from collections import defaultdict

from megalodon import megalodon_helper as mh, mods


VERBOSE = False
TRUE_TEXT_VALUES = set(('y', 'yes', 't', 'true', 'on', '1'))
FALSE_TEXT_VALUES = set(('n', 'no', 'f', 'false', 'off', '0'))


def output_mods_data(all_mod_llrs, all_can_llrs, out_fn):
    if VERBOSE:
        sys.stderr.write('Merging modified base data\n')
    all_mod_bases = list(set(all_mod_llrs.keys()).intersection(
        all_can_llrs.keys()))
    if len(set(all_mod_llrs.keys()).difference(all_mod_bases)) > 0:
        sys.stderr.write(
            'WARNING: Modified base(s) found in modified dataset which were ' +
            'not found in canonical dataset: {}\n'.format(','.join(
                set(all_mod_llrs.keys()).difference(all_mod_bases))))
    if len(set(all_can_llrs.keys()).difference(all_mod_bases)) > 0:
        sys.stderr.write(
            'WARNING: Modified base(s) found in modified dataset which were ' +
            'not found in canonical dataset: {}\n'.format(','.join(
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
    for pos, mods_pos_llrs in mods_db.iter_pos_scores():
        for mod_base, mod_pos_llrs in mods_pos_llrs.items():
            print(mod_bases)
            all_llrs[mod_base].append(mod_pos_llrs)
    all_llrs = dict((mod_base, np.concatenate(mod_llrs))
                    for mod_base, mod_llrs in all_llrs.items())
    return all_llrs


def get_samp_stats_w_ground_truth(mods_db_fn, gt_mod_pos, gt_can_pos):
    all_mod_llrs = defaultdict(list)
    all_can_llrs = defaultdict(list)
    mods_db = mods.ModsDb(mods_db_fn, pos_index_in_memory=False)
    for (chrm, _, pos), mods_pos_llrs in mods_db.iter_pos_scores(
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


def text_to_bool(val):
    lower_val = val.lower()
    if lower_val in TRUE_TEXT_VALUES:
        return True
    elif lower_val in FALSE_TEXT_VALUES:
        return False
    raise mh.MegaError('Invalid boolean string encountered: "{}".'.format(val))


def parse_ground_truth_data(gt_data_fn):
    gt_mod_pos = set()
    gt_can_pos = set()
    with open(gt_data_fn) as fp:
        for line in fp:
            chrm, pos, is_mod = line.strip().split(',')
            if text_to_bool(is_mod):
                gt_mod_pos.add((chrm, int(pos)))
            else:
                gt_can_pos.add((chrm, int(pos)))
    if VERBOSE:
        sys.stderr.write((
            'Loaded ground truth data with {} modified sites and {} ' +
            'canonical sites.\n').format(len(gt_mod_pos) , len(gt_can_pos)))
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
                                mh.PR_MOD_NAME)
    if args.ground_truth_data is not None:
        if VERBOSE:
            sys.stderr.write('Parsing ground truth data\n')
        gt_mod_pos, gt_can_pos = parse_ground_truth_data(
            args.ground_truth_data)
        if VERBOSE:
            sys.stderr.write('Reading sample data\n')
        all_mod_llrs, all_can_llrs = get_samp_stats_w_ground_truth(
            db_fn, gt_mod_pos, gt_can_pos)
    else:
        if VERBOSE:
            sys.stderr.write('Reading ground truth modified sample\n')
        all_mod_llrs = get_samp_stats(db_fn)
        if VERBOSE:
            sys.stderr.write('Reading ground truth canonical sample\n')
        all_can_llrs = get_samp_stats(mh.get_megalodon_fn(
            args.control_megalodon_results_dir, mh.PR_MOD_NAME))
    if VERBOSE:
        mod_summary = [
            (mod, len(all_mod_llrs[mod]) if mod in all_mod_llrs else 0,
             len(all_can_llrs[mod]) if mod in all_can_llrs else 0)
            for mod in set(all_mod_llrs).union(all_can_llrs)]
        sys.stderr.write(
            'Data summary:\n\tmod\tmod_N\tcan_N\n' + '\n'.join(
                '\t' + '\t'.join(map(str, x)) for x in mod_summary) + '\n')
    output_mods_data(all_mod_llrs, all_can_llrs, args.out_filename)


if __name__ == '__main__':
    main()
