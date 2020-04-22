import sys
import argparse
import numpy as np

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


def get_parser():
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
        '--strand-specific-sites', action='store_true',
        help='Sites in --ground-truth-data are strand-specific')
    parser.add_argument(
        '--out-filename', default='mod_calibration_statistics.npz',
        help='Output filename for text summary. Default: %(default)s')
    parser.add_argument(
        '--quiet', action='store_true',
        help='Suppress progress information.')

    return parser


def main():
    args = get_parser().parse_args()
    global VERBOSE
    VERBOSE = not args.quiet

    if args.ground_truth_data is None and \
       args.control_megalodon_results_dir is None:
        sys.stderr.write(
            '***** ERROR ***** Must provide either ' +
            '--control-megalodon-results-dir or --ground-truth-data')
        sys.exit()

    db_fn = mh.get_megalodon_fn(args.megalodon_results_dir,
                                mh.PR_MOD_NAME)
    if args.ground_truth_data is not None:
        if VERBOSE:
            sys.stderr.write('Parsing ground truth data\n')
        gt_mod_pos, gt_can_pos = mh.parse_ground_truth_file(
            args.ground_truth_data, include_strand=args.strand_specific_sites)
        if VERBOSE:
            sys.stderr.write((
                'Loaded ground truth data with {} modified sites and {} ' +
                'canonical sites.\n').format(len(gt_mod_pos), len(gt_can_pos)))
            sys.stderr.write(
                'Reading ground truth modified base statistics from ' +
                'database\n')
        all_mod_llrs, all_can_llrs = mods.extract_stats_at_valid_sites(
            db_fn, [gt_mod_pos, gt_can_pos],
            include_strand=args.strand_specific_sites)
    else:
        if VERBOSE:
            sys.stderr.write(
                'Reading ground truth modified base statistics from ' +
                'database\n')
        all_mod_llrs = mods.extract_all_stats(db_fn)
        if VERBOSE:
            sys.stderr.write(
                'Reading ground truth modified base statistics from ' +
                'canonical sample database\n')
            sys.stderr.write('Reading ground truth \n')
        all_can_llrs = mods.extract_all_stats(mh.get_megalodon_fn(
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
