import sys
import numpy as np

from megalodon import logging, megalodon_helper as mh, mods
from ._extras_parsers import get_parser_calibrate_generate_modified_bases_stats


LOGGER = logging.get_logger()


def output_mods_data(
        all_mod_llrs, all_can_llrs, mod_base_set, exclude_mod_bases, out_fn):
    LOGGER.info('Merging modified base data')
    all_mod_bases = list(set(all_mod_llrs.keys()).intersection(
        all_can_llrs.keys()))
    if len(set(all_mod_llrs.keys()).difference(all_mod_bases)) > 0:
        LOGGER.warning(
            'Modified base(s) found in modified dataset which were not ' +
            'found in canonical dataset: {}'.format(','.join(
                set(all_mod_llrs.keys()).difference(all_mod_bases))))
    if len(set(all_can_llrs.keys()).difference(all_mod_bases)) > 0:
        LOGGER.warning(
            'Modified base(s) found in modified dataset which were ' +
            'not found in canonical dataset: {}'.format(','.join(
                set(all_mod_llrs.keys()).difference(all_mod_bases))))
    if mod_base_set is not None:
        all_mod_bases = list(set(all_mod_bases).intersection(mod_base_set))
        if len(all_mod_bases) == 0:
            LOGGER.error((
                'No modified bases to process.\n\tModified bases from ' +
                'results: {}\n\tModified base set: {}').format(
                    ','.join(all_mod_bases), ','.join(mod_base_set)))
    if exclude_mod_bases is not None:
        all_mod_bases = list(set(all_mod_bases).difference(exclude_mod_bases))
        if len(all_mod_bases) == 0:
            LOGGER.error((
                'No modified bases to process.\n\tModified bases from ' +
                'results: {}\n\tExcluded modified bases: {}').format(
                    ','.join(all_mod_bases), ','.join(exclude_mod_bases)))
    mod_base_stats = {mods.GT_ALL_MOD_BASE_STR: all_mod_bases}
    for mod_base in all_mod_bases:
        mod_base_stats[mods.GT_MOD_LLR_STR.format(
            mod_base)] = all_mod_llrs[mod_base]
        mod_base_stats[mods.GT_CAN_LLR_STR.format(
            mod_base)] = all_can_llrs[mod_base]
    np.savez(out_fn, **mod_base_stats)


def _main(args):
    logging.init_logger(quiet=args.quiet)

    if args.ground_truth_data is None and \
       args.control_megalodon_results_dir is None:
        LOGGER.error(
            'Must provide either --control-megalodon-results-dir or ' +
            '--ground-truth-data')
        sys.exit()

    db_fn = mh.get_megalodon_fn(args.megalodon_results_dir,
                                mh.PR_MOD_NAME)
    if args.ground_truth_data is not None:
        LOGGER.info('Parsing ground truth data')
        gt_mod_pos, gt_can_pos = mh.parse_ground_truth_file(
            args.ground_truth_data, include_strand=args.strand_specific_sites)
        LOGGER.info((
            'Loaded ground truth data with {} modified sites and {} ' +
            'canonical sites.').format(len(gt_mod_pos), len(gt_can_pos)))
        LOGGER.info(
            'Reading ground truth modified base statistics from ' +
            'database.')
        all_mod_llrs, all_can_llrs = mods.extract_stats_at_valid_sites(
            db_fn, [gt_mod_pos, gt_can_pos],
            include_strand=args.strand_specific_sites)
    else:
        LOGGER.info(
            'Reading ground truth modified base statistics from ' +
            'database')
        all_mod_llrs = mods.extract_all_stats(db_fn)
        LOGGER.info(
            'Reading ground truth modified base statistics from ' +
            'canonical sample database')
        all_can_llrs = mods.extract_all_stats(mh.get_megalodon_fn(
            args.control_megalodon_results_dir, mh.PR_MOD_NAME))

    mod_summary = [
        (mod, len(all_mod_llrs[mod]) if mod in all_mod_llrs else 0,
         len(all_can_llrs[mod]) if mod in all_can_llrs else 0)
        for mod in set(all_mod_llrs).union(all_can_llrs)]
    LOGGER.info(
        'Data summary:\n\tmod\tmod_N\tcan_N\n' + '\n'.join(
            '\t' + '\t'.join(map(str, x)) for x in mod_summary))
    output_mods_data(
        all_mod_llrs, all_can_llrs, args.modified_bases_set,
        args.exclude_modified_bases, args.out_filename)


if __name__ == '__main__':
    _main(get_parser_calibrate_generate_modified_bases_stats().parse_args())
