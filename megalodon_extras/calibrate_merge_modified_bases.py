import sys

import numpy as np

from megalodon import calibration, megalodon_helper as mh
from ._extras_parsers import get_parser_calibrate_merge_modified_bases


def _main(args):
    mh.prep_out_fn(args.out_filename, args.overwrite)

    calib_data = np.load(args.modified_base_calibration_files[-1])
    stratify_type = str(calib_data[calibration.STRAT_TYPE_TXT])
    num_calib_vals = np.int(calib_data[calibration.SMOOTH_NVALS_TXT])
    mod_calibs = {}
    for mod_base in calib_data[calibration.MOD_BASES_TXT]:
        mod_calibs[mod_base] = (
            calib_data[mod_base + calibration.LLR_RANGE_SUFFIX].copy(),
            calib_data[mod_base + calibration.CALIB_TABLE_SUFFIX].copy())
    for mod_calib_fn in args.modified_base_calibration_files[-2::-1]:
        calib_data = np.load(mod_calib_fn)
        assert stratify_type == str(calib_data[calibration.STRAT_TYPE_TXT])
        assert num_calib_vals == np.int(
            calib_data[calibration.SMOOTH_NVALS_TXT])
        for mod_base in calib_data[calibration.MOD_BASES_TXT]:
            # overwrite calibration data with files passed earlier
            mod_calibs[mod_base] = (
                calib_data[mod_base + calibration.LLR_RANGE_SUFFIX].copy(),
                calib_data[mod_base + calibration.CALIB_TABLE_SUFFIX].copy())

    save_kwargs = {}
    for mod_base, (mod_llr_range, mod_calib) in mod_calibs.items():
        save_kwargs[mod_base + calibration.LLR_RANGE_SUFFIX] = mod_llr_range
        save_kwargs[mod_base + calibration.CALIB_TABLE_SUFFIX] = mod_calib

    # save calibration table for reading into mod calibration table
    sys.stderr.write('Saving calibrations to file.\n')
    mod_bases = list(mod_calibs.keys())
    np.savez(
        args.out_filename,
        stratify_type=stratify_type,
        smooth_nvals=num_calib_vals,
        mod_bases=mod_bases,
        **save_kwargs)


if __name__ == '__main__':
    _main(get_parser_calibrate_merge_modified_bases().parse_args())
