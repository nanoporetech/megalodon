from collections import defaultdict

import numpy as np

from megalodon import logging, mods
from ._extras_parsers import get_parser_calibrate_merge_modified_bases_stats


LOGGER = logging.get_logger()


def _main(args):
    logging.init_logger()

    fn_mod_base_llrs = defaultdict(lambda: ([], []))
    for llr_fn in args.modified_base_calibration_stats_files:
        llrs_data = np.load(llr_fn)
        for mod_base in llrs_data[mods.GT_ALL_MOD_BASE_STR]:
            fn_mod_base_llrs[mod_base][0].append(
                llrs_data[mods.GT_MOD_LLR_STR.format(mod_base)]
            )
            fn_mod_base_llrs[mod_base][1].append(
                llrs_data[mods.GT_CAN_LLR_STR.format(mod_base)]
            )

    mod_base_stats = {mods.GT_ALL_MOD_BASE_STR: list(fn_mod_base_llrs)}
    for mod_base, (mod_llrs, can_llrs) in fn_mod_base_llrs.items():
        mod_base_stats[mods.GT_MOD_LLR_STR.format(mod_base)] = np.concatenate(
            mod_llrs
        )
        mod_base_stats[mods.GT_CAN_LLR_STR.format(mod_base)] = np.concatenate(
            can_llrs
        )
    np.savez(args.out_filename, **mod_base_stats)


if __name__ == "__main__":
    _main(get_parser_calibrate_merge_modified_bases_stats().parse_args())
