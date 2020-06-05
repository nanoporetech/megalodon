from collections import defaultdict

import numpy as np
from tqdm import tqdm

from megalodon import mods, megalodon_helper as mh
from ._extras_parsers import get_parser_modified_bases_estimate_threshold


def _main(args):
    mods_db = mods.ModsDb(
        mh.get_megalodon_fn(args.megalodon_results_dir, mh.PR_MOD_NAME))
    scores = []
    num_pos = (mods_db.get_num_uniq_mod_pos()
               if args.num_positions is None else args.num_positions)
    for n_pos, (pos_id, pos_chrm, strand, pos) in tqdm(
            enumerate(mods_db.iter_pos()), total=num_pos, smoothing=0):
        pr_mod_stats = mods_db.get_pos_stats(
            (pos_id, pos_chrm, strand, pos), return_uuids=True)
        mod_type_stats = defaultdict(dict)
        for r_stats in pr_mod_stats:
            mod_type_stats[r_stats.read_id][r_stats.mod_base] = r_stats.score
        for r_mod_stats in mod_type_stats.values():
            mod_lps = np.array(list(r_mod_stats.values()))
            with np.errstate(divide='ignore'):
                can_lp = np.log1p(-np.exp(mod_lps).sum())
            for mod_base, mod_lp in r_mod_stats.items():
                if mod_base != args.mod_base:
                    continue
                scores.append(can_lp - mod_lp)

        if args.num_positions is not None and n_pos >= args.num_positions:
            break

    scores = np.array(scores)
    frac_mod = args.fraction_modified
    if frac_mod is None:
        thresh_vals = np.percentile(scores, (args.mod_percentile,
                                             100 - args.mod_percentile))
        thresh_val = np.abs(thresh_vals).min()
        n_can = np.greater_equal(scores, thresh_val).sum()
        n_mod = np.less_equal(scores, -thresh_val).sum()
        frac_mod = n_mod / (n_mod + n_can)
        print('Fraction mod: {}'.format(frac_mod))
    llr_thresh = np.percentile(scores, frac_mod * 100)
    print('Threshold: {}'.format(llr_thresh))


if __name__ == '__main__':
    _main(get_parser_modified_bases_estimate_threshold().parse_args())
