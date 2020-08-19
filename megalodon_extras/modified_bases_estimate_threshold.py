import numpy as np
from tqdm import tqdm

from megalodon import logging, mods, megalodon_helper as mh
from ._extras_parsers import get_parser_modified_bases_estimate_threshold


LOGGER = logging.get_logger()


def _main(args):
    logging.init_logger()

    LOGGER.info('Loading database position statistics')
    mods_db = mods.ModsDb(
        mh.get_megalodon_fn(args.megalodon_results_dir, mh.PR_MOD_NAME))
    db_mods = set(mod_base for mod_base, _ in mods_db.get_mod_long_names())
    if args.mod_base not in db_mods:
        raise mh.MegaError('Target modified base not found in mods database.')

    scores = []
    bar = tqdm(total=args.num_statistics, smoothing=0)
    for (chrm, strand, pos), mod_llrs in mods_db.iter_pos_scores(
            convert_pos=True, compute_llrs=True):
        for mod_base, reads_llrs in mod_llrs.items():
            if mod_base != args.mod_base:
                continue
            bar.update(len(reads_llrs))
            scores.extend(reads_llrs)
        if args.num_statistics is not None and bar.n >= args.num_statistics:
            break

    LOGGER.info('Esitmating fraction of modified bases')
    scores = np.array(scores)
    frac_mod = args.fraction_modified
    if frac_mod is None:
        thresh_vals = np.percentile(
            scores, (args.mod_percentile, 100 - args.mod_percentile))
        thresh_val = np.abs(thresh_vals).min()
        n_can = np.greater_equal(scores, thresh_val).sum()
        n_mod = np.less_equal(scores, -thresh_val).sum()
        frac_mod = n_mod / (n_mod + n_can)
        print('Fraction mod: {}'.format(frac_mod))
    llr_thresh = np.percentile(scores, frac_mod * 100)
    print('Threshold: {}'.format(llr_thresh))


if __name__ == '__main__':
    _main(get_parser_modified_bases_estimate_threshold().parse_args())
