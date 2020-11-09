import numpy as np
from tqdm import tqdm

from megalodon import mods, megalodon_helper as mh
from ._extras_parsers import get_parser_per_read_text_modified_bases


def _main(args):
    mods_db = mods.ModsDb(
        mh.get_megalodon_fn(args.megalodon_results_dir, mh.PR_MOD_NAME),
        in_mem_dbid_to_uuid=True)
    mods_txt_fp = open(
        mh.get_megalodon_fn(args.megalodon_results_dir, mh.PR_MOD_TXT_NAME)
        if args.out_filename is None else args.out_filename, 'w')
    mods_txt_fp.write('\t'.join(mods_db.text_field_names) + '\n')
    rec_tmplt = '\t'.join('{}' for _ in mods_db.text_field_names) + '\n'
    bar = tqdm(desc='Processing Per-read Data', unit='per-read calls',
               total=mods_db.get_num_uniq_stats(), smoothing=0,
               dynamic_ncols=True)
    for (chrm, strand, pos), pos_lps in mods_db.iter_pos_scores(
            convert_pos=True):
        bar.update(len(pos_lps))
        str_strand = mh.int_strand_to_str(strand)
        mod_out_text = ''
        prev_dbid = None
        mod_bs, r_lps = [], []
        for read_dbid, mod_dbid, lp in sorted(pos_lps):
            if prev_dbid != read_dbid and prev_dbid is not None:
                uuid = mods_db.get_uuid(prev_dbid)
                # compute and store log likelihood ratios
                with np.errstate(divide='ignore'):
                    can_lp = np.log1p(-np.exp(r_lps).sum())
                for mod_b, r_lp in zip(mod_bs, r_lps):
                    mod_out_text += rec_tmplt.format(
                        uuid, chrm, str_strand, pos, r_lp, can_lp, mod_b)
                mod_bs, r_lps = [], []
            prev_dbid = read_dbid
            mod_bs.append(mods_db.get_mod_base(mod_dbid))
            r_lps.append(lp)
        uuid = mods_db.get_uuid(prev_dbid)
        # compute and store log likelihood ratios
        with np.errstate(divide='ignore'):
            can_lp = np.log1p(-np.exp(r_lps).sum())
        for mod_b, r_lp in zip(mod_bs, r_lps):
            mod_out_text += rec_tmplt.format(
                uuid, chrm, str_strand, pos, r_lp, can_lp, mod_b)
        mods_txt_fp.write(mod_out_text)
    mods_txt_fp.close()


if __name__ == '__main__':
    _main(get_parser_per_read_text_modified_bases().parse_args())
