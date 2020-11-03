from collections import namedtuple

import pysam
from tqdm import tqdm

from megalodon import backends, logging, mods, megalodon_helper as mh
from ._extras_parsers import get_parser_modified_bases_split_calls_by_motif


LOGGER = logging.get_logger()

MOTIF_INFO = namedtuple('MOTIF_INFO', (
    'bases_before', 'bases_after', 'raw_motif', 'motif', 'rc_motif'))


########################
# data table functions #
########################

def split_data(in_mods_db, out_mods_dbs, ref):
    LOGGER.info('Inserting modified base data')
    bar = tqdm(desc='Inserting Data', unit='per-read calls',
               total=in_mods_db.get_num_uniq_stats(), smoothing=0,
               dynamic_ncols=True)
    curr_ref_seq = curr_chrm = None
    # TODO multiprocess over contigs (need to implement iteration over range
    # of pos_dbids via chrm string)
    for pos_dbid, pos_mod_data in in_mods_db.iter_pos_scores():
        bar.update(len(pos_mod_data))
        chrm, strand, pos = in_mods_db.get_pos(pos_dbid)
        if chrm != curr_chrm:
            curr_chrm = chrm
            curr_ref_seq = ref.fetch(chrm)
        for out_mods_db, motif_info in out_mods_dbs:
            motif_match = (
                motif_info.motif.match(
                    curr_ref_seq[pos - motif_info.bases_before:
                                 pos + motif_info.bases_after + 1])
                if strand == 1 else
                motif_info.rc_motif.match(
                    curr_ref_seq[pos - motif_info.bases_after:
                                 pos + motif_info.bases_before + 1]))
            if motif_match is not None:
                pos_insert_data = [(lp, pos_dbid, mod_dbid, read_dbid)
                                   for read_dbid, mod_dbid, lp in pos_mod_data]
                out_mods_db.insert_batch_data(pos_insert_data)
                break
    bar.close()


##########
# motifs #
##########

def parse_motifs(raw_motifs):
    motifs = []
    for raw_motif, bases_before in raw_motifs:
        bases_before = int(bases_before)
        bases_after = len(raw_motif) - bases_before - 1
        motif = mh.compile_motif_pat(raw_motif)
        rc_motif = mh.compile_rev_comp_motif_pat(raw_motif)
        motifs.append(MOTIF_INFO(
            bases_before=bases_before, bases_after=bases_after,
            raw_motif=raw_motif, motif=motif, rc_motif=rc_motif))

    return motifs


########
# main #
########

def _main(args):
    logging.init_logger(
        args.megalodon_directory, out_suffix=args.output_suffix)

    # parse motifs
    motifs = parse_motifs(args.motif)
    # open indexed FASTA reference
    ref = pysam.FastaFile(args.reference)

    LOGGER.info('Extracting mods and chrms from input database')
    in_mods_db = mods.ModsDb(mh.get_megalodon_fn(
        args.megalodon_directory, mh.PR_MOD_NAME))
    alphabet, _, mod_long_names = in_mods_db.get_alphabet_info()
    ref_names_and_lens = list(zip(*in_mods_db.iter_chrms()))[1:]
    LOGGER.info('Extracting read uuid table')
    in_uuids = [uuid for _, uuid in in_mods_db.iter_uuids()]

    LOGGER.info('Opening new per-read modified base statistics databases')
    model_info = backends.DetachedModelInfo(
        alphabet=alphabet, mod_long_names=mod_long_names)
    out_mods_dbs = []
    for motif_info in motifs:
        out_dir = '{}.{}_{}'.format(args.output_prefix, motif_info.raw_motif,
                                    motif_info.bases_before)
        mh.mkdir(out_dir, overwrite=False)
        mods_info = mods.ModInfo(model_info, out_dir=out_dir)
        mods.init_mods_db(mods_info, ref_names_and_lens)
        out_mods_dbs.append((
            mods.ModsDb(mods_info.mods_db_fn, read_only=False),
            motif_info))
        out_mods_dbs[-1][0].insert_uuids(in_uuids)
        out_mods_dbs[-1][0].commit()

    # commit so read uuids are available to worker processes
    LOGGER.info('Inserting per-read calls from input databases')
    split_data(in_mods_db, out_mods_dbs, ref)

    # TOOD do this in separate processes
    LOGGER.info(
        'Creating data covering indices for efficient iteration by position')
    for out_mods_db, _ in out_mods_dbs:
        out_mods_db.create_data_covering_index()
        out_mods_db.commit()
        out_mods_db.close()
        LOGGER.info('Finished indexing {}'.format(out_mods_db.fn))


if __name__ == '__main__':
    _main(get_parser_modified_bases_split_calls_by_motif().parse_args())
