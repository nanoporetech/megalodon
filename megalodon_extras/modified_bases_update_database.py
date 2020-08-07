import sqlite3
from time import time

from tqdm import tqdm

from megalodon import logging, mods
from ._extras_parsers import get_parser_modified_bases_update_database


DEBUG = False
N_DEBUG = 50000000

INSERT_BATCH_SIZE = 10000

LOGGER = logging.get_logger()


def get_read_id(uuid, read_ids, new_db):
    try:
        read_id = read_ids[uuid]
    except KeyError:
        new_db.cur.execute('INSERT INTO read (uuid) VALUES (?)', (uuid,))
        read_id = new_db.cur.lastrowid
        read_ids[uuid] = read_id
    return read_id, read_ids


def insert_data(new_db, insert_batch):
    new_db.cur.executemany('INSERT INTO data VALUES (?,?,?,?)', insert_batch)


def fill_mods(old_cur, new_db):
    read_ids = {}
    n_recs = old_cur.execute('SELECT MAX(rowid) FROM mods').fetchone()[0]
    old_cur.execute('SELECT * FROM mods')
    insert_batch = []
    for i, (uuid, chrm, strand, pos, score, mod_base, motif, motif_pos,
            raw_motif) in tqdm(enumerate(old_cur), total=n_recs, smoothing=0,
                               dynamic_ncols=True):
        if DEBUG and i > N_DEBUG:
            break
        read_id, read_ids = get_read_id(uuid, read_ids, new_db)
        pos_id = new_db.get_pos_id_or_insert(chrm, strand, pos)
        mod_base_id = new_db.get_mod_base_id_or_insert(
            mod_base, motif, motif_pos, raw_motif)
        insert_batch.append((score, pos_id, mod_base_id, read_id))
        if len(insert_batch) >= INSERT_BATCH_SIZE:
            insert_data(new_db, insert_batch)
            insert_batch = []

    if len(insert_batch) >= 0:
        insert_data(new_db, insert_batch)


def fill_refs(old_cur, new_db):
    old_cur.execute('SELECT DISTINCT chrm FROM mods')
    for ref_name, in old_cur:
        new_db.insert_chrm(ref_name)
    new_db.create_chrm_index()


def _main(args):
    raise NotImplementedError(
        'The previous version of this script updated version 0 to ' +
        'version 1. Updgreade to version 2 not yet implemented.')
    logging.init_logger()
    old_db = sqlite3.connect(args.old_db)
    old_cur = old_db.cursor()
    new_db = mods.ModsDb(args.new_db, read_only=False)

    LOGGER.info('Reading/loading reference record names.')
    fill_refs(old_cur, new_db)

    LOGGER.info('Reading/loading modified base scores.')
    fill_mods(old_cur, new_db)

    if not DEBUG:
        new_db.create_mod_index()
        t0 = time()
        LOGGER.info('Creating positions index.')
        new_db.create_pos_index()
        t1 = time()
        LOGGER.info('Took {} seconds.'.format(t1 - t0))
        LOGGER.info('Creating scores position index.')
        new_db.create_data_covering_index()
        LOGGER.info('Took {} seconds.'.format(time() - t1))
    new_db.close()


if __name__ == '__main__':
    _main(get_parser_modified_bases_update_database().parse_args())
