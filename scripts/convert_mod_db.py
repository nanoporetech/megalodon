import sys
import sqlite3
import argparse
from tqdm import tqdm
from time import time

from megalodon import logging, megalodon_helper as mh, mods

DEBUG = False
N_DEBUG = 10000000

INSERT_BATCH_SIZE = 100000

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'old_db',
        help='Megalodon version 0.1 modified base data base.')
    parser.add_argument(
        '--new-db', default='megalodon_mods.db',
        help='Output data base name. Should replace ' +
        'per_read_modified_base_calls.db in megalodon results directory in ' +
        'order to process further. Default: %(default)s')

    return parser

def get_read_id(uuid, read_ids, new_db):
    try:
        read_id = read_ids[uuid]
    except KeyError:
        new_db.cur.execute('INSERT INTO read (uuid) VALUES (?)', (uuid,))
        read_id = new_db.cur.lastrowid
        read_ids[uuid] = read_id
    return read_id, read_ids

def get_pos_id(chrm, strand, pos, poss, new_db, chrms):
    chrm_id = chrms[chrm]
    try:
        pos_id = poss[(chrm_id, strand, pos)]
    except KeyError:
        new_db.cur.execute('INSERT INTO pos (pos_chrm, pos, strand) ' +
                           'VALUES (?, ?, ?)', (chrm_id, pos, strand))
        pos_id = new_db.cur.lastrowid
        poss[(chrm_id, strand, pos)] = pos_id
    return pos_id, poss

def insert_data(new_db, insert_batch):
    new_db.cur.executemany('INSERT INTO data VALUES (?,?,?,?)', insert_batch)
    return

def fill_mods(old_cur, new_db, chrms):
    read_ids, poss = {}, {}
    n_recs = old_cur.execute('SELECT MAX(rowid) FROM mods').fetchone()[0]
    old_cur.execute('SELECT * FROM mods')
    insert_batch = []
    for i, (uuid, chrm, strand, pos, score, mod_base, motif, motif_pos,
         raw_motif) in tqdm(enumerate(old_cur), total=n_recs, smoothing=0,
                            dynamic_ncols=True):
        if DEBUG and i > N_DEBUG: break
        read_id, read_ids = get_read_id(uuid, read_ids, new_db)
        pos_id, poss = get_pos_id(chrm, strand, pos, poss, new_db, chrms)
        mod_base_id = new_db.get_mod_base_id_or_insert(
            mod_base, motif, motif_pos, raw_motif)
        insert_batch.append((score, pos_id, mod_base_id, read_id))
        if len(insert_batch) >= INSERT_BATCH_SIZE:
            insert_data(new_db, insert_batch)
            insert_batch = []

    if len(insert_batch) >= 0:
        insert_data(new_db, insert_batch)

    return

def fill_refs(old_cur, new_db):
    old_cur.execute('SELECT DISTINCT chrm FROM mods')
    chrms = {}
    for ref_name, in old_cur:
        chrms[ref_name] = new_db.add_chrm(ref_name)
    return chrms

def main():
    args = get_parser().parse_args()

    old_db = sqlite3.connect(args.old_db)
    old_cur = old_db.cursor()
    new_db = mods.ModsDb(args.new_db, read_only=False)

    sys.stderr.write('Reading/loading reference record names.\n')
    chrms = fill_refs(old_cur, new_db)

    sys.stderr.write('Reading/loading modified base scores.\n')
    fill_mods(old_cur, new_db, chrms)

    if not DEBUG:
        t0 = time()
        sys.stderr.write('Creating positions index.\n')
        new_db.create_pos_index()
        t1 = time()
        sys.stderr.write('Took {} seconds.\n'.format(t1 - t0))
        sys.stderr.write('Creating scores position index.\n')
        new_db.create_data_pos_index()
        sys.stderr.write('Took {} seconds.\n'.format(time() - t1))
    new_db.close()

    return

if __name__ == '__main__':
    main()
