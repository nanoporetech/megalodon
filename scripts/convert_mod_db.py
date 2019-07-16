import sys
import sqlite3
import argparse
from tqdm import tqdm

from megalodon import logging, megalodon_helper as mh, mods

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

def fill_mods(old_cur, new_db):
    read_ids = {}
    n_recs = old_cur.execute('SELECT COUNT(*) FROM mods').fetchone()[0]
    old_cur.execute('SELECT * FROM mods')
    for (uuid, chrm, strand, pos, score, mod_base, motif, motif_pos,
         raw_motif) in tqdm(old_cur, total=n_recs, smoothing=0):
        try:
            read_id = read_ids[uuid]
        except KeyError:
            new_db.cur.execute('INSERT INTO read (uuid) VALUES (?)', (uuid,))
            read_id = new_db.cur.lastrowid
            read_ids[uuid] = read_id
        pos_id = new_db.get_pos_id_or_insert(chrm, strand, pos)
        mod_base_id = new_db.get_mod_base_id_or_insert(
            mod_base, motif, motif_pos, raw_motif)
        new_db.cur.execute('INSERT INTO data VALUES (?,?,?,?)',
                           (score, pos_id, mod_base_id, read_id))
    return

def fill_refs(old_cur, new_db):
    old_cur.execute('SELECT DISTINCT chrm FROM mods')
    for ref_name, in old_cur:
        new_db.add_chrm(ref_name)
    return

def main():
    args = get_parser().parse_args()

    old_db = sqlite3.connect(args.old_db)
    old_cur = old_db.cursor()
    new_db = mods.ModsDb(args.new_db, read_only=False)

    sys.stderr.write('Reading/loading reference record names.\n')
    fill_refs(old_cur, new_db)

    sys.stderr.write('Reading/loading modified base scores.\n')
    fill_mods(old_cur, new_db)

    new_db.create_data_pos_index()
    new_db.close()

    return

if __name__ == '__main__':
    main()
