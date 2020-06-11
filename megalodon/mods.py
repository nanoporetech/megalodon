import os
import re
import sys
import pysam
import queue
import sqlite3
import datetime
from array import array
from multiprocessing.connection import wait
from collections import defaultdict, namedtuple, OrderedDict

import numpy as np

from megalodon import (
    calibration, decode, logging, mapping, megalodon_helper as mh)
from megalodon._version import MEGALODON_VERSION


AGG_INFO = namedtuple('AGG_INFO', ('method', 'binary_threshold'))
DEFAULT_AGG_INFO = AGG_INFO(mh.MOD_BIN_THRESH_NAME, None)

ANNOT_MODS = namedtuple('ANNOT_MODS', ('mod_seq', 'mod_qual'))

# slots in ground truth mods numpy arrays for calibration
GT_ALL_MOD_BASE_STR = 'all_mod_bases'
GT_MOD_LLR_STR = '{}_mod_llrs'
GT_CAN_LLR_STR = '{}_can_llrs'

FIXED_VCF_MI = [
    'INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
    'INFO=<ID=SN,Number=1,Type=String,Description="Strand">',
    'FORMAT=<ID=VALID_DP,Number=1,Type=Integer,' +
    'Description="Valid Read Depth">',
]
MOD_MI_TMPLTS = [
    'FORMAT=<ID={0},Number=1,Type=Float,Description=' +
    '"{1} Modified Base Proportion">']
FORMAT_LOG_PROB_MI = (
    'FORMAT=<ID=LOG_PROBS,Number=A,Type=String,' +
    'Description="Per-read log10 likelihoods for modified ' +
    'bases (semi-colon separated)">')

SAMPLE_NAME = 'SAMPLE'
MOD_MAP_RG_ID = '1'
MOD_MAP_MAX_QUAL = 40

POS_IDX_CHNG_ERR_MSG = (
    'Inserting chromosomes forces change in in-memory index size. Please ' +
    'insert all chromosomes at database initialization.')

OUT_BUFFER_LIMIT = 10000

_PROFILE_MODS_QUEUE = False
_PROFILE_MODS_AUX = False

LOGGER = logging.get_logger()


###########
# Mods DB #
###########

class ModsDb(object):
    # note foreign key constraint is not applied here as this
    # drastically reduces efficiency. Namely the search for pos_id
    # when inserting into the data table forces a scan of a very
    # large table or maintainance of a very large pos table index
    # both of which slow data base speed
    # thus foreign key constraint must be handled by the class
    db_tables = OrderedDict((
        ('chrm', OrderedDict((
            ('chrm_id', 'INTEGER PRIMARY KEY'),
            ('chrm', 'TEXT'),
            ('chrm_len', 'INTEGER')))),
        ('pos', OrderedDict((
            ('pos_id', 'INTEGER PRIMARY KEY'),
            ('pos_chrm', 'INTEGER'),
            ('strand', 'INTEGER'),
            ('pos', 'INTEGER')))),
        ('mod', OrderedDict((
            ('mod_id', 'INTEGER PRIMARY KEY'),
            ('mod_base', 'TEXT'),
            ('motif', 'TEXT'),
            ('motif_pos', 'INTEGER'),
            ('raw_motif', 'TEXT')))),
        ('mod_long_names', OrderedDict((
            ('mod_base', 'TEXT'),
            ('mod_long_name', 'TEXT')))),
        ('read', OrderedDict((
            ('read_id', 'INTEGER PRIMARY KEY'),
            ('uuid', 'TEXT')))),
        ('data', OrderedDict((
            ('score', 'FLOAT'),
            ('score_pos', 'INTEGER'),
            ('score_mod', 'INTEGER'),
            ('score_read', 'INTEGER'))))))

    pos_mem_dt = np.uint32
    pos_mem_max = np.iinfo(pos_mem_dt).max

    # namedtuple for returning mods from a single position
    mod_data = namedtuple('mod_data', [
        'read_id', 'chrm', 'strand', 'pos', 'score', 'mod_base', 'motif',
        'motif_pos', 'raw_motif'])
    text_field_names = (
        'read_id', 'chrm', 'strand', 'pos', 'mod_log_prob',
        'can_log_prob', 'mod_base', 'motif')

    def __init__(self, fn, init_db_tables=False, read_only=True, db_safety=1,
                 in_mem_chrm_to_dbid=False, in_mem_dbid_to_chrm=False,
                 in_mem_pos_to_dbid=False, in_mem_dbid_to_pos=False,
                 in_mem_mod_to_dbid=False, in_mem_dbid_to_mod=False,
                 in_mem_uuid_to_dbid=False, in_mem_dbid_to_uuid=False,
                 force_uint32_pos_to_dbid=False):
        """ Interface to database containing modified base statistics.

        Default settings are for optimal read_only performance.
        """
        self.fn = mh.resolve_path(fn)
        self.init_db_tables = init_db_tables
        self.read_only = read_only
        self.in_mem_chrm_to_dbid = in_mem_chrm_to_dbid or in_mem_pos_to_dbid
        self.in_mem_pos_to_dbid = in_mem_pos_to_dbid
        self.in_mem_mod_to_dbid = in_mem_mod_to_dbid
        self.in_mem_uuid_to_dbid = in_mem_uuid_to_dbid
        self.in_mem_dbid_to_chrm = in_mem_dbid_to_chrm or in_mem_dbid_to_pos
        self.in_mem_dbid_to_pos = in_mem_dbid_to_pos
        self.in_mem_dbid_to_mod = in_mem_dbid_to_mod
        self.in_mem_dbid_to_uuid = in_mem_dbid_to_uuid
        self.force_uint32 = force_uint32_pos_to_dbid

        if self.read_only:
            if not os.path.exists(self.fn):
                LOGGER.error((
                    'Modified base per-read database file ({}) does ' +
                    'not exist.').format(self.fn))
                raise mh.MegaError('Invalid mods DB filename.')
            self.db = sqlite3.connect('file:' + self.fn + '?mode=ro', uri=True)
        else:
            self.db = sqlite3.connect(self.fn, timeout=mh.SQLITE_TIMEOUT)

        # initialize main cursor
        self.cur = self.db.cursor()

        if self.init_db_tables:
            # create tables
            for tbl_name, tbl in self.db_tables.items():
                try:
                    self.cur.execute("CREATE TABLE {} ({})".format(
                        tbl_name, ','.join((
                            '{} {}'.format(*ft) for ft in tbl.items()))))
                except sqlite3.OperationalError:
                    raise mh.MegaError(
                        'Modified bases database already exists. Either ' +
                        'provide location for new database or open in ' +
                        'read_only mode.')
        else:
            self.check_tables_init()

        if self.read_only:
            # use memory mapped file access
            self.cur.execute('PRAGMA mmap_size = {}'.format(
                mh.MEMORY_MAP_LIMIT))
        else:
            if db_safety < 2:
                # set asynchronous mode to off for max speed
                self.cur.execute('PRAGMA synchronous = OFF')
            if db_safety < 1:
                # set no rollback mode
                self.cur.execute('PRAGMA journal_mode = OFF')

        # load requested in memory indices
        self.load_in_mem_chrm()
        self.load_in_mem_pos()
        self.load_in_mem_mod()
        self.load_in_mem_uuid()

    def check_tables_init(self):
        missing_tables = []
        for tbl_name in self.db_tables:
            if len(self.cur.execute(
                    'SELECT name FROM sqlite_master ' +
                    'WHERE type="table" AND name=?',
                    (tbl_name, )).fetchall()) == 0:
                missing_tables.append(tbl_name)
        if len(missing_tables) > 0:
            raise mh.MegaError(
                'Per-read modified base database tables not initialized. '
                'Missing tables: {}'.format(', '.join(missing_tables)))

    def check_data_covering_index_exists(self):
        if len(self.cur.execute(
                'SELECT name FROM sqlite_master WHERE type="index" AND name=?',
                ('data_cov_idx', )).fetchall()) == 0:
            raise mh.MegaError('Data covering index does not exist.')

    #########################
    # getter data functions #
    #########################
    def get_chrm_dbid(self, chrm):
        """ Get database ID.

        Args:
            chrm (str): Chromosome name

        Returns:
            Database ID (int).
        """
        try:
            if self.in_mem_chrm_to_dbid:
                dbid = self.chrm_to_dbid[chrm]
            else:
                dbid = self.cur.execute(
                    'SELECT chrm_id FROM chrm WHERE chrm=?',
                    (chrm,)).fetchone()[0]
        except (TypeError, KeyError):
            raise mh.MegaError('chrm "{}" not found in database.'.format(chrm))
        return dbid

    def get_chrm(self, chrm_dbid):
        """ Get chromosome name and length

        Args:
            chrm_dbid (int): Chromosome database ID

        Returns:
            Chromosome name (str) and length (int)
        """
        try:
            if self.in_mem_dbid_to_chrm:
                return self.dbid_to_chrm[chrm_dbid]
            return self.cur.execute(
                'SELECT chrm, chrm_len FROM chrm WHERE chrm_id=?',
                (chrm_dbid, )).fetchall()[0]
        except (TypeError, KeyError):
            raise mh.MegaError('Reference record (chromosome ID) not found ' +
                               'in mods database.')

    def get_mod_base_data(self, mod_dbid):
        """ Get modified base information from database

        Args:
            mod_dbid (int): Modified base database ID

        Returns:
            Modified base single letter code (str), motif searched pattern
                (str), relative modified base position in motif (int; 0-based),
                raw motif observed from reference (str)
        """
        try:
            if self.in_mem_dbid_to_mod:
                mod_base_data = self.dbid_to_mod[mod_dbid]
            else:
                mod_base_data = self.cur.execute(
                    'SELECT mod_base, motif, motif_pos, raw_motif FROM mod ' +
                    'WHERE mod_id=?', (mod_dbid,)).fetchall()[0]
        except (TypeError, KeyError):
            raise mh.MegaError('Modified base not found in mods database.')
        return mod_base_data

    def get_uuid(self, read_dbid):
        """ Get read UUID from database.

        Args:
            read_dbid (int): Database read UUID ID

        Returns:
            UUID (int): Universal read identifier
        """
        try:
            if self.in_mem:
                uuid = self.dbid_to_uuid[read_dbid]
            else:
                uuid = self.cur.execute(
                    'SELECT uuid FROM read WHERE read_id=?',
                    (read_dbid, )).fetchall()[0][0]
        except (TypeError, KeyError):
            raise mh.MegaError('Read ID not found in vars database.')
        return uuid

    def get_pos(self, pos_dbid):
        """ Get position information from database

        Args:
            pos_dbid (int): Position database ID

        Returns:
            Chromosome name (str), strand (int; +1=forward strand, -1=reverse
                strand), position (int; 0-based)
        """
        try:
            if self.in_mem_dbid_to_pos:
                chrm_dbid, strand, pos = self.dbid_to_pos[pos_dbid]
            else:
                chrm_dbid, strand, pos = self.cur.execute(
                    'SELECT pos_chrm, strand, pos FROM pos WHERE pos_id=?',
                    (pos_dbid,)).fetchone()
        except (TypeError, KeyError):
            raise mh.MegaError('Position not found in database.')
        chrm, _ = self.get_chrm(chrm_dbid)
        return chrm, strand, pos

    def get_pos_stats(self, pos_data, return_uuids=False):
        """ Get all statistics mapped to a reference position. The data
        covering index should be created to increase the access speed for
        this function.

        Args:
            pos_data (tuple): Positon database ID (int), chromosome database id
                (int), strand (int; +1=forward strand, -1=reverse strand),
                position 0-based (int)
            return_uuids (bool): Whether to return database read ids (defalt)
                or UUIDs.

        Returns:
            List containing megalodon.mods.ModsDb.mod_data objects mapped to
                specified reference position.
        """
        read_id_conv = self.get_uuid if return_uuids else lambda x: x
        # these attributes are specified in self.iter_pos
        pos_dbid, chrm_dbid, strand, pos = pos_data
        chrm, _ = self.get_chrm(chrm_dbid)
        return [
            self.mod_data(read_id_conv(read_dbid), chrm, strand, pos, score,
                          *self.get_mod_base_data(mod_dbid))
            for read_dbid, score, mod_dbid in self.cur.execute(
                    'SELECT score_read, score, score_mod FROM data ' +
                    'WHERE score_pos=?', (pos_dbid, )).fetchall()]

    def get_num_uniq_chrms(self):
        """ Get number of chromosomes/contigs stored in the database
        """
        num_chrms = self.cur.execute(
            'SELECT MAX(chrm_id) FROM chrm').fetchone()[0]
        if num_chrms is None:
            return 0
        return num_chrms

    def get_num_uniq_mods(self):
        """ Get the number of unique modified base database entries. Including
        searched and observed motif.
        """
        num_mod_bases = self.cur.execute(
            'SELECT MAX(mod_id) FROM mod').fetchone()[0]
        if num_mod_bases is None:
            return 0
        return num_mod_bases

    def get_num_uniq_reads(self):
        """ Get number of unique reads in database
        """
        num_reads = self.cur.execute(
            'SELECT MAX(read_id) FROM read').fetchone()[0]
        if num_reads is None:
            return 0
        return num_reads

    def get_num_uniq_mod_pos(self):
        """ Get number of unique positions in database
        """
        num_pos = self.cur.execute('SELECT MAX(pos_id) FROM pos').fetchone()[0]
        if num_pos is None:
            return 0
        return num_pos

    def get_num_uniq_stats(self):
        """ Get number of per-read modified base statistics stored in database
        """
        num_stats = self.cur.execute(
            'SELECT MAX(rowid) FROM data').fetchone()[0]
        if num_stats is None:
            return 0
        return num_stats

    ##################################
    # load in-memory index functions #
    ##################################
    def load_in_mem_chrm(self):
        if not (self.in_mem_chrm_to_dbid or self.in_mem_dbid_to_chrm):
            return
        dbid_chrm = self.cur.execute(
            'SELECT chrm_id, chrm, chrm_len FROM chrm').fetchall()
        if self.in_mem_chrm_to_dbid:
            self.chrm_to_dbid = dict(
                (chrm, dbid) for dbid, chrm, _ in dbid_chrm)
        if self.in_mem_dbid_to_chrm:
            self.dbid_to_chrm = dict(
                (dbid, (chrm, chrm_len)) for dbid, chrm, chrm_len in dbid_chrm)

    def check_in_mem_pos_size(self, chrm_lens):
        if 2 * sum(chrm_lens) > self.pos_mem_max:
            if self.force_uint32:
                LOGGER.warning('Forcing uint32')
                return
            if len(self.pos_to_dbid) > 0:
                raise mh.MegaError(POS_IDX_CHNG_ERR_MSG)
            self.pos_mem_dt = np.uint64
            self.pos_mem_max = np.iinfo(self.pos_mem_dt).max
            if 2 * sum(chrm_lens) > self.pos_mem_max:
                raise mh.MegaError(
                    'Cannot store modified base positions in memory for ' +
                    'a genome this size. Please specify ' +
                    '--mod-positions-on-disk.')

    def load_in_mem_pos(self):
        if not (self.in_mem_pos_to_dbid or self.in_mem_dbid_to_pos):
            return
        dbid_pos = self.cur.execute(
            'SELECT pos_id, pos_chrm, strand, pos FROM pos').fetchall()
        if self.in_mem_pos_to_dbid:
            self.pos_to_dbid = {}
            self.check_in_mem_pos_size([
                chrm_len for _, _, chrm_len in self.iter_chrms()])
            # position to database id is stored as a dictionary of numpy arrays
            # this was done since standard dictionaries cause memory errors
            # for large genomes potentially toward the end of a long run.
            self.pos_to_dbid = dict(
                ((chrm_dbid, strand),
                 np.full(chrm_len, self.pos_mem_max, self.pos_mem_dt))
                for chrm_dbid, _, chrm_len in self.iter_chrms()
                for strand in (1, -1))
            for dbid, chrm_dbid, strand, pos in dbid_pos:
                self.pos_to_dbid[(chrm_dbid, strand)][pos] = dbid
        if self.in_mem_dbid_to_pos:
            self.dbid_to_pos = dict(
                (dbid, (chrm_dbid, strand, pos))
                for dbid, chrm_dbid, strand, pos in dbid_pos)

    def load_in_mem_mod(self):
        if not (self.in_mem_mod_to_dbid or self.in_mem_dbid_to_mod):
            return
        dbid_mod = self.cur.execute(
            'SELECT mod_id, mod_base, motif, motif_pos, raw_motif ' +
            'FROM mod').fetchall()
        if self.in_mem_mod_to_dbid:
            self.mod_to_dbid = dict(
                ((mod_base, motif, motif_pos, raw_motif), dbid)
                for dbid, mod_base, motif, motif_pos, raw_motif in dbid_mod)
        if self.in_mem_dbid_to_mod:
            self.dbid_to_mod = dict(
                (dbid, (mod_base, motif, motif_pos, raw_motif))
                for dbid, mod_base, motif, motif_pos, raw_motif in dbid_mod)

    def load_in_mem_uuid(self):
        if not (self.in_mem_uuid_to_dbid or self.in_mem_dbid_to_uuid):
            return
        dbid_uuid = self.cur.execute(
            'SELECT read_id, uuid FROM read').fetchall()
        if self.in_mem_uuid_to_dbid:
            self.uuid_to_dbid = dict((uuid, dbid) for dbid, uuid in dbid_uuid)
        if self.in_mem_dbid_to_uuid:
            self.dbid_to_uuid = dict(dbid_uuid)

    #########################
    # insert data functions #
    #########################
    def get_chrm_dbid_or_insert(self, chrm, chrm_len):
        try:
            chrm_dbid = self.get_chrm_dbid(chrm)
        except mh.MegaError:
            self.cur.execute('INSERT INTO chrm (chrm, chrm_len) VALUES (?,?)',
                             (chrm, chrm_len))
            chrm_dbid = self.cur.lastrowid
            if self.in_mem_chrm_to_dbid:
                self.chrm_to_dbid[chrm] = chrm_dbid
            if self.in_mem_dbid_to_chrm:
                self.dbid_to_chrm[chrm_dbid] = (chrm, chrm_len)
            if self.in_mem_pos_to_dbid:
                if sum(chrm_len for _, _, chrm_len in
                       self.iter_chrms()) * 2 > self.pos_mem_max:
                    raise mh.MegaError(POS_IDX_CHNG_ERR_MSG)
                self.pos_to_dbid[(chrm_dbid, 1)] = np.full(
                    chrm_len, self.pos_mem_max, self.pos_mem_dt)
                self.pos_to_dbid[(chrm_dbid, -1)] = np.full(
                    chrm_len, self.pos_mem_max, self.pos_mem_dt)
        return chrm_dbid

    def insert_chrms(self, ref_names_and_lens):
        next_chrm_id = self.get_num_uniq_chrms() + 1
        self.cur.executemany('INSERT INTO chrm (chrm, chrm_len) VALUES (?,?)',
                             zip(*ref_names_and_lens))
        chrm_dbids = list(range(next_chrm_id,
                                next_chrm_id + len(ref_names_and_lens[0])))
        if self.in_mem_chrm_to_dbid:
            self.chrm_to_dbid.update(zip(
                ref_names_and_lens[0], chrm_dbids))
        if self.in_mem_dbid_to_chrm:
            self.dbid_to_chrm.update(zip(chrm_dbids, zip(*ref_names_and_lens)))

        # Add position index numpy arrays for new chromosomes if stored in mem
        if self.in_mem_pos_to_dbid:
            self.check_in_mem_pos_size(ref_names_and_lens[1])
            self.pos_to_dbid.update(
                ((chrm_i + next_chrm_id, strand), np.full(
                    chrm_len, self.pos_mem_max, self.pos_mem_dt))
                for chrm_i, chrm_len in enumerate(ref_names_and_lens[1])
                for strand in (1, -1))

    def get_pos_dbid_or_insert(self, chrm_dbid, strand, pos):
        try:
            if self.in_mem_pos_to_dbid:
                pos_dbid = int(self.pos_to_dbid[(chrm_dbid, strand)][pos])
                if pos_dbid == self.pos_mem_max:
                    raise TypeError
            else:
                pos_dbid = self.cur.execute(
                    'SELECT pos_id FROM pos WHERE pos_chrm=? AND strand=? ' +
                    'AND pos=?', (chrm_dbid, strand, pos)).fetchall()[0][0]
        except TypeError:
            self.cur.execute(
                'INSERT INTO pos (pos_chrm, strand, pos) VALUES (?,?,?)',
                (chrm_dbid, strand, pos))
            pos_dbid = self.cur.lastrowid
            if self.in_mem_pos_to_dbid:
                self.pos_to_dbid[(chrm_dbid, strand)][pos] = pos_dbid
            if self.in_mem_dbid_to_pos:
                self.dbid_to_pos[pos_dbid] = (chrm_dbid, strand, pos)
        return pos_dbid

    def get_pos_dbids_or_insert(self, r_uniq_pos, chrm_dbid, strand):
        """ Get position database IDs. If positions are not found in the
        database they will be inserted.

        Args:
            r_uniq_pos (list): Unique positions (int).
            chrm_dbid (int): Chromosome database ID.
            strand (int): +1 = forward strandl -1 = reverse strand

        Returns:
            List of position database IDs for each entry in r_uniq_pos.
        """
        def add_pos_to_db(pos_to_add):
            # separate function for profiling
            self.cur.executemany(
                'INSERT INTO pos (pos_chrm, strand, pos) VALUES (?,?,?)',
                ((chrm_dbid, strand, int(pos)) for pos in pos_to_add))

        if len(r_uniq_pos) == 0:
            return []

        if self.in_mem_pos_to_dbid:
            cs_pos_to_dbid = self.pos_to_dbid[(chrm_dbid, strand)]
            # extract positions that have not yet been observed
            pos_to_add = [
                r_uniq_pos[to_add_idx] for to_add_idx in np.where(np.equal(
                    cs_pos_to_dbid[r_uniq_pos], self.pos_mem_max))[0]]
        else:
            cs_pos_to_dbid = dict(
                pos_and_dbid for pos in r_uniq_pos
                for pos_and_dbid in self.cur.execute(
                        'SELECT pos, pos_id FROM pos ' +
                        'WHERE pos_chrm=? AND strand=? AND pos=?',
                        (chrm_dbid, strand, pos)).fetchall())
            pos_to_add = tuple(r_uniq_pos.difference(cs_pos_to_dbid))

        # insert positions and update in memory indices
        if len(pos_to_add) > 0:
            next_pos_id = self.get_num_uniq_mod_pos() + 1
            add_pos_to_db(pos_to_add)
            pos_dbids = list(range(next_pos_id, next_pos_id + len(pos_to_add)))
            if self.in_mem_pos_to_dbid:
                cs_pos_to_dbid[pos_to_add] = np.array(
                    pos_dbids, dtype=self.pos_mem_dt)
            else:
                cs_pos_to_dbid.update(zip(pos_to_add, pos_dbids))
            if self.in_mem_dbid_to_pos:
                self.dbid_to_pos.update(
                    (pos_dbid, (chrm_dbid, strand, pos))
                    for pos_dbid, pos in zip(pos_dbids, pos_to_add))

        # return pos_id for each entry in r_uniq_pos
        return [int(cs_pos_to_dbid[pos]) for pos in r_uniq_pos]

    def insert_mod_long_names(self, mod_long_names):
        self.cur.executemany(
            'INSERT INTO mod_long_names (mod_base, mod_long_name) ' +
            'VALUES (?,?)', mod_long_names)

    def get_mod_base_dbid_or_insert(
            self, mod_base, motif, motif_pos, raw_motif):
        try:
            if self.in_mem_mod_to_dbid:
                mod_dbid = self.mod_to_dbid[
                    (mod_base, motif, motif_pos, raw_motif)]
            else:
                mod_dbid = self.cur.execute(
                    'SELECT mod_id FROM mod WHERE mod_base=? AND motif=? ' +
                    'AND motif_pos=? AND raw_motif=?',
                    (mod_base, motif, motif_pos, raw_motif)).fetchall()[0][0]
        except (TypeError, KeyError):
            self.cur.execute(
                'INSERT INTO mod (mod_base, motif, motif_pos, raw_motif) ' +
                'VALUES (?,?,?,?)', (mod_base, motif, motif_pos, raw_motif))
            mod_dbid = self.cur.lastrowid
            if self.in_mem_mod_to_dbid:
                self.mod_to_dbid[
                    (mod_base, motif, motif_pos, raw_motif)] = mod_dbid
            if self.in_mem_dbid_to_mod:
                self.dbid_to_mod[mod_dbid] = (
                    mod_base, motif, motif_pos, raw_motif)
        return mod_dbid

    def get_mod_base_ids_or_insert(self, r_uniq_mod_bases):
        """ Get modified base database IDs. If modified bases are not found in
        the database they will be inserted.

        Args:
            r_uniq_mod_bases (list): List of modified base descriptions. Each
                entry must be composed of 1) modified base single letter code
                (str), 2) motif searched pattern (str), 3) relative modified
                base position in motif (int; 0-based), and 4) raw motif
                observed from reference (str)

        Returns:
            List with modified base ID for each entry in r_uniq_mod_bases
        """
        if len(r_uniq_mod_bases) == 0:
            return []

        if self.in_mem_mod_to_dbid:
            mod_to_dbid = self.mod_to_dbid
        else:
            mod_to_dbid = dict(
                (mod_data, mod_dbid[0])
                for mod_data in r_uniq_mod_bases for mod_dbid in
                self.cur.execute(
                    'SELECT mod_id FROM mod WHERE mod_base=? AND motif=? ' +
                    'AND motif_pos=? AND raw_motif=?', mod_data).fetchall())

        mod_bases_to_add = tuple(set(r_uniq_mod_bases).difference(mod_to_dbid))

        if len(mod_bases_to_add) > 0:
            next_mod_base_id = self.get_num_uniq_mods() + 1
            self.cur.executemany(
                'INSERT INTO mod (mod_base, motif, motif_pos, raw_motif) ' +
                'VALUES (?,?,?,?)', mod_bases_to_add)
            mod_dbids = list(range(next_mod_base_id,
                                   next_mod_base_id + len(mod_bases_to_add)))
            # update either extracted entries from DB or in memory index
            mod_to_dbid.update(zip(mod_bases_to_add, mod_dbids))
            if self.in_mem_dbid_to_mod:
                self.dbid_to_mod.update(zip(mod_dbids, mod_bases_to_add))

        return [mod_to_dbid[mod_data] for mod_data in r_uniq_mod_bases]

    def get_read_dbid_or_insert(self, uuid):
        """ Get database ID for a read uuid. If value is not found in the
        database it will be inserted.

        Args:
            uuid (str): Unique read identifier

        Returns:
            Database ID (int)
        """
        try:
            if self.in_mem_uuid_to_dbid:
                read_dbid = self.uuid_to_dbid[uuid]
            else:
                read_dbid = self.cur.execute(
                    'SELECT read_id FROM read WHERE uuid=?',
                    (uuid, )).fetchone()[0]
        except (TypeError, KeyError):
            self.cur.execute('INSERT INTO read (uuid) VALUES (?)', (uuid, ))
            read_dbid = self.cur.lastrowid
            if self.in_mem_dbid_to_uuid:
                self.uuid_to_dbid[uuid] = read_dbid
        return read_dbid

    def get_read_dbids_or_insert(self, uuids):
        """ Get database IDs for a list of uuids. If values are not found in
        the database they will be inserted.

        Args:
            uuids (list): Unique read identifiers (str)

        Returns:
            List of database IDs (int)
        """
        if self.in_mem_uuid_to_dbid:
            uuid_to_dbid = self.uuid_to_dbid
        else:
            uuid_to_dbid = dict(
                (uuid, uuid_dbid[0]) for uuid in uuids for uuid_dbid in
                self.cur.execute('SELECT read_id FROM read WHERE uuid=?',
                                 uuid).fetchall())

        uuids_to_add = tuple(set(uuids).difference(uuid_to_dbid))

        if len(uuids_to_add) > 0:
            next_read_dbid = self.get_num_uniq_reads() + 1
            self.cur.executemany(
                'INSERT INTO read (uuid) VALUES (?)',
                ((uuid, ) for uuid in uuids_to_add))
            read_dbids = list(range(next_read_dbid,
                                    next_read_dbid + len(uuids_to_add)))
            # update either extracted entries from DB or in memory index
            uuid_to_dbid.update(zip(uuids_to_add, read_dbids))
            if self.in_mem_dbid_to_uuid:
                self.dbid_to_uuid.update(zip(read_dbids, uuids_to_add))

        return [uuid_to_dbid[uuid] for uuid in uuids]

    def insert_read_uuid(self, uuid):
        """ Insert unique read identifier into database.

        Args:
            uuid (str): Unique read identifier

        Returns:
            Database ID.
        """
        self.cur.execute('INSERT INTO read (uuid) VALUES (?)', (uuid,))
        return self.cur.lastrowid

    def insert_read_data(self, r_insert_data):
        """ Insert data from a single read.

        Args:
            r_insert_data (list): Each entry should contain 4 elements:
                1) modified base log probability (float), 2) position database
                ID 3) modified base database ID, and 4) read database ID.
        """
        self.cur.executemany('INSERT INTO data VALUES (?,?,?,?)',
                             r_insert_data)

    def insert_data(self, score, pos_id, mod_base_id, read_id):
        self.cur.execute(
            'INSERT INTO data (score, score_pos, score_mod, score_read) ' +
            'VALUES (?,?,?,?)', (score, pos_id, mod_base_id, read_id))
        return self.cur.lastrowid

    #########################
    # table index functions #
    #########################
    def create_chrm_index(self):
        self.cur.execute('CREATE UNIQUE INDEX chrm_idx ON chrm(chrm)')

    def create_pos_index(self):
        self.cur.execute('CREATE UNIQUE INDEX pos_idx ON pos' +
                         '(pos_chrm, strand, pos)')

    def create_mod_index(self):
        self.cur.execute('CREATE UNIQUE INDEX mod_idx ON ' +
                         'mod(mod_base, motif, motif_pos, raw_motif)')

    def create_data_covering_index(self):
        # TODO add progress info to this step.
        self.cur.execute('CREATE INDEX data_cov_idx ON data(' +
                         'score_pos, score_read, score_mod, score)')

    ##################
    # data iterators #
    ##################
    def iter_pos_scores(self, return_pos=False):
        """ Iterate log likelihood ratios (log(P_can / P_mod)) by position.
        Yield tuples of position (either id or converted genomic coordiante)
        and dictionary with mod base single letter keys pointing to list of
        log likelihood ratios.

        Note this function iterates over the index created by
        create_data_covering_index, so should be very fast.

        If return_pos is set to True, position data will be returned in format
        (chrm, strand, pos); else database position id will be returned
        """
        def extract_pos_llrs(pos_lps):
            mod_llrs = defaultdict(list)
            for read_pos_lps in pos_lps.values():
                mt_lps = np.array([lp for _, lp in read_pos_lps])
                with np.errstate(divide='ignore'):
                    can_lp = np.log1p(-np.exp(mt_lps).sum())
                for mod_dbid, lp in read_pos_lps:
                    mod_llrs[mod_dbid].append(can_lp - lp)
            return dict(mod_llrs)

        self.check_data_covering_index_exists()
        # set function to transform pos_id
        extract_pos = self.get_pos if return_pos else lambda x: x
        # use local cursor since extracting pos or mod might use class cursor
        local_cursor = self.db.cursor()
        local_cursor.execute(
            'SELECT score_pos, score_mod, score_read, score FROM data ' +
            'ORDER BY score_pos')
        # initialize variables with first value
        prev_pos, mod_dbid, read_id, lp = local_cursor.fetchone()
        pos_lps = defaultdict(list)
        for curr_pos, mod_dbid, read_id, lp in local_cursor:
            if curr_pos != prev_pos:
                yield extract_pos(prev_pos), extract_pos_llrs(pos_lps)
                pos_lps = defaultdict(list)
            pos_lps[read_id].append((self.get_mod_base_data(mod_dbid)[0], lp))
            prev_pos = curr_pos
        yield extract_pos(prev_pos), extract_pos_llrs(pos_lps)

    def get_all_chrm_and_lens(self):
        """ Get chromosome names and lengths

        Returns:
            Tuple with two lists. First is chromosome names (str) and second
            is chromosome lengths (int).
        """
        try:
            self.cur.execute('SELECT chrm, chrm_len FROM chrm')
            return tuple(zip(*self.cur))
        except sqlite3.OperationalError:
            raise mh.MegaError(
                'Old megalodon database scheme detected. Please re-run ' +
                'megalodon processing or downgrade megalodon installation.')

    def iter_chrms(self):
        """ Iterate over chromosomes from database

        Yields:
            Tuple containing 1) database ID (int), 2) chromosome name
                (str), and 3) length (int)
        """
        if self.in_mem_dbid_to_chrm:
            for chrm_dbid, (chrm, chrm_len) in self.dbid_to_chrm.items():
                yield chrm_dbid, chrm, chrm_len
        else:
            # use local cursor since other processing might use class cursor
            local_cursor = self.db.cursor()
            local_cursor.execute('SELECT chrm_id, chrm, chrm_len FROM chrm')
            for chrm_dbid, chrm, chrm_len in local_cursor:
                yield chrm_dbid, chrm, chrm_len

    def iter_mod_bases(self):
        """ Iterate over modified base information from database

        Yields:
            Tuple containing 1) database ID, 2) modified base single letter
                code (str), 3) motif searched pattern (str), 4) relative
                modified base position in motif (int; 0-based), 5) raw motif
                observed from reference (str)
        """
        if self.in_mem_dbid_to_mod:
            for mod_dbid, (mod_base, motif, motif_pos,
                           raw_motif) in self.dbid_to_mod.items():
                yield mod_dbid, mod_base, motif, motif_pos, raw_motif
        elif self.in_mem_mod_to_dbid:
            for (mod_base, motif, motif_pos,
                 raw_motif), mod_dbid in self.mod_to_dbid.items():
                yield mod_dbid, mod_base, motif, motif_pos, raw_motif
        else:
            # use local cursor since other processing might use class cursor
            local_cursor = self.db.cursor()
            local_cursor.execute(
                'SELECT mod_id, mod_base, motif, motif_pos, raw_motif ' +
                'FROM mod')
            for (mod_dbid, mod_base, motif, motif_pos,
                 raw_motif) in local_cursor:
                yield mod_dbid, mod_base, motif, motif_pos, raw_motif

    def get_mod_long_names(self):
        """ Get modified base long names

        Returns:
            List of tuples containing 1) modified base single letter code and
                2) modified base long name
        """
        return self.cur.execute(
            'SELECT mod_base, mod_long_name FROM mod_long_names').fetchall()

    def iter_uuids(self):
        """ Iterate over UUIDs from database

        Yields:
            Tuple containing 1) database ID, 2) UUID
        """
        if self.in_mem_dbid_to_uuid:
            for read_dbid, uuid in self.dbid_to_uuid.items():
                yield read_dbid, uuid
        elif self.in_mem_dbid_to_uuid:
            for uuid, read_dbid in self.dbid_to_uuid.items():
                yield read_dbid, uuid
        else:
            # use local cursor since other processing might use class cursor
            local_cursor = self.db.cursor()
            local_cursor.execute('SELECT read_id, uuid FROM read')
            for read_dbid, uuid in local_cursor:
                yield read_dbid, uuid

    def iter_pos(self):
        if self.in_mem_dbid_to_pos:
            for pos_dbid, (chrm_dbid, strand, pos) in self.dbid_to_pos.items():
                yield pos_dbid, chrm_dbid, strand, pos
        elif self.in_mem_pos_to_dbid:
            for (chrm_dbid, strand), cs_pos in self.pos_to_dbid.items():
                valid_cs_pos = np.where(np.not_equal(
                    cs_pos, self.pos_mem_max))[0]
                for pos, pos_dbid in zip(valid_cs_pos, cs_pos[valid_cs_pos]):
                    yield pos_dbid, chrm_dbid, strand, pos
        else:
            # use local cursor since other processing might use class cursor
            local_cursor = self.db.cursor()
            local_cursor.execute(
                'SELECT pos_id, pos_chrm, strand, pos FROM pos')
            for pos_dbid, chrm_dbid, strand, pos in local_cursor:
                yield pos_dbid, chrm_dbid, strand, pos

    def iter_data(self):
        # use local cursor since other processing might use class cursor
        local_cursor = self.db.cursor()
        local_cursor.execute(
            'SELECT score, uuid, mod_base, motif, motif_pos, raw_motif, ' +
            'strand, pos, chrm, chrm_len ' +
            'FROM data ' +
            'INNER JOIN read ON data.score_read = read.read_id ' +
            'INNER JOIN mod ON data.score_mod = mod.mod_id ' +
            'INNER JOIN pos ON data.score_pos = pos.pos_id ' +
            'INNER JOIN chrm ON pos.pos_chrm = chrm.chrm_id')
        for data in local_cursor:
            yield data

    def close(self):
        self.db.commit()
        self.db.close()


##################
# Mod DB Helpers #
##################

def extract_all_stats(mods_db_fn):
    """ Extract all log-likelihood ratios (log(P_can / P_mod)) from a mods
    database.

    Returns:
        Dictionary with mod base single letter code keys and numpy array of
        log likelihood ratio values.
    """
    all_llrs = defaultdict(list)
    mods_db = ModsDb(mods_db_fn)
    for _, mods_pos_llrs in mods_db.iter_pos_scores():
        for mod_base, mod_pos_llrs in mods_pos_llrs.items():
            all_llrs[mod_base].append(mod_pos_llrs)
    all_llrs = dict((mod_base, np.concatenate(mod_llrs))
                    for mod_base, mod_llrs in all_llrs.items())
    return all_llrs


def extract_stats_at_valid_sites(
        mods_db_fn, valid_sites_sets, include_strand=True):
    """ Extract all log-likelihood ratios (log(P_can / P_mod)) from a mods
    database at set of valid sites.

    Args:
        mods_db_fn: Modified base database filename
        valid_sites_sets: List of sets containing valid positions. Either
            (chrm, pos) or (chrm, strand, pos); strand should be +/-1
        include_strand: Boolean value indicating whether positions include
            strand.

    Returns:
        List of dictionaries with single letter modified base code keys and
            numpy array of log likelihood ratio values. The list matches the
            order of the valid_sites_sets argument.
    """
    all_stats = [defaultdict(list) for _ in valid_sites_sets]
    # load database with positions in memory to avoid disk reads
    mods_db = ModsDb(mods_db_fn, in_mem_dbid_to_pos=True)
    for (chrm, strand, pos), mods_pos_llrs in mods_db.iter_pos_scores(
            return_pos=True):
        site_key = (chrm, strand, pos) if include_strand else (chrm, pos)
        for sites_i, valid_sites in enumerate(valid_sites_sets):
            if site_key in valid_sites:
                for mod_base, mod_pos_llrs in mods_pos_llrs.items():
                    all_stats[sites_i][mod_base].append(mod_pos_llrs)
    r_all_stats = [
        dict((mod_base, np.concatenate(mod_llrs))
             for mod_base, mod_llrs in stats_i.items())
        for stats_i in all_stats]
    return r_all_stats


########################
# Reference Mod Markup #
########################

def annotate_all_mods(r_start, ref_seq, r_mod_scores, strand, mods_info):
    """ Annotate reference sequence with called modified bases.

    Args:
        r_start (int): Reference start position for this read.
        ref_seq (str): Read-centric reference sequence corresponding to
            this read.
        r_mod_scores (list): Per-reference position modified base calls, as
            returned from megalodon.mods.call_read_mods.
        strand (int): 1 for forward strand -1 for reverse strand
        mods_info (mods.ModInfo): Object containing information about modified
            base processing

    Returns:
        mods.ANNOT_MODS object annotated with all modified bases.

    Note: Reference sequence is in read orientation and mod calls are in
    genome coordiates.
    """
    all_mods_seq, all_mods_qual = [], []
    prev_pos = 0
    if strand == -1:
        ref_seq = ref_seq[::-1]
    for mod_pos, mod_lps, mod_bases, _, _, _ in sorted(r_mod_scores):
        can_lp = np.log1p(-np.exp(mod_lps).sum())
        if can_lp - mod_lps.max() > mods_info.mod_thresh:
            base_lp = can_lp
            base = ref_seq[mod_pos - r_start]
        else:
            most_prob_mod = np.argmax(mod_lps)
            base_lp = mod_lps[most_prob_mod]
            base = mod_bases[most_prob_mod]
        if mods_info.map_base_conv is not None:
            # convert base for bisulfite-like output
            base = base.translate(mods_info.map_base_conv)
        all_mods_seq.append(
            ref_seq[prev_pos:mod_pos - r_start] + base)
        all_mods_qual.extend(
            [MOD_MAP_MAX_QUAL, ] * (mod_pos - r_start - prev_pos) +
            [min(mh.log_prob_to_phred(base_lp, False), MOD_MAP_MAX_QUAL)])
        prev_pos = mod_pos - r_start + 1

    all_mods_seq.append(ref_seq[prev_pos:])
    all_mods_seq = ''.join(all_mods_seq)
    all_mods_qual.extend([MOD_MAP_MAX_QUAL] * len(ref_seq[prev_pos:]))
    all_mods_qual = list(map(int, all_mods_qual))
    if strand == -1:
        all_mods_seq = all_mods_seq[::-1]
        all_mods_qual = all_mods_qual[::-1]

    return ANNOT_MODS(all_mods_seq, all_mods_qual)


def annotate_mods_per_mod(r_start, ref_seq, r_mod_scores, strand, mods_info):
    """ Annotate reference sequence with called modified bases. Produce one
    mods.ANNOT_MODS output for each modified base in mods_info.mod_long_names.

    Args:
        r_start (int): Reference start position for this read.
        ref_seq (str): Read-centric reference sequence corresponding to
            this read.
        r_mod_scores (list): Per-reference position modified base calls, as
            returned from megalodon.mods.call_read_mods.
        strand (int): 1 for forward strand -1 for reverse strand
        mods_info (mods.ModInfo): Object containing information about modified
            base processing

    Returns:
        Dictionary with mod_base single letter code pointing to
        mods.ANNOT_MODS with that modified base annotated.

    Note: Reference sequence is in read orientation and mod calls are in
    genome coordiates.
    """
    # seq, qual and prev_pos for each mod
    per_mod_data = dict((mod_base, [[], [], 0])
                        for mod_base, _ in mods_info.mod_long_names)
    if strand == -1:
        ref_seq = ref_seq[::-1]
    for mod_pos, mod_lps, mod_bases, _, _, _ in sorted(r_mod_scores):
        can_lp = np.log1p(-np.exp(mod_lps).sum())
        # annotate per-mod sequences and qualities
        for mod_idx, mod_base in enumerate(mod_bases):
            # called canonical
            if can_lp - mod_lps[mod_idx] > mods_info.mod_thresh:
                base_lp = can_lp
                base = ref_seq[mod_pos - r_start]
            else:
                base_lp = mod_lps[mod_idx]
                base = mod_base
            if mods_info.map_base_conv is not None:
                # convert base for bisulfite-like output
                base = base.translate(mods_info.map_base_conv)
            per_mod_data[mod_base][0].append(
                ref_seq[per_mod_data[mod_base][2]:mod_pos - r_start] + base)
            per_mod_data[mod_base][1].extend(
                [MOD_MAP_MAX_QUAL, ] *
                (mod_pos - r_start - per_mod_data[mod_base][2]) +
                [min(mh.log_prob_to_phred(base_lp, False), MOD_MAP_MAX_QUAL)])
            per_mod_data[mod_base][2] = mod_pos - r_start + 1

    per_mod_ret_data = {}
    for mod_base, (mod_seq, mod_qual, mod_prev_pos) in per_mod_data.items():
        mod_seq.append(ref_seq[mod_prev_pos:])
        mod_seq = ''.join(mod_seq)
        mod_qual.extend([MOD_MAP_MAX_QUAL] * len(ref_seq[mod_prev_pos:]))
        mod_qual = list(map(int, mod_qual))
        if strand == -1:
            mod_seq = mod_seq[::-1]
            mod_qual = mod_qual[::-1]
        per_mod_ret_data[mod_base] = ANNOT_MODS(mod_seq, mod_qual)

    return per_mod_ret_data


########################
# Per-read Mod Scoring #
########################

def score_mod_seq(
        tpost, seq, mod_cats, can_mods_offsets,
        tpost_start=0, tpost_end=None, all_paths=False):
    """ Score a section of log transition posteriors against a proposed sequence
    using a global mapping.

    Args:
        tpost `ndarray`: Log transition posteriors to be scored
        seq: `ndarray`: Integer encoded proposed sequence
        mod_cats `ndarray`: Integer encoded proposed modified base labels
        can_mods_offsets: `ndarray`: Offset into modbase transition for each
            canonical base
        tpost_start (int): Start position within post (Default: 0)
        tpost_end (int): end position within post (Default: full posterior)
        all_paths (bool): Produce the forwards all paths score
            (default: Viterbi best path)

    Returns:
        Float score representing probability of proposed sequence
    """
    seq = seq.astype(np.uintp)
    if tpost_end is None:
        tpost_end = tpost.shape[0]

    return decode.score_mod_seq(
        tpost, seq, mod_cats, can_mods_offsets, tpost_start, tpost_end,
        all_paths)


def call_read_mods(
        r_ref_pos, r_ref_seq, ref_to_block, r_post, mods_info, mod_pos_conn,
        mod_sig_map_q, sig_map_res, signal_reversed, uuid):
    def iter_motif_sites():
        search_ref_seq = r_ref_seq[::-1] if signal_reversed else r_ref_seq
        ref_seq_len = len(r_ref_seq)
        max_pos = ref_seq_len - mods_info.edge_buffer
        for motif, rel_pos, mod_bases, raw_motif in mods_info.all_mod_motifs:
            for motif_match in motif.finditer(search_ref_seq):
                m_pos = motif_match.start() + rel_pos
                if signal_reversed:
                    m_pos = ref_seq_len - m_pos - 1
                if m_pos < mods_info.edge_buffer:
                    continue
                if m_pos > max_pos:
                    break
                yield m_pos, mod_bases, motif_match.group(), rel_pos, raw_motif

    def get_pos_mod_read_dbids(r_mod_scores):
        # send unique positions and mod bases to connection to extract
        # database ids
        r_uniq_pos = sorted(set(pms[0] for pms in r_mod_scores))
        r_uniq_mod_bases = sorted(set(
            (mod_base, motif, motif_pos, raw_motif)
            for _, _, mod_bases, motif, motif_pos, raw_motif in r_mod_scores
            for mod_base in mod_bases))
        mod_pos_conn.send((
            r_uniq_pos, r_uniq_mod_bases, r_ref_pos.chrm, r_ref_pos.strand,
            uuid))
        r_pos_dbids, r_mod_dbids, read_dbid = mod_pos_conn.recv()
        # enumerate data to be added to main data table
        r_pos_dbids = dict(zip(r_uniq_pos, r_pos_dbids))
        r_mod_dbids = dict(zip(r_uniq_mod_bases, r_mod_dbids))
        return r_pos_dbids, r_mod_dbids, read_dbid

    # call all mods overlapping this read
    r_mod_scores = []
    # ignore when one or more mod_llrs is -inf (or close enough for exp)
    # occurs in compute_log_probs function, but more efficient to seterr
    # at this higher level
    with np.errstate(divide='ignore', over='ignore'):
        for (pos, mod_bases, ref_motif, rel_pos,
             raw_motif) in iter_motif_sites():
            pos_bb, pos_ab = min(mods_info.mod_context_bases, pos), min(
                mods_info.mod_context_bases, len(r_ref_seq) - pos - 1)
            try:
                pos_ref_seq = mh.seq_to_int(
                    r_ref_seq[pos - pos_bb:pos + pos_ab + 1])
            except mh.MegaError:
                ref_pos = (r_ref_pos.start + pos if r_ref_pos.strand == 1 else
                           r_ref_pos.start + len(r_ref_seq) - pos - 1)
                LOGGER.debug(
                    'Invalid sequence encountered calling modified base ' +
                    'at {}:{}'.format(r_ref_pos.chrm, ref_pos))
                continue
            pos_can_mods = np.zeros_like(pos_ref_seq)

            blk_start, blk_end = (ref_to_block[pos - pos_bb],
                                  ref_to_block[pos + pos_ab])
            if blk_end - blk_start < (mods_info.mod_context_bases * 2) + 1:
                # no valid mapping over large inserted query bases
                # i.e. need as many "events/strides" as bases for valid mapping
                continue

            loc_can_score = score_mod_seq(
                r_post, pos_ref_seq, pos_can_mods, mods_info.can_mods_offsets,
                blk_start, blk_end, mods_info.mod_all_paths)
            if loc_can_score is None:
                raise mh.MegaError('Score computation error (memory error)')

            calib_llrs = []
            for mod_base in mod_bases:
                pos_mod_mods = pos_can_mods.copy()
                pos_mod_mods[pos_bb] = mods_info.str_to_int_mod_labels[
                    mod_base]
                loc_mod_score = score_mod_seq(
                    r_post, pos_ref_seq, pos_mod_mods,
                    mods_info.can_mods_offsets, blk_start, blk_end,
                    mods_info.mod_all_paths)
                if loc_mod_score is None:
                    raise mh.MegaError(
                        'Score computation error (memory error)')

                # calibrate llr scores
                calib_llrs.append(mods_info.calibrate_llr(
                    loc_can_score - loc_mod_score, mod_base))

            # due to calibration mutli-mod log likelihoods could result in
            # inferred negative reference likelihood, so re-normalize here
            loc_mod_lps = calibration.compute_log_probs(np.array(calib_llrs))

            if (r_ref_pos.strand == 1 and not signal_reversed) or (
                    r_ref_pos.strand == -1 and signal_reversed):
                mod_ref_pos = r_ref_pos.start + pos
            else:
                mod_ref_pos = r_ref_pos.end - pos - 1
            r_mod_scores.append((
                mod_ref_pos, loc_mod_lps, mod_bases, ref_motif, rel_pos,
                raw_motif))

    all_mods_seq = per_mod_seqs = None
    if mods_info.do_ann_all_mods:
        # ignore divide around full annotate_mods call to avoid overhead
        # on many calls to errstate
        with np.errstate(divide='ignore'):
            all_mods_seq = annotate_all_mods(
                r_ref_pos.start, r_ref_seq, r_mod_scores,
                r_ref_pos.strand, mods_info)
    if mods_info.do_ann_per_mod:
        with np.errstate(divide='ignore'):
            per_mod_seqs = annotate_mods_per_mod(
                r_ref_pos.start, r_ref_seq, r_mod_scores,
                r_ref_pos.strand, mods_info)

    # send mod annotated seqs to signal mapping queue if requested
    if mod_sig_map_q is not None and sig_map_res.pass_filts:
        # import locally so that import of mods module does not require
        # taiyaki install (required for signal_mapping module)
        from megalodon import signal_mapping
        if sig_map_res.ref_out_info.annotate_mods:
            invalid_chars = set(all_mods_seq.mod_seq).difference(
                sig_map_res.ref_out_info.alphabet)
            if len(invalid_chars) > 0:
                raise mh.MegaError(
                    'Inavlid charcters found in mapped signal sequence: ' +
                    '({})'.format(''.join(invalid_chars)))
            # replace reference sequence with mod annotated sequence
            sig_map_res = sig_map_res._replace(ref_seq=all_mods_seq.mod_seq)

        mod_sig_map_q.put(signal_mapping.get_remapping(*sig_map_res[1:]))

    r_pos_dbids, r_mod_dbids, read_dbid = get_pos_mod_read_dbids(r_mod_scores)
    r_insert_data = [
        (mod_lp, r_pos_dbids[mod_pos],
         r_mod_dbids[(mod_base, ref_motif, rel_pos, raw_motif)], read_dbid)
        for mod_pos, mod_lps, mod_bases, ref_motif, rel_pos, raw_motif in
        r_mod_scores for mod_lp, mod_base in zip(mod_lps, mod_bases)]

    mod_out_text = None
    if mods_info.write_mods_txt:
        txt_tmplt = '\t'.join('{}' for _ in ModsDb.text_field_names)
        mod_out_text = ''
        for (pos, mod_lps, mod_bases, ref_motif, rel_pos,
             raw_motif) in r_mod_scores:
            with np.errstate(divide='ignore'):
                can_lp = np.log1p(-np.exp(mod_lps).sum())
            mod_out_text += '\n'.join((
                txt_tmplt.format(
                    uuid, r_ref_pos.chrm, r_ref_pos.strand, pos, mod_lp,
                    can_lp, mod_base, '{}:{}'.format(raw_motif, rel_pos))
                for mod_lp, mod_base in zip(mod_lps, mod_bases))) + '\n'

    return r_insert_data, all_mods_seq, per_mod_seqs, mod_out_text


#######################
# Per-read Mod Output #
#######################

def init_mods_db(mods_db_fn, db_safety, ref_names_and_lens, mods_info):
    mods_db = ModsDb(
        mods_db_fn, db_safety=db_safety, read_only=False, init_db_tables=True)
    mods_db.insert_chrms(ref_names_and_lens)
    mods_db.insert_mod_long_names(mods_info.mod_long_names)
    mods_db.db.commit()
    LOGGER.debug((
        'mod_getter: in_mem_indices:\n\t' +
        'chrm -> dbid : {}\n\tdbid -> chrm : {}\n\t' +
        'pos  -> dbid : {}\n\tdbid -> pos  : {}\n\t' +
        'mod  -> dbid : {}\n\tdbid -> mod  : {}\n\t' +
        'uuid -> dbid : {}\n\tdbid -> uuid : {}').format(
            mods_db.in_mem_chrm_to_dbid, mods_db.in_mem_dbid_to_chrm,
            mods_db.in_mem_pos_to_dbid, mods_db.in_mem_dbid_to_pos,
            mods_db.in_mem_mod_to_dbid, mods_db.in_mem_dbid_to_mod,
            mods_db.in_mem_uuid_to_dbid, mods_db.in_mem_dbid_to_uuid))
    mods_db.close()


def _get_mods_queue(
        mods_q, mods_conn, mods_db_fn, db_safety, ref_names_and_lens, ref_fn,
        mods_txt_fn, pr_refs_fn, ref_out_info, mod_map_fns, map_fmt,
        mods_info):
    def write_mod_alignment(
            read_id, mod_seq, mod_quals, chrm, strand, r_st, fp):
        a = pysam.AlignedSegment()
        a.query_name = read_id
        a.flag = 0 if strand == 1 else 16
        a.reference_id = fp.get_tid(chrm)
        a.reference_start = r_st
        a.template_length = len(mod_seq)
        a.mapping_quality = MOD_MAP_MAX_QUAL
        a.set_tags([('RG', MOD_MAP_RG_ID)])

        # convert to reference based sequence
        if strand == -1:
            mod_seq = mh.revcomp(mod_seq)
            mod_quals = mod_quals[::-1]
        a.query_sequence = mod_seq
        a.query_qualities = array('B', mod_quals)
        a.cigartuples = [(0, len(mod_seq)), ]
        fp.write(a)

    def store_mod_call(
            r_insert_data, mod_out_text, all_mods_seq, per_mod_seqs, read_id,
            chrm, strand, r_start, ref_seq, read_len, q_st, q_en, cigar,
            been_warned):
        try:
            mods_db.insert_read_data(r_insert_data)
            mods_db.db.commit()
        except Exception as e:
            if not been_warned:
                LOGGER.warning(
                    'Error inserting modified base scores into database. ' +
                    'See log debug output for error details.')
                been_warned = True
            import traceback
            LOGGER.debug(
                'Error inserting modified base scores into database: ' +
                str(e) + '\n' + traceback.format_exc())

        if mods_txt_fp is not None and len(r_mod_scores) > 0:
            mods_txt_fp.write(mod_out_text)
        if (pr_refs_fn is not None and mapping.read_passes_filters(
                ref_out_info, read_len, q_st, q_en, cigar)):
            pr_refs_fp.write('>{}\n{}\n'.format(read_id, all_mods_seq.mod_seq))
        if mod_map_fns is not None:
            for mod_base, _ in mods_info.mod_long_names:
                mod_seq, mod_qual = per_mod_seqs[mod_base]
                write_mod_alignment(
                    read_id, mod_seq, mod_qual, chrm, strand, r_start,
                    mod_map_fps[mod_base])

        return been_warned

    been_warned = False

    mods_db = ModsDb(mods_db_fn, db_safety=db_safety, read_only=False)

    if mods_txt_fn is None:
        mods_txt_fp = None
    else:
        mods_txt_fp = open(mods_txt_fn, 'w')
        mods_txt_fp.write('\t'.join(mods_db.text_field_names) + '\n')
    if pr_refs_fn is not None:
        pr_refs_fp = open(pr_refs_fn, 'w')
    if mod_map_fns is not None:
        try:
            w_mode = mh.MAP_OUT_WRITE_MODES[map_fmt]
        except KeyError:
            raise mh.MegaError('Invalid mapping output format')
        header = {
            'HD': {'VN': '1.4'},
            'SQ': [{'LN': ref_len, 'SN': ref_name}
                   for ref_name, ref_len in sorted(zip(*ref_names_and_lens))],
            'RG': [{'ID': MOD_MAP_RG_ID, 'SM': SAMPLE_NAME}, ]}
        mod_map_fps = dict((
            (mod_base, pysam.AlignmentFile(
                mod_map_fn + map_fmt, w_mode, header=header,
                reference_filename=ref_fn))
            for mod_base, mod_map_fn in mod_map_fns))

    LOGGER.debug('mod_getter: init complete')

    while True:
        try:
            # note strand is +1 for fwd or -1 for rev
            (r_mod_scores, all_mods_seq, per_mod_seqs, mod_out_text), (
                read_id, chrm, strand, r_start, ref_seq, read_len, q_st, q_en,
                cigar) = mods_q.get(block=True, timeout=1)
        except queue.Empty:
            if mods_conn.poll():
                break
            continue
        try:
            been_warned = store_mod_call(
                r_mod_scores, mod_out_text, all_mods_seq, per_mod_seqs,
                read_id, chrm, strand, r_start, ref_seq, read_len, q_st, q_en,
                cigar, been_warned)
        except Exception as e:
            LOGGER.debug('Error processing mods output for read: ' +
                         '{}\nError type: {}'.format(read_id, str(e)))

    while not mods_q.empty():
        (r_mod_scores, all_mods_seq, per_mod_seqs, mod_out_text), (
            read_id, chrm, strand, r_start, ref_seq, read_len, q_st, q_en,
            cigar) = mods_q.get(block=False)
        try:
            been_warned = store_mod_call(
                r_mod_scores, mod_out_text, all_mods_seq, per_mod_seqs,
                read_id, chrm, strand, r_start, ref_seq, read_len, q_st, q_en,
                cigar, been_warned)
        except Exception as e:
            LOGGER.debug('Error processing mods output for read: ' +
                         '{}\nError type: {}'.format(read_id, str(e)))

    if mods_txt_fp is not None:
        mods_txt_fp.close()
    if pr_refs_fn is not None:
        pr_refs_fp.close()
    if mod_map_fns is not None:
        for mod_map_fp in mod_map_fps.values():
            mod_map_fp.close()

    mods_db.create_data_covering_index()
    mods_db.close()


if _PROFILE_MODS_QUEUE:
    _get_mods_queue_wrapper = _get_mods_queue

    def _get_mods_queue(*args):
        import cProfile
        cProfile.runctx('_get_mods_queue_wrapper(*args)', globals(), locals(),
                        filename='mods_getter_queue.prof')


def _mod_aux_table_inserts(mod_db_fn, db_safety, pos_in_mem, mod_pos_conns):
    mods_db = ModsDb(
        mod_db_fn, db_safety=db_safety, in_mem_pos_to_dbid=pos_in_mem,
        in_mem_dbid_to_chrm=True, in_mem_mod_to_dbid=True, read_only=False)
    # loop over connections to read worker processes until all have been
    # exhausted
    while mod_pos_conns:
        for r in wait(mod_pos_conns):
            try:
                conn_res = r.recv()
            except EOFError:
                mod_pos_conns.remove(r)
            else:
                r_uniq_pos, r_uniq_mod_bases, chrm, strand, uuid = conn_res
                chrm_dbid = mods_db.get_chrm_dbid(chrm)
                r_pos_dbids = mods_db.get_pos_dbids_or_insert(
                    r_uniq_pos, chrm_dbid, strand)
                r_mod_dbids = mods_db.get_mod_base_ids_or_insert(
                    r_uniq_mod_bases)
                read_dbid = mods_db.insert_read_uuid(uuid)
                mods_db.db.commit()
                r.send((r_pos_dbids, r_mod_dbids, read_dbid))

    mods_db.close()


if _PROFILE_MODS_AUX:
    _mod_aux_table_inserts_wrapper = _mod_aux_table_inserts

    def _mod_aux_table_inserts(*args):
        import cProfile
        cProfile.runctx('_mod_aux_table_inserts_wrapper(*args)',
                        globals(), locals(),
                        filename='mods_aux_inserts.prof')


############
# Mod Info #
############

class ModInfo(object):
    def distinct_bases(self, b1, b2):
        return len(set(mh.SINGLE_LETTER_CODE[b1]).intersection(
            mh.SINGLE_LETTER_CODE[b2])) == 0

    def distinct_motifs(self):
        if len(self.all_mod_motifs) in (0, 1):
            return True
        for n1, (_, rel_pos1, _, raw_motif1) in enumerate(self.all_mod_motifs):
            for (_, rel_pos2, _, raw_motif2) in self.all_mod_motifs[n1 + 1:]:
                # compute overlapping positions relative to modified position
                bb = min(rel_pos1, rel_pos2)
                ab = min(len(raw_motif1) - rel_pos1 - 1,
                         len(raw_motif2) - rel_pos2 - 1)
                if all(not self.distinct_bases(b1, b2) for b1, b2 in zip(
                        raw_motif1[rel_pos1 - bb:rel_pos1 + ab + 1],
                        raw_motif2[rel_pos2 - bb:rel_pos2 + ab + 1])):
                    return False
        return True

    def _parse_map_base_conv(self):
        if self.map_base_conv_raw is None:
            self.map_base_conv = None
        else:
            # create conversion dicationary for translate function
            from_bases = [ord(fb) for fb, tb in self.map_base_conv_raw]
            to_bases = [ord(tb) for fb, tb in self.map_base_conv_raw]
            self.map_base_conv = dict(zip(from_bases, to_bases))

    def _parse_mod_motifs(self, all_mod_motifs_raw):
        # note only works for mod_refactor models currently
        self.all_mod_motifs = []
        if all_mod_motifs_raw is None or len(all_mod_motifs_raw) == 0:
            for can_base, mod_bases in self.can_base_mods.items():
                self.all_mod_motifs.append((
                    re.compile(can_base), 0, mod_bases, can_base))
        else:
            # parse detection motifs
            for mod_bases, raw_motif, pos in all_mod_motifs_raw:
                for mod_base in mod_bases:
                    assert mod_base in self.str_to_int_mod_labels, (
                        'Modified base label ({}) not found in model ' +
                        'alphabet ({}).').format(
                            mod_base, list(self.str_to_int_mod_labels.keys()))
                pos = int(pos)
                can_base = next(
                    can_base for can_base, can_mods in
                    self.can_base_mods.items() if mod_bases[0] in can_mods)
                assert (can_base == raw_motif[pos]), (
                    'Invalid modified base motif. Raw motif modified ' +
                    'position ({}) base ({}) does not match ' +
                    'collapsed alphabet value ({}).').format(
                        pos, raw_motif[pos], can_base)
                motif = re.compile(''.join(
                    '[{}]'.format(mh.SINGLE_LETTER_CODE[letter])
                    for letter in raw_motif))
                self.all_mod_motifs.append((motif, pos, mod_bases, raw_motif))

            if not self.distinct_motifs():
                raise mh.MegaError(
                    'One provided motif can be found within another motif. ' +
                    'Only distinct sets of motifs are accepted')

    def __init__(
            self, model_info, all_mod_motifs_raw=None, mod_all_paths=False,
            write_mods_txt=None, mod_context_bases=None,
            do_output_mods=False, mods_calib_fn=None,
            mod_output_fmts=[mh.MOD_BEDMETHYL_NAME],
            edge_buffer=mh.DEFAULT_EDGE_BUFFER, pos_index_in_memory=True,
            agg_info=DEFAULT_AGG_INFO, mod_thresh=0.0, do_ann_all_mods=False,
            do_ann_per_mod=False, map_base_conv=None):
        # this is pretty hacky, but these attributes are stored here as
        # they are generally needed alongside other modbase info
        # don't want to pass all of these parameters around individually though
        # as this would make function signatures too complicated
        self.mod_all_paths = mod_all_paths
        self.write_mods_txt = write_mods_txt
        self.mod_context_bases = mod_context_bases
        self.do_output_mods = do_output_mods
        self.mod_long_names = model_info.mod_long_names
        self.calib_table = calibration.ModCalibrator(mods_calib_fn)
        self.mod_output_fmts = mod_output_fmts
        self.edge_buffer = edge_buffer
        self.pos_index_in_memory = pos_index_in_memory
        self.agg_info = agg_info
        self.mod_thresh = mod_thresh
        self.do_ann_all_mods = do_ann_all_mods
        self.do_ann_per_mod = do_ann_per_mod
        self.map_base_conv_raw = map_base_conv

        self.alphabet = model_info.can_alphabet
        self.ncan_base = len(self.alphabet)
        try:
            self.alphabet = self.alphabet.decode()
        except AttributeError:
            pass
        if model_info.is_cat_mod:
            # TODO also output "(alt to C)" for each mod
            LOGGER.info(
                'Using canonical alphabet {} and modified bases {}.'.format(
                    self.alphabet, ' '.join(
                        '{}={}'.format(*mod_b)
                        for mod_b in model_info.mod_long_names)))
        else:
            LOGGER.info(
                'Using canonical alphabet {}.'.format(self.alphabet))

        self.nbase = len(self.alphabet)
        self.n_can_state = (self.ncan_base + self.ncan_base) * (
            self.ncan_base + 1)
        if model_info.is_cat_mod:
            self.nmod_base = model_info.n_mods
            self.can_base_mods = model_info.can_base_mods
            self.can_mods_offsets = model_info.can_indices
            self.str_to_int_mod_labels = model_info.str_to_int_mod_labels
            assert (
                model_info.output_size - self.n_can_state ==
                self.nmod_base + 1), (
                    'Alphabet ({}) and model number of modified bases ({}) ' +
                    'do not agree.').format(
                        self.alphabet,
                        model_info.output_size - self.n_can_state - 1)
        else:
            self.nmod_base = 0
            self.can_base_mods = {}
            self.can_mods_offsets = None
            self.str_to_int_mod_labels = None

        # parse mod motifs or use "swap" base if no motif provided
        self._parse_mod_motifs(all_mod_motifs_raw)
        self._parse_map_base_conv()

    def calibrate_llr(self, score, mod_base):
        return self.calib_table.calibrate_llr(score, mod_base)


#################
# modVCF Writer #
#################

class ModSite(object):
    """ Modified base site for entry into Mod Writers.
    Currently only handles a single sample.
    """
    def __init__(
            self, chrom, pos, strand, ref_seq, mod_bases,
            id='.', qual='.', filter='.', info=None, sample_dict=None,
            ref_mod_pos=0, mod_props=None):
        self.strand = strand
        self.mod_bases = mod_bases
        self.mod_props = mod_props
        mod_seqs = ','.join((
            ref_seq[:ref_mod_pos] + mod_base + ref_seq[ref_mod_pos + 1:]
            for mod_base in mod_bases))

        # attributes must match header text from ModVcfWriter
        self.chrom = chrom
        self.pos = int(pos)
        self.id = str(id)
        self.ref = ref_seq.upper()
        self.alt = mod_seqs
        self.qual = qual
        self.filter = str(filter)

        # info and gentype data fields
        if info is None:
            info = {}
        if mh.STRAND_FIELD_NAME not in info:
            info[mh.STRAND_FIELD_NAME] = strand
        self.info_dict = info
        if sample_dict is None:
            sample_dict = OrderedDict()
        self.sample_dict = sample_dict

        if self.mod_props is not None:
            self.add_mod_props(self.mod_props)

    @property
    def _sorted_format_keys(self):
        sorted_keys = sorted(self.sample_dict.keys())
        if 'LOG_PROBS' in sorted_keys:
            # move log probs to end of format field for easier
            # human readability
            sorted_keys.append(sorted_keys.pop(sorted_keys.index('LOG_PROBS')))
        return sorted_keys

    @property
    def format(self):
        return ':'.join(map(str, self._sorted_format_keys))

    @property
    def sample(self):
        return ':'.join((str(self.sample_dict[k])
                         for k in self._sorted_format_keys))

    @property
    def info(self):
        str_tags = []
        for key, value in self.info_dict.items():
            # If key is of type 'Flag', print only key, else 'key=value'
            if value is True:
                str_tags.append(key)
            else:
                if isinstance(value, (tuple, list)):
                    value = ','.join(map(str, value))
                str_tags.append('{}={}'.format(key, value))
        return ':'.join(str_tags)

    def add_tag(self, tag, value=None):
        self.info_dict[tag] = value

    def add_sample_field(self, tag, value=None):
        self.sample_dict[tag] = value

    def add_mod_props(self, mod_props):
        with np.errstate(divide='ignore'):
            can_pl = -10 * np.log10(1 - sum(mod_props.values()))
        self.qual = '{:.0f}'.format(
            np.abs(np.around(np.minimum(can_pl, mh.MAX_PL_VALUE))))
        for mod_name, mod_prop in mod_props.items():
            self.add_sample_field(mod_name, '{:.4f}'.format(mod_prop))
        return

    def get_coverage(self, default_value=0):
        if 'VALID_DP' in self.sample_dict:
            return self.sample_dict['VALID_DP']
        elif 'DP' in self.sample_dict:
            return self.sample_dict['DP']
        elif 'DP' in self.info_dict:
            return self.info_dict['DP']
        return default_value

    def __eq__(self, mod2):
        return (self.chrm, self.pos, self.strand) == (
            mod2.chrm, mod2.pos, mod2.strand)

    def __ne__(self, var2):
        return (self.chrm, self.pos, self.strand) != (
            var2.chrm, var2.pos, var2.strand)

    def __lt__(self, var2):
        return (self.chrm, self.pos, self.strand) < (
            var2.chrm, var2.pos, var2.strand)

    def __le__(self, var2):
        return (self.chrm, self.pos, self.strand) <= (
            var2.chrm, var2.pos, var2.strand)

    def __gt__(self, var2):
        return (self.chrm, self.pos, self.strand) > (
            var2.chrm, var2.pos, var2.strand)

    def __ge__(self, var2):
        return (self.chrm, self.pos, self.strand) >= (
            var2.chrm, var2.pos, var2.strand)


class ModVcfWriter(object):
    """ modVCF writer class
    """
    version_options = set(['4.2', ])

    def __init__(
            self, basename, mods, mode='w',
            header=('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
                    'INFO', 'FORMAT', 'SAMPLE'),
            extra_meta_info=FIXED_VCF_MI, version='4.2', ref_fn=None,
            ref_names_and_lens=None, write_mod_lp=False,
            buffer_limit=OUT_BUFFER_LIMIT):
        self.basename = basename
        self.mods = mods
        self.mode = mode
        self.buffer_limit = buffer_limit
        self.buffer = []
        self.header = header
        if version not in self.version_options:
            raise ValueError('version must be one of {}'.format(
                self.version_options))
        self.version = version
        contig_mis = [] if ref_names_and_lens is None else [
            mh.CONTIG_MI.format(ref_name, ref_len)
            for ref_name, ref_len in sorted(zip(*ref_names_and_lens))]
        self.meta = [
            mh.VCF_VERSION_MI.format(self.version),
            mh.FILE_DATE_MI.format(
                datetime.date.today().strftime("%Y%m%d")),
            mh.SOURCE_MI.format(MEGALODON_VERSION),
            mh.REF_MI.format(ref_fn)] + contig_mis + extra_meta_info + [
                mod_tmplt.format(*mod_name) for mod_name in self.mods
                for mod_tmplt in MOD_MI_TMPLTS]
        if write_mod_lp:
            self.meta.append(FORMAT_LOG_PROB_MI)
        self.filename = '{}.{}'.format(
            self.basename, mh.MOD_OUTPUT_EXTNS[mh.MOD_VCF_NAME])
        self.handle = open(self.filename, self.mode, encoding='utf-8')
        self.handle.write('\n'.join('##' + line for line in self.meta) + '\n')
        self.handle.write('#' + '\t'.join(self.header) + '\n')

    def write_mod_site(self, mod_site):
        elements = [getattr(mod_site, field.lower()) for field in self.header]
        elements = ['.' if e == '' else e for e in elements]
        # VCF POS field is 1-based
        elements[self.header.index('POS')] += 1
        self.buffer.append('\t'.join(map(str, elements)))
        if len(self.buffer) > self.buffer_limit:
            self.handle.write('\n'.join(self.buffer) + '\n')
            self.buffer = []

    def close(self):
        if len(self.buffer) > 0:
            self.handle.write('\n'.join(self.buffer) + '\n')
        self.handle.close()


class ModBedMethylWriter(object):
    """ bedMethyl writer class

    Note that the bedMethyl format cannot store more than one modification
    type, so multiple file handles will be opened.
    """
    def __init__(
            self, basename, mods, mode='w', buffer_limit=OUT_BUFFER_LIMIT):
        self.basename = basename
        self.mods = mods
        self.mod_short_names, self.mod_long_names = zip(*self.mods)
        self.mode = mode
        self.buffer_limit = buffer_limit
        self.buffers = dict(
            (mod_short_name, []) for mod_short_name, _ in self.mods)
        self.handles = dict(
            (mod_short_name, open('{}.{}.{}'.format(
                self.basename, mod_long_name,
                mh.MOD_OUTPUT_EXTNS[mh.MOD_BEDMETHYL_NAME]),
                                  self.mode, encoding='utf-8'))
            for mod_short_name, mod_long_name in self.mods)

    def write_mod_site(self, mod_site):
        for mod_base, mod_prop in mod_site.mod_props.items():
            if mod_base not in self.mod_short_names:
                mh.warning('Invalid modified base encountered during ' +
                           'bedMethyl output.')
                continue

            cov = mod_site.get_coverage()
            self.buffers[mod_base].append(
                ('{chrom}\t{pos}\t{end}\t.\t{score}\t{strand}\t{pos}' +
                 '\t{end}\t0,0,0\t{cov}\t{perc}').format(
                     chrom=mod_site.chrom, pos=mod_site.pos,
                     end=mod_site.pos + 1, strand=mod_site.strand, cov=cov,
                     score=min(int(cov), 1000), perc=int(mod_prop * 100)))
            if len(self.buffers[mod_base]) > self.buffer_limit:
                self.handles[mod_base].write(
                    '\n'.join(self.buffers[mod_base]) + '\n')
                self.buffers[mod_base] = []

    def close(self):
        for mod_base, handle in self.handles.items():
            if len(self.buffers[mod_base]) > 0:
                handle.write(
                    '\n'.join(self.buffers[mod_base]) + '\n')
            handle.close()


class ModWigWriter(object):
    """ Modified base wiggle variableStep writer class

    Note that the wiggle/bedgraph format cannot store more than one
    modification type or multiple strands, so multiple file handles will
    be opened.
    """
    def __init__(
            self, basename, mods, mode='w',
            strands={'+': 'fwd_strand', '-': 'rev_strand'}):
        self.basename = basename
        self.mods = mods
        self.mods_lookup = dict(mods)
        self.mod_short_names, self.mod_long_names = zip(*self.mods)
        self.mode = mode
        self.strands = strands

        self.mod_sites_data = dict(
            ((mod_short_name, strand), defaultdict(list))
            for mod_short_name, mod_long_name in self.mods
            for strand, strand_name in strands.items())

    def write_mod_site(self, mod_site):
        if mod_site.strand not in self.strands:
            mh.warning('Invalid strand encountered during wiggle output.')
            return
        for mod_base, mod_prop in mod_site.mod_props.items():
            if mod_base not in self.mod_short_names:
                mh.warning('Invalid modified base encountered during ' +
                           'wiggle output.')
                continue
            self.mod_sites_data[(mod_base, mod_site.strand)][
                mod_site.chrom].append((mod_site.pos, mod_prop))

    def close(self):
        # write all data on close since all data is required to write
        # wiggle format
        for (mod_base, strand), all_cs_mod_sites in \
              self.mod_sites_data.items():
            with open('{}.{}.{}.{}'.format(
                    self.basename, self.mods_lookup[mod_base],
                    self.strands[strand],
                    mh.MOD_OUTPUT_EXTNS[mh.MOD_WIG_NAME]),
                      self.mode, encoding='utf-8') as wig_fp:
                # write header
                track_name = ('Modified Base {} Proportion Modified ' +
                              '({})').format(
                                  self.mods_lookup[mod_base],
                                  self.strands[strand])
                wig_fp.write((
                    'track type=wiggle_0 name="{0}" ' +
                    'description="{0}"\n').format(track_name))
                for chrom, cs_mod_sites in all_cs_mod_sites.items():
                    wig_fp.write(
                        'variableStep chrom={} span=1\n'.format(chrom))
                    wig_fp.write('\n'.join((
                        '{} {}'.format(pos, mod_prop)
                        for pos, mod_prop in sorted(cs_mod_sites))) + '\n')


##########################
# Mods Aggregation Class #
##########################

class AggMods(mh.AbstractAggregationClass):
    """ Class to assist in database queries for per-site aggregation of
    modified base calls over reads.

    Warning, setting pos_index_in_memory for a large database will drastically
    increase the startup time and memory usage (~8GB for a human genome at
    CpG sites).
    """
    def __init__(self, mods_db_fn, agg_info=DEFAULT_AGG_INFO,
                 write_mod_lp=False, load_uuid_index_in_memory=False,
                 no_indices_in_mem=False):
        # open as read only database (default)
        if no_indices_in_mem:
            self.mods_db = ModsDb(mods_db_fn)
        if load_uuid_index_in_memory:
            self.mods_db = ModsDb(
                mods_db_fn, in_mem_dbid_to_chrm=True, in_mem_dbid_to_mod=True,
                in_mem_dbid_to_uuid=True)
        else:
            self.mods_db = ModsDb(
                mods_db_fn, in_mem_dbid_to_chrm=True, in_mem_dbid_to_mod=True)
        self.n_uniq_mods = None
        assert agg_info.method in mh.MOD_AGG_METHOD_NAMES
        self.agg_method = agg_info.method
        self.binary_thresh = agg_info.binary_threshold
        self.write_mod_lp = write_mod_lp
        self._mod_long_names = self.mods_db.get_mod_long_names()

    def get_mod_long_names(self):
        if self._mod_long_names is None:
            self._mod_long_names = self.mods_db.get_mod_long_names()
        return self._mod_long_names

    def num_uniq(self):
        if self.n_uniq_mods is None:
            self.n_uniq_mods = self.mods_db.get_num_uniq_mod_pos()
        return self.n_uniq_mods

    def iter_uniq(self):
        # fill queue with full position information to make
        # workers avoid the ordered pos data extraction
        for q_val in self.mods_db.iter_pos():
            yield q_val

    def est_binary_thresh(self, pos_scores):
        valid_cov = 0
        mod_types = set(mt for read_mods in pos_scores.values()
                        for mt in read_mods.keys())
        mods_cov = dict((mt, 0) for mt in mod_types)
        for read_pos_lps in pos_scores.values():
            mt_lps = np.array(list(read_pos_lps.values()))
            with np.errstate(divide='ignore'):
                can_lp = np.log1p(-np.exp(mt_lps).sum())
            if can_lp > mt_lps.max():
                if np.exp(can_lp) > self.binary_thresh:
                    valid_cov += 1
            else:
                if np.exp(mt_lps.max()) > self.binary_thresh:
                    valid_cov += 1
                    mods_cov[list(read_pos_lps.keys())[np.argmax(mt_lps)]] += 1

        if valid_cov == 0:
            return mods_cov, valid_cov
        mods_props = OrderedDict(sorted(
            (mod_type, mod_cov / valid_cov)
            for mod_type, mod_cov in mods_cov.items()))
        return mods_props, valid_cov

    def est_em_prop(
            self, pos_scores, max_iters=15, conv_tol=0.001, min_prop=0.001):
        """Estimate proportion of modified bases at a position via EM
        computation
        """
        max_prop = 1.0 - min_prop

        def clip(mix_props):
            """Clip proportions to specified range, while maintaining sum to 1
            """
            lower_clip_idx = np.less(mix_props, min_prop)
            upper_clip_idx = np.greater(mix_props, max_prop)
            unclip_idx = np.logical_not(
                np.logical_or(lower_clip_idx, upper_clip_idx))
            if unclip_idx.sum() == mix_props.shape[0]:
                return mix_props
            if unclip_idx.sum() > 0:
                # compute amount to be added to unclipped probabilties
                unclip_diff = ((-min_prop * lower_clip_idx.sum()) + (
                    min_prop * upper_clip_idx.sum())) / unclip_idx.sum()
                mix_props[unclip_idx] += unclip_diff
            mix_props[lower_clip_idx] = min_prop
            mix_props[upper_clip_idx] = max_prop
            return mix_props

        def unclip(mix_props):
            """Set proportions at or outside clip range to 0 or 1
            """
            lower_clip_idx = np.less_equal(mix_props, min_prop)
            upper_clip_idx = np.greater_equal(mix_props, max_prop)
            unclip_idx = np.logical_not(
                np.logical_or(lower_clip_idx, upper_clip_idx))
            if unclip_idx.sum() == mix_props.shape[0]:
                return mix_props
            if unclip_idx.sum() > 0:
                unclip_diff = (
                    mix_props[lower_clip_idx].sum() +
                    (mix_props[upper_clip_idx] - 1.0).sum()) / unclip_idx.sum()
                mix_props[unclip_idx] += unclip_diff
            mix_props[lower_clip_idx] = 0.0
            mix_props[upper_clip_idx] = 1.0
            return mix_props

        mod_types = sorted(set(mt for read_mods in pos_scores.values()
                               for mt in read_mods.keys()))
        mt_probs = np.exp(np.array([[r_mods[mt] for mt in mod_types]
                                    for r_mods in pos_scores.values()]))
        can_probs = 1 - mt_probs.sum(1)
        all_probs = np.column_stack([can_probs, mt_probs])
        curr_mix_props = clip(
            np.bincount(all_probs.argmax(axis=1),
                        minlength=all_probs.shape[1]) / all_probs.shape[0])
        for _ in range(max_iters):
            prev_mix_props = curr_mix_props.copy()
            curr_mix_props = all_probs * curr_mix_props
            curr_mix_props = clip(np.mean(curr_mix_props.transpose() /
                                          curr_mix_props.sum(1), axis=1))
            if np.abs(curr_mix_props - prev_mix_props).max() < conv_tol:
                break

        mods_props = OrderedDict(zip(mod_types, unclip(curr_mix_props)[1:]))
        return mods_props, len(pos_scores)

    def compute_mod_stats(self, mod_pos, agg_method=None, valid_read_ids=None):
        if agg_method is None:
            agg_method = self.agg_method
        if agg_method not in mh.MOD_AGG_METHOD_NAMES:
            raise NotImplementedError(
                'No modified base proportion estimation method: {}'.format(
                    agg_method))

        pr_mod_stats = self.mods_db.get_pos_stats(
            mod_pos, return_uuids=valid_read_ids is not None)
        mod_type_stats = defaultdict(dict)
        for r_stats in pr_mod_stats:
            if valid_read_ids is not None and \
               r_stats.read_id not in valid_read_ids:
                continue
            mod_type_stats[r_stats.read_id][r_stats.mod_base] = r_stats.score
        total_cov = len(mod_type_stats)
        if total_cov == 0:
            raise mh.MegaError('No valid reads cover modified base location')
        if agg_method == mh.MOD_BIN_THRESH_NAME:
            mod_props, valid_cov = self.est_binary_thresh(mod_type_stats)
        elif agg_method == mh.MOD_EM_NAME:
            mod_props, valid_cov = self.est_em_prop(mod_type_stats)

        r0_stats = pr_mod_stats[0]
        strand = mh.int_strand_to_str(r0_stats.strand)
        mod_site = ModSite(
            chrom=r0_stats.chrm, pos=r0_stats.pos, strand=strand,
            ref_seq=r0_stats.motif, ref_mod_pos=r0_stats.motif_pos,
            mod_bases=list(mod_props.keys()), mod_props=mod_props)
        mod_site.add_tag('DP', '{}'.format(total_cov))
        mod_site.add_sample_field('DP', '{}'.format(total_cov))
        mod_site.add_sample_field('VALID_DP', '{}'.format(int(valid_cov)))

        if self.write_mod_lp:
            mods_lps = [[] for _ in mod_props]
            for read_mod_scores in mod_type_stats.values():
                try:
                    for mod_i, mod_lp in enumerate([
                            read_mod_scores[mod_type]
                            for mod_type in mod_props]):
                        mods_lps[mod_i].append(mod_lp)
                except KeyError:
                    continue
            mod_site.add_sample_field('LOG_PROBS', ','.join(
                ';'.join('{:.2f}'.format(lp) for lp in mod_i_lps)
                for mod_i_lps in mods_lps))

        return mod_site

    def close(self):
        self.mods_db.close()


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
