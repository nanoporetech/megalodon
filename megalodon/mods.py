import os
import re
import sys
import array
import pysam
import queue
import sqlite3
import datetime
import traceback
from collections import defaultdict, namedtuple, OrderedDict

import numpy as np

from megalodon import (
    calibration, decode, logging, mapping, megalodon_helper as mh)
from megalodon._version import MEGALODON_VERSION


LOGGER = logging.get_logger()

AGG_INFO = namedtuple('AGG_INFO', (
    'method', 'binary_threshold', 'expit_k', 'expit_x0'))
AGG_INFO.__new__.__defaults__ = (mh.DEFAULT_MOD_LK, mh.DEFAULT_MOD_L0)
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
BEDMETHYL_TMPLT = ('{chrom}\t{pos}\t{end}\t.\t{score}\t{strand}\t{pos}' +
                   '\t{end}\t0,0,0\t{cov}\t{perc:.1f}')

SAMPLE_NAME = 'SAMPLE'
MOD_MAP_RG_ID = '1'
MOD_MAP_INVALID_BASE_LP = -2
MOD_MAP_MAX_QUAL = 40

POS_IDX_CHNG_ERR_MSG = (
    'Inserting chromosomes forces change in in-memory index size. Please ' +
    'insert all chromosomes at database initialization.')
# trailing newlines allow message to be seen even after dynamic progress output
DB_TIMEOUT_ERR_MSG = (
    'Modified base database {} insert failed, likely due to timeout.\n' +
    'Potential fixes: move output to fast disk, set `--database-safety 0`, ' +
    'increase --mod-database-timeout.\nFuture failures will be logged ' +
    'without warnings, but will trigger re-establishment of database ' +
    'connection so the source issue should be resolved.\n{}' + ('\n' * 6))

OUT_BUFFER_LIMIT = 10000

_PROFILE_MODS_QUEUE = False


###########
# Mods DB #
###########

class ModsDb:
    """ Interface to the SQLite database containing per-read modified base
    statistics. When first creating a new modified base database it is highly
    recommended to use the  `mods.init_mods_db` helper command.

    Once a database is initialized the default object initialization opens a
    read-only connection. Setting `read_only=False` opens the connection for
    writing modified base statistics.
    """
    # TODO annotate available attributes/lookup tables after initialization

    # Note foreign key constraint is not applied here as this
    # drastically reduces efficiency. Thus foreign key constraint must be
    # handled by the class.
    db_tables = OrderedDict((
        ('chrm', OrderedDict((
            ('chrm_id', 'INTEGER PRIMARY KEY'),
            ('chrm', 'TEXT'),
            ('chrm_len', 'INTEGER')))),
        ('mod_long_names', OrderedDict((
            ('mod_id', 'INTEGER PRIMARY KEY'),
            ('mod_base', 'TEXT'),
            ('can_base', 'TEXT'),
            ('mod_long_name', 'TEXT')))),
        ('read', OrderedDict((
            ('read_id', 'INTEGER PRIMARY KEY'),
            ('uuid', 'TEXT')))),
        ('data', OrderedDict((
            ('score', 'FLOAT'),
            ('score_pos', 'INTEGER'),
            ('score_mod', 'INTEGER'),
            ('score_read', 'INTEGER'))))))

    # namedtuple for returning mods from a single position
    mod_data = namedtuple('mod_data', [
        'read_id', 'chrm', 'strand', 'pos', 'score', 'mod_base'])
    text_field_names = (
        'read_id', 'chrm', 'strand', 'pos', 'mod_log_prob',
        'can_log_prob', 'mod_base')

    def __init__(self, fn, init_db_tables=False, read_only=True, db_safety=0,
                 in_mem_uuid_to_dbid=False, in_mem_dbid_to_uuid=False,
                 mod_db_timeout=mh.DEFAULT_MOD_DATABASE_TIMEOUT):
        """ Interface to database containing modified base statistics.

        Default settings are for read_only performance with in-memory indices
        for chromosomes, and modified bases. Reference position lookup is
        computed from the in-memory chromosome index.

        If read database IDs or read UUIDs are to be accessed repeatedly it is
        strongly suggested that this value be loaded at database initialization
        by setting in_mem_uuid_to_dbid or in_mem_dbid_to_uuid to True.
        """
        self.fn = mh.resolve_path(fn)
        self.init_db_tables = init_db_tables
        self.read_only = read_only
        self.db_safety = db_safety
        self.in_mem_uuid_to_dbid = in_mem_uuid_to_dbid
        self.in_mem_dbid_to_uuid = in_mem_dbid_to_uuid
        self.db_timeout = mod_db_timeout

        # establish connection and initialize main cursor
        self.establish_db_conn()
        self.set_cursor()

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

        if not self.init_db_tables:
            # load requested in memory indices
            self.load_in_mem_chrm()
            self.load_in_mem_mod()
            self.load_in_mem_uuid()

    def establish_db_conn(self):
        if self.read_only:
            if not os.path.exists(self.fn):
                LOGGER.error((
                    'Modified base per-read database file ({}) does ' +
                    'not exist.').format(self.fn))
                raise mh.MegaError('Invalid mods DB filename.')
            self.db = sqlite3.connect('file:' + self.fn + '?mode=ro', uri=True)
            # use memory mapped file access
            self.db.execute('PRAGMA mmap_size = {}'.format(
                mh.MEMORY_MAP_LIMIT))
        else:
            self.db = sqlite3.connect(self.fn, timeout=self.db_timeout)
            self.db.execute('PRAGMA page_size = {}'.format(
                mh.SQLITE_PAGE_SIZE))
            self.db.execute('PRAGMA max_page_count = {}'.format(
                mh.SQLITE_MAX_PAGE_COUNT))
            self.db.execute("PRAGMA temp_store_directory = '{}'".format(
                os.path.dirname(self.fn)))
            if self.db_safety < 2:
                # set asynchronous mode to off for max speed
                self.db.execute('PRAGMA synchronous = OFF')
            if self.db_safety < 1:
                # set no rollback mode
                self.db.execute('PRAGMA journal_mode = OFF')
            self.db.execute('PRAGMA cache_size = {}'.format(
                -mh.SQLITE_CACHE_SIZE))
            self.db.execute('PRAGMA threads = {}'.format(
                mh.SQLITE_THREADS))

    def set_cursor(self):
        self.cur = self.db.cursor()

    def commit(self):
        self.db.commit()

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

    ##################################
    # load in-memory index functions #
    ##################################

    def _load_chrm_offsets(self, ref_names_and_lens):
        """ Positions are encoded as integer values ordered first by chromomes,
        then position, then strand. Forward strand (encoded as `1` as in
        `mappy`) positions are all even values and reverse strand (encoded as
        `-1`) are all odd values. Chromosomes are ordered as passed in, but
        expected in LooseVersion (as in unix `sort -V`) provided via
        megalodon_helper.RefName.
        """
        # store dict from chrm to chrm_len
        self._chrm_len_lookup = dict(
            (chrm_name, int(chrm_len))
            for chrm_name, chrm_len in zip(*ref_names_and_lens))
        self.chrm_names = ref_names_and_lens[0]
        self.chrm_lens = [int(cl) for cl in ref_names_and_lens[1]]

        self.num_chrms = len(self.chrm_names)
        # create data structures to convert (chrm, strand, pos) <--> pos_dbid
        # Store offsets to start of each chromosome
        self._chrm_offsets = np.insert(
            np.cumsum(self.chrm_lens[:-1]) * 2, 0, 0)
        # dictionary from chrm to dbid offsets
        self._chrm_offset_lookup = dict(zip(
            self.chrm_names, map(int, self._chrm_offsets)))

    def load_in_mem_chrm(self):
        ref_names_and_lens = list(zip(*self.cur.execute(
            'SELECT chrm, chrm_len FROM chrm').fetchall()))
        self._load_chrm_offsets(ref_names_and_lens)

    def _create_mod_lookups(self, mod_db_info):
        self.dbid_to_mod = dict(
            (dbid, mod_base) for dbid, mod_base, _, _ in mod_db_info)
        self.mod_to_dbid = dict(
            (mod_base, dbid) for dbid, mod_base, _, _ in mod_db_info)
        self.mod_to_can = dict(
            (mod_base, can_base) for _, mod_base, can_base, _ in mod_db_info)
        self.mod_to_long_name = dict(
            (mod_base, mln) for _, mod_base, _, mln in mod_db_info)

    def load_in_mem_mod(self):
        mod_db_info = self.cur.execute(
            'SELECT mod_id, mod_base, can_base, mod_long_name ' +
            'FROM mod_long_names').fetchall()
        self._create_mod_lookups(mod_db_info)

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
    # getter data functions #
    #########################

    @staticmethod
    def _get_strand_offset(strand):
        if strand == 1:
            return 0
        elif strand == -1:
            return 1
        raise mh.MegaError(
            'Invalid strand passed to DB (expected 1 or -1): {}'.format(
                strand))

    def get_pos_dbid(self, chrm, strand, pos):
        """ Compute database ID for position
        """
        if self.init_db_tables:
            raise mh.MegaError(
                'Cannot extract position database ID from connection opened ' +
                'for initialization.')
        if pos >= self._chrm_len_lookup[chrm]:
            raise mh.MegaError((
                'Attempt to extract position past the end of a chromosome.' +
                ' {}:{}').format(chrm, pos))
        strand_offset = self._get_strand_offset(strand)
        return self._chrm_offset_lookup[chrm] + (pos * 2) + strand_offset

    def get_pos_dbids(self, r_uniq_pos, chrm, strand):
        """ Get position database IDs. If positions are not found in the
        database they will be inserted.

        Args:
            r_uniq_pos (list): Unique positions (int).
            chrm (str): Chromosome name.
            strand (int): +1 = forward strandl -1 = reverse strand

        Returns:
            List of position database IDs for each entry in r_uniq_pos.
        """
        if self.init_db_tables:
            raise mh.MegaError(
                'Cannot extract position database IDs from connection ' +
                'opened for initialization.')
        if max(r_uniq_pos) >= self._chrm_len_lookup[chrm]:
            raise mh.MegaError((
                'Attempt to extract position past the end of a chromosome.' +
                ' {}:{}').format(chrm, max(r_uniq_pos)))
        strand_offset = self._get_strand_offset(strand)
        cs_offset = self._chrm_offset_lookup[chrm]
        return [cs_offset + (pos * 2) + strand_offset for pos in r_uniq_pos]

    def get_pos(self, pos_dbid):
        """ Get position from database ID

        Args:
            pos_dbid (int): Position database ID

        Returns:
            Chromosome name (str), strand (int; +1=forward strand, -1=reverse
                strand), position (int; 0-based)
        """
        if self.init_db_tables:
            raise mh.MegaError(
                'Cannot extract position from connection opened for ' +
                'initialization.')
        chrm_idx = np.searchsorted(self._chrm_offsets, pos_dbid, 'right') - 1
        chrm = self.chrm_names[chrm_idx]
        strand = 1 if pos_dbid % 2 == 0 else -1
        # note conversion to int so that sqlite3 does not store np.int32 bytes
        pos = (pos_dbid - int(self._chrm_offsets[chrm_idx])) // 2
        return chrm, strand, pos

    def get_mod_base(self, mod_dbid):
        """ Get modified base from database ID.

        Args:
            mod_dbid (int): Modified base database ID

        Returns:
            Modified base single letter code (str)
        """
        if self.init_db_tables:
            raise mh.MegaError(
                'Cannot extract modified base from connection opened for ' +
                'initialization.')
        try:
            return self.dbid_to_mod[mod_dbid]
        except KeyError:
            raise mh.MegaError(
                'Modified base ID not found in database: {}'.format(mod_dbid))

    def get_mod_base_dbid(self, mod_base):
        if self.init_db_tables:
            raise mh.MegaError(
                'Cannot extract modified base database ID from connection ' +
                'opened for initialization.')
        try:
            return self.mod_to_dbid[mod_base]
        except KeyError:
            raise mh.MegaError(
                'Modified base not found in database: {}'.format(mod_base))

    def get_alphabet_info(self):
        mod_base_to_can = dict()
        for mod_base, can_base, mln in self.get_full_mod_data():
            mod_base_to_can[mod_base] = (can_base, mln)

        # extract canonical bases associated with modified base
        can_bases = set(can_base for can_base, _ in mod_base_to_can.values())
        # determine first valid canonical alphabet compatible with database
        can_alphabet = None
        for v_alphabet in mh.VALID_ALPHABETS:
            if len(can_bases.difference(v_alphabet)) == 0:
                can_alphabet = v_alphabet
                break
        if can_alphabet is None:
            LOGGER.error(
                ('Mods database does not contain valid canonical ' +
                 'bases ({})').format(''.join(sorted(can_bases))))
            raise mh.MegaError('Invalid alphabet.')

        # compute full output alphabet and ordered modified base long names
        can_base_to_mods = dict(
            (can_base, [
                (mod_base, mln)
                for mod_base, (mcan_base, mln) in mod_base_to_can.items()
                if mcan_base == can_base])
            for can_base in can_alphabet)
        alphabet, collapse_alphabet = '', ''
        mod_long_names = []
        for can_base in can_alphabet:
            alphabet += can_base
            collapse_alphabet += can_base
            for mod_base, mln in can_base_to_mods[can_base]:
                alphabet += mod_base
                collapse_alphabet += can_base
                mod_long_names.append(mln)

        return alphabet, collapse_alphabet, mod_long_names

    def get_uuid(self, read_dbid):
        """ Get read UUID from database.

        Args:
            read_dbid (int): Database read UUID ID

        Returns:
            UUID (int): Universal read identifier
        """
        try:
            if self.in_mem_dbid_to_uuid:
                uuid = self.dbid_to_uuid[read_dbid]
            else:
                uuid = self.cur.execute(
                    'SELECT uuid FROM read WHERE read_id=?',
                    (read_dbid, )).fetchall()[0][0]
        except (TypeError, KeyError):
            raise mh.MegaError(
                'Read ID not found in database: {}'.format(read_dbid))
        return uuid

    def get_read_dbid(self, uuid):
        """ Get database ID given read UUID

        Args:
            uuid (int): Universal read identifier

        Returns:
            read_dbid (int): Database read UUID ID
        """
        try:
            if self.in_mem_uuid_to_dbid:
                read_dbid = self.uuid_to_dbid[uuid]
            else:
                read_dbid = self.cur.execute(
                    'SELECT read_id FROM read WHERE uuid=?',
                    (uuid, )).fetchall()[0][0]
        except (TypeError, KeyError):
            raise mh.MegaError(
                'UUID not found in database: {}'.format(uuid))
        return read_dbid

    def get_pos_stats(
            self, pos_dbid, return_uuids=False, get_without_index=False):
        """ Get all statistics mapped to a reference position. The data
        covering index should be created to increase the access speed for
        this function.

        Args:
            pos_dbid (int): Positon database ID (int)
            return_uuids (bool): Whether to return database read ids (defalt)
                or UUIDs.
            get_without_index (bool): Force extraction without index. This
                requires iterating through the entire data table.

        Returns:
            List containing megalodon.mods.ModsDb.mod_data objects mapped to
                specified reference position.
        """
        if (not get_without_index) and self.check_data_covering_index_exists():
            raise mh.MegaError(
                'Cannot extract position statistics without covering index.')
        read_id_conv = self.get_uuid if return_uuids else lambda x: x
        self.cur.execute(
            'SELECT score_read, score, score_mod FROM data WHERE score_pos=?',
            (pos_dbid, ))
        chrm, strand, pos = self.get_pos(pos_dbid)
        return [
            self.mod_data(read_id_conv(read_dbid), chrm, strand, pos, score,
                          self.get_mod_base(mod_dbid))
            for read_dbid, score, mod_dbid in self.cur]

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
            'SELECT MAX(mod_id) FROM mod_long_names').fetchone()[0]
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

    def get_num_uniq_stats(self):
        """ Get number of per-read modified base statistics stored in database
        """
        num_stats = self.cur.execute(
            'SELECT MAX(rowid) FROM data').fetchone()[0]
        if num_stats is None:
            return 0
        return num_stats

    #########################
    # insert data functions #
    #########################

    def insert_chrms(self, ref_names_and_lens):
        # chromosomes must be entered only once for a database
        if len(self.cur.execute('SELECT * FROM chrm').fetchall()) > 0:
            raise mh.MegaError(
                'Chromosomes/contigs have already been set for this database.')
        # use version sort for chromosomes/contigs
        s_ref_names_and_lens = list(zip(*sorted(
            list(zip(*ref_names_and_lens)),
            key=lambda x: mh.RefName(x[0]))))
        # save chrms to database
        self.cur.executemany('INSERT INTO chrm (chrm, chrm_len) VALUES (?,?)',
                             zip(*s_ref_names_and_lens))
        # add save to internal data structure to determine database ids for
        # positions
        self._load_chrm_offsets(s_ref_names_and_lens)

    def insert_mod_long_names(self, mod_long_names, mod_base_to_can):
        # modified bases must be entered only once for a database
        if len(self.cur.execute(
                'SELECT * FROM mod_long_names').fetchall()) > 0:
            raise mh.MegaError(
                'Modified bases have already been set for this database.')
        insert_mod_info = [(mod_base, mod_base_to_can[mod_base], mln)
                           for mod_base, mln in mod_long_names]
        self.cur.executemany(
            'INSERT INTO mod_long_names (mod_base, can_base, mod_long_name) ' +
            'VALUES (?,?,?)', insert_mod_info)
        self._create_mod_lookups([
            (mod_dbid, *mb_info)
            for mod_dbid, mb_info in enumerate(insert_mod_info)])

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
            if self.in_mem_uuid_to_dbid:
                self.dbid_to_uuid[read_dbid] = uuid
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
            # self.uuid_to_dbid updated automatically by assignment above

        return [uuid_to_dbid[uuid] for uuid in uuids]

    def insert_uuid(self, uuid):
        """ Insert unique read identifier into database.

        Args:
            uuid (str): Unique read identifier

        Returns:
            Database ID.
        """
        self.cur.execute('INSERT INTO read (uuid) VALUES (?)', (uuid,))
        read_dbid = self.cur.lastrowid
        if self.in_mem_dbid_to_uuid:
            self.uuid_to_dbid[uuid] = read_dbid
        if self.in_mem_uuid_to_dbid:
            self.dbid_to_uuid[read_dbid] = uuid
        return read_dbid

    def insert_uuids(self, uuids):
        """ Insert unique read identifiers into database.

        Args:
            uuids (list): List of unique read identifiers

        Returns:
            Database IDs.
        """
        uuids = list(uuids)
        next_read_dbid = self.get_num_uniq_reads() + 1
        self.cur.executemany(
            'INSERT INTO read (uuid) VALUES (?)', ((uuid, ) for uuid in uuids))
        read_dbids = list(range(next_read_dbid, next_read_dbid + len(uuids)))
        if self.in_mem_dbid_to_uuid:
            self.dbid_to_uuid.update(zip(read_dbids, uuids))
        if self.in_mem_uuid_to_dbid:
            self.uuid_to_dbid.update(zip(uuids, read_dbids))
        return read_dbids

    def insert_read_data(self, r_insert_data, read_dbid):
        """ Insert data from a single read.

        Args:
            r_insert_data (list): Each entry should contain 4 elements:
                1) modified base log probability (float), 2) position database
                ID, and 3) modified base database ID.
            read_dbid (int): Read database ID
        """
        self.cur.executemany(
            'INSERT INTO data VALUES (?,?,?,?)',
            ((*score_pos_mod, read_dbid) for score_pos_mod in r_insert_data))

    def insert_batch_data(self, b_insert_data):
        """ Insert batch data

        Args:
            b_insert_data (list): Each entry should contain 4 elements:
                1) modified base log probability (float), 2) position database
                ID 3) modified base database ID, and 4) read database ID.
        """
        self.cur.executemany('INSERT INTO data VALUES (?,?,?,?)',
                             b_insert_data)

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

    def create_mod_index(self):
        self.cur.execute('CREATE UNIQUE INDEX mod_idx ON ' +
                         'mod_long_names(mod_base)')

    def create_data_covering_index(self):
        # TODO add progress info to this step.
        self.cur.execute('CREATE INDEX data_cov_idx ON data(' +
                         'score_pos, score_read, score_mod, score)')

    ##################
    # data iterators #
    ##################

    def iter_pos_scores(
            self, convert_pos=False, compute_llrs=False, pos_range=None):
        """ Iterate over scores grouped by position. Default arguments iterate
        over raw database values for maximal speed. Yeilds a tuple of
        1) position and 2) statistics

        If convert_pos it False (default), yield position database id, if
        convert_pos is True yield (chrm, strand, pos)

        If compute_llrs is False (default), yield lists containing tuples of
        1) read_id, 2) mod database ID and 3) log probability
        If compute_llrs is True, yield dictionary of mod_base keys to list log
        likelihood ratios.

        Note this function iterates over the index created by
        create_data_covering_index, so should be very fast.

        If pos_range is provided, the database query will restrict the
        extracted sites to a specific range. This parameter should consist of
        the contig name (str), start (int) and end (int) coordinates.
        """
        def extract_pos_llrs(pos_lps):
            mod_llrs = dict((self.get_mod_base(mod_dbid), [])
                            for mod_dbid in set(list(zip(*pos_lps))[1]))
            prev_dbid = None
            mod_bs, r_lps = [], []
            for read_dbid, mod_dbid, lp in sorted(pos_lps):
                if prev_dbid != read_dbid and prev_dbid is not None:
                    # compute and store log likelihood ratios
                    with np.errstate(divide='ignore'):
                        can_lp = np.log1p(-np.exp(r_lps).sum())
                    for mod_b, r_lp in zip(mod_bs, r_lps):
                        mod_llrs[mod_b].append(can_lp - r_lp)
                    mod_bs, r_lps = [], []
                prev_dbid = read_dbid
                mod_bs.append(self.get_mod_base(mod_dbid))
                r_lps.append(lp)
            # compute and store last log likelihood ratios
            with np.errstate(divide='ignore'):
                can_lp = np.log1p(-np.exp(np.array(r_lps)).sum())
            for mod_b, r_lp in zip(mod_bs, r_lps):
                mod_llrs[mod_b].append(can_lp - r_lp)

            return mod_llrs

        self.check_data_covering_index_exists()

        pos_func = self.get_pos if convert_pos else lambda x: x
        stat_func = extract_pos_llrs if compute_llrs else lambda x: x

        pos_lps = list()
        # use local cursor since extracting pos or mod might use class cursor
        local_cursor = self.db.cursor()
        if pos_range is None:
            local_cursor.execute(
                'SELECT score_pos, score_mod, score_read, score FROM data ' +
                'ORDER BY score_pos')
        else:
            # determine pos database ID range
            chrm, pos_st, pos_en = pos_range
            if pos_st < 0:
                pos_st = 0
            if pos_en >= self._chrm_len_lookup[chrm]:
                pos_en = self._chrm_len_lookup[chrm] - 1
            pos_dbid_range = (self.get_pos_dbid(chrm, 1, pos_st),
                              self.get_pos_dbid(chrm, -1, pos_en))
            # restrict query to specified range
            local_cursor.execute(
                'SELECT score_pos, score_mod, score_read, score FROM data ' +
                'WHERE score_pos BETWEEN ? AND ? ORDER BY score_pos',
                pos_dbid_range)
        # initialize variables with first value
        first_score = local_cursor.fetchone()
        # if no scores are stored break out of iterator
        if first_score is None:
            return
        prev_pos, mod_dbid, read_dbid, lp = first_score
        pos_lps.append((read_dbid, mod_dbid, lp))
        for curr_pos, mod_dbid, read_dbid, lp in local_cursor:
            if curr_pos != prev_pos:
                yield pos_func(prev_pos), stat_func(pos_lps)
                pos_lps = list()
            pos_lps.append((read_dbid, mod_dbid, lp))
            prev_pos = curr_pos
        yield pos_func(prev_pos), stat_func(pos_lps)

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
        for chrm_dbid, chrm, chrm_len in zip(
                range(self.num_chrms), self.chrm_names, self.chrm_lens):
            yield chrm_dbid, chrm, chrm_len

    def iter_mod_bases(self):
        """ Iterate over modified base information from database

        Yields:
            Tuple containing 1) database ID, and 2) modified base single letter
                code (str)
        """
        if self.init_db_tables:
            raise mh.MegaError(
                'Cannot iterate modified bases from connection ' +
                'opened for initialization.')
        for mod_dbid, mod_base in self.dbid_to_mod.items():
            yield mod_dbid, mod_base

    def get_mod_long_names(self):
        """ Get modified base long names

        Returns:
            List of tuples containing 1) modified base single letter code and
                2) modified base long name
        """
        return self.cur.execute(
            'SELECT mod_base, mod_long_name FROM mod_long_names').fetchall()

    def get_full_mod_data(self):
        """ Get all modified base data

        Returns:
            List of tuples containing 1) modified base single letter code,
                2) canonical base and 3) modified base long name
        """
        return self.cur.execute(
            'SELECT mod_base, can_base, mod_long_name ' +
            'FROM mod_long_names').fetchall()

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

    def iter_data(self):
        # use local cursor since other processing might use class cursor
        local_cursor = self.db.cursor()
        local_cursor.execute(
            'SELECT score, uuid, mod_base, score_pos FROM data ' +
            'INNER JOIN read ON data.score_read = read.read_id ' +
            'INNER JOIN mod_long_names ' +
            'ON data.score_mod = mod_long_names.mod_id')
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
    for _, mods_pos_llrs in mods_db.iter_pos_scores(compute_llrs=True):
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
    mods_db = ModsDb(mods_db_fn)
    for (chrm, strand, pos), mods_pos_llrs in mods_db.iter_pos_scores(
            convert_pos=True, compute_llrs=True):
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

def annotate_all_mods(
        r_start, ref_seq, r_mod_scores, strand, mods_info,
        per_site_thresh=None):
    """ Annotate reference sequence with called modified bases.

    Args:
        r_start (int): Reference start position for this read.
        ref_seq (str): Read-centric reference sequence corresponding to
            this read.
        r_mod_scores (list): Per-reference position read modified base calls
            including postiion (on chrm/strand), modbase log probs and modified
            base single letter codes
        strand (int): 1 for forward strand -1 for reverse strand
        mods_info (mods.ModInfo): Object containing information about modified
            base processing
        per_site_thresh (np.ndarray): Score thresholds per position,
            score := log(P_can/P_mod).

    Returns:
        mods.ANNOT_MODS object annotated with all modified bases.

    Note: Reference sequence is in read orientation and mod calls are in
    genome coordiates.
    """
    all_mods_seq, all_mods_qual = [], []
    prev_pos = 0
    if strand == -1:
        ref_seq = ref_seq[::-1]
    for mod_pos, mod_lps, mod_bases in sorted(r_mod_scores):
        if mod_lps is None:
            base_lp = MOD_MAP_INVALID_BASE_LP
            base = ref_seq[mod_pos - r_start]
        else:
            can_lp = np.log1p(-np.exp(mod_lps).sum())
            mod_thresh = (mods_info.mod_thresh if per_site_thresh is None else
                          per_site_thresh[mod_pos])
            if can_lp - mod_lps.max() > mod_thresh:
                base_lp = can_lp
                base = ref_seq[mod_pos - r_start]
            else:
                most_prob_mod = np.argmax(mod_lps)
                base_lp = mod_lps[most_prob_mod]
                base = mod_bases[most_prob_mod]
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


def annotate_mods_per_mod(
        r_start, ref_seq, r_mod_scores, strand, mods_info,
        per_site_thresh=None):
    """ Annotate reference sequence with called modified bases. Produce one
    mods.ANNOT_MODS output for each modified base in mods_info.mod_long_names.

    Args:
        r_start (int): Reference start position for this read.
        ref_seq (str): Read-centric reference sequence corresponding to
            this read.
        r_mod_scores (list): Per-reference position read modified base calls
            including postiion (on chrm/strand), modbase log probs and modified
            base single letter codes
        strand (int): 1 for forward strand -1 for reverse strand
        mods_info (mods.ModInfo): Object containing information about modified
            base processing
        per_site_thresh (np.ndarray): Score thresholds per position,
            score := log(P_can/P_mod).

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
    for mod_pos, mod_lps, mod_bases in sorted(r_mod_scores):
        if mod_lps is not None:
            can_lp = np.log1p(-np.exp(mod_lps).sum())
        # annotate per-mod sequences and qualities
        for mod_idx, mod_base in enumerate(mod_bases):
            if mod_lps is None:
                base_lp = MOD_MAP_INVALID_BASE_LP
                base = ref_seq[mod_pos - r_start]
            else:
                mod_thresh = (mods_info.mod_thresh if per_site_thresh is None
                              else per_site_thresh[mod_pos])
                # called canonical
                if can_lp - mod_lps[mod_idx] > mod_thresh:
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


def format_mm_ml_tags(r_start, ref_seq, r_mod_scores, strand, mods_info):
    """ Format Mm and Ml tags for BAM output. See
    https://github.com/samtools/hts-specs/pull/418 for format details.

    Args:
        r_start (int): Reference start position for this read.
        ref_seq (str): Read-centric reference sequence corresponding to
            this read.
        r_mod_scores (list): Per-reference position read modified base calls
            including postiion (on chrm/strand), modbase log probs and modified
            base single letter codes
        strand (int): 1 for forward strand -1 for reverse strand
        mods_info (mods.ModInfo): Object containing information about modified
            base processing

    Returns:
        Mm string tag and Ml array tag
    """

    per_mod_probs = defaultdict(list)
    for mod_pos, mod_lps, mod_bases in sorted(r_mod_scores):
        # mod_lps is set to None if invalid sequence is encountered or too
        # few events are found around a mod
        if mod_lps is None:
            continue
        for mod_lp, mod_base in zip(mod_lps, mod_bases):
            mod_prob = np.exp(mod_lp)
            if mod_prob < mods_info.map_min_prob:
                continue
            read_pos = mod_pos - r_start
            if strand == -1:
                read_pos = len(ref_seq) - 1 - read_pos
            per_mod_probs[mod_base].append((read_pos, mod_prob))

    mm_tag, ml_tag = '', array.array('B')
    for mod_base, pos_probs in per_mod_probs.items():
        mod_poss, probs = zip(*sorted(pos_probs))
        can_base = mods_info.mod_base_to_can[mod_base]
        # compute modified base positions relative to the running total of the
        # associated canonical base
        can_base_mod_poss = np.cumsum(
            [1 if b == can_base else 0 for b in ref_seq])[
                np.array(mod_poss)] - 1
        mm_tag += '{}+{}{};'.format(
            can_base, mh.convert_legacy_mods(mod_base),
            ''.join(',{}'.format(d) for d in np.diff(np.insert(
                can_base_mod_poss, 0, -1)) - 1))
        # extract mod scores and scale to 0-255 range
        scaled_probs = np.floor(np.array(probs) * 256)
        # last interval includes prob=1
        scaled_probs[scaled_probs == 256] = 255
        ml_tag.extend(scaled_probs.astype(np.uint8))

    return mm_tag, ml_tag


def get_mod_annotated_seqs(
        mods_info, sig_map_res, r_ref_pos, r_mod_scores, r_ref_seq):
    all_mods_seq = per_mod_seqs = per_site_thresh = None
    if (mods_info.do_ann_all_mods or mods_info.do_output.mod_map) and \
       sig_map_res is not None and \
       sig_map_res.ref_out_info.per_site_threshs is not None:
        try:
            per_site_thresh = sig_map_res.ref_out_info.per_site_threshs[(
                r_ref_pos.chrm, r_ref_pos.strand)]
        except KeyError:
            cov_pos = '' if len(r_mod_scores) == 0 else ','.join(
                map(str, sorted(list(zip(*r_mod_scores))[0])))
            LOGGER.debug((
                '{} PerSiteThreshContigNotFound {}:{} {}').format(
                    sig_map_res.read_id, r_ref_pos.chrm,
                    '+' if r_ref_pos.strand == 1 else '-', cov_pos))
    try:
        if mods_info.do_ann_all_mods:
            # ignore divide around full annotate_mods call to avoid overhead
            # on many calls to errstate
            with np.errstate(divide='ignore'):
                all_mods_seq = annotate_all_mods(
                    r_ref_pos.start, r_ref_seq, r_mod_scores,
                    r_ref_pos.strand, mods_info, per_site_thresh)
        if mods_info.do_output.mod_map:
            if mods_info.map_emulate_bisulfite:
                with np.errstate(divide='ignore'):
                    per_mod_seqs = annotate_mods_per_mod(
                        r_ref_pos.start, r_ref_seq, r_mod_scores,
                        r_ref_pos.strand, mods_info, per_site_thresh)
            else:
                per_mod_seqs = format_mm_ml_tags(
                    r_ref_pos.start, r_ref_seq, r_mod_scores,
                    r_ref_pos.strand, mods_info)
    except KeyError:
        cov_pos = '' if len(r_mod_scores) == 0 else ','.join(
            map(str, sorted(list(zip(*r_mod_scores))[0])))
        LOGGER.debug((
            '{} PerSiteThreshSiteNotFound {}:{} {}').format(
                sig_map_res.read_id, r_ref_pos.chrm,
                '+' if r_ref_pos.strand == 1 else '-', cov_pos))

    return all_mods_seq, per_mod_seqs


def send_signal_mapping(
        mod_sig_map_q, sig_map_res, all_mods_seq, failed_reads_q, fast5_fn):
    # send mod annotated seqs to signal mapping queue if requested
    if mod_sig_map_q is not None and sig_map_res.pass_filts:
        is_valid_mapping = True
        # import locally so that import of mods module does not require
        # taiyaki install (required for signal_mapping module)
        from megalodon import signal_mapping
        if sig_map_res.ref_out_info.do_output.mod_sig_maps:
            invalid_chars = set(all_mods_seq.mod_seq).difference(
                sig_map_res.ref_out_info.alphabet_info.alphabet)
            if len(invalid_chars) > 0:
                is_valid_mapping = False
                if failed_reads_q is not None:
                    fail_msg = (
                        'Invalid charcters for signal mapping found in ' +
                        'mapped sequence: ({})').format(''.join(invalid_chars))
                    # Send invalid character code to failed reads queue
                    failed_reads_q.put(tuple(mh.READ_STATUS(
                        is_err=True, do_update_prog=True, err_type=fail_msg,
                        fast5_fn=fast5_fn)))
            else:
                # replace reference sequence with mod annotated sequence
                sig_map_res = sig_map_res._replace(
                    ref_seq=all_mods_seq.mod_seq)

        if is_valid_mapping:
            mod_sig_map_q.put(signal_mapping.get_remapping(*sig_map_res[1:]))


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
        r_ref_pos, r_ref_seq, ref_to_block, r_post, mods_info, mod_sig_map_q,
        sig_map_res, signal_reversed, uuid, failed_reads_q, fast5_fn):
    def iter_motif_sites():
        search_ref_seq = r_ref_seq[::-1] if signal_reversed else r_ref_seq
        for motif, rel_pos, mod_bases, raw_motif in mods_info.all_mod_motifs:
            for motif_match in motif.finditer(search_ref_seq):
                m_pos = motif_match.start() + rel_pos
                if signal_reversed:
                    m_pos = len(r_ref_seq) - m_pos - 1
                yield m_pos, mod_bases, rel_pos, raw_motif

    def filter_mod_score(r_mod_scores):
        # remove uncalled sites and sites too close to the edge of a read
        return [pms for pms in r_mod_scores
                if pms[1] is not None and
                r_ref_pos.start + mods_info.edge_buffer < pms[0] <
                r_ref_pos.end - mods_info.edge_buffer]

    # call all mods overlapping this read
    r_mod_scores = []
    # ignore when one or more mod_llrs is -inf (or close enough for exp)
    # occurs in compute_log_probs function, but more efficient to seterr
    # at this higher level
    with np.errstate(divide='ignore', over='ignore'):
        for pos, mod_bases, rel_pos, raw_motif in iter_motif_sites():
            if (r_ref_pos.strand == 1 and not signal_reversed) or (
                    r_ref_pos.strand == -1 and signal_reversed):
                mod_ref_pos = r_ref_pos.start + pos
            else:
                mod_ref_pos = r_ref_pos.end - pos - 1
            pos_bb, pos_ab = min(mods_info.mod_context_bases, pos), min(
                mods_info.mod_context_bases, len(r_ref_seq) - pos - 1)
            try:
                pos_ref_seq = mh.seq_to_int(
                    r_ref_seq[pos - pos_bb:pos + pos_ab + 1])
            except mh.MegaError:
                # Add None score for per-read annotation (to be filtered)
                r_mod_scores.append((mod_ref_pos, None, mod_bases))
                LOGGER.debug('InvalidSequence {}:{}'.format(
                    r_ref_pos.chrm, mod_ref_pos))
                continue

            pos_can_mods = np.zeros_like(pos_ref_seq)
            blk_start, blk_end = (ref_to_block[pos - pos_bb],
                                  ref_to_block[pos + pos_ab])
            if blk_end - blk_start < (mods_info.mod_context_bases * 2) + 1:
                # need as many "events/strides" as bases for valid mapping
                # Add None scores for per-read annotation (to be filtered)
                r_mod_scores.append((mod_ref_pos, None, mod_bases))
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
            r_mod_scores.append((mod_ref_pos, loc_mod_lps, mod_bases))

    all_mods_seq, per_mod_seqs = get_mod_annotated_seqs(
        mods_info, sig_map_res, r_ref_pos, r_mod_scores, r_ref_seq)
    r_mod_scores = filter_mod_score(r_mod_scores)
    send_signal_mapping(
        mod_sig_map_q, sig_map_res, all_mods_seq, failed_reads_q, fast5_fn)
    if mods_info.do_output.db:
        r_insert_data = [
            (mod_lp, mods_info.get_pos_dbid(
                r_ref_pos.chrm, r_ref_pos.strand, mod_pos),
             mods_info.get_mod_base_dbid(mod_base))
            for mod_pos, mod_lps, mod_bases in r_mod_scores
            for mod_lp, mod_base in zip(mod_lps, mod_bases)]
    else:
        r_insert_data = None

    mod_out_text = None
    if mods_info.do_output.text:
        str_strand = mh.int_strand_to_str(r_ref_pos.strand)
        txt_tmplt = '\t'.join('{}' for _ in ModsDb.text_field_names) + '\n'
        mod_out_text = []
        for pos, mod_lps, mod_bases in r_mod_scores:
            with np.errstate(divide='ignore'):
                can_lp = np.log1p(-np.exp(mod_lps).sum())
            for mod_lp, mod_base in zip(mod_lps, mod_bases):
                mod_out_text.append(txt_tmplt.format(
                    uuid, r_ref_pos.chrm, str_strand, pos, mod_lp, can_lp,
                    mod_base))
        mod_out_text = ''.join(mod_out_text)

    return r_insert_data, all_mods_seq, per_mod_seqs, mod_out_text


#######################
# Per-read Mod Output #
#######################

def init_mods_db(mods_info, ref_names_and_lens):
    """ Initialize a new modified bases database.

    Args:
        mods_info (mods.ModInfo): ModInfo object minimally containing
            mods_db_fn, mod_long_names and mod_base_to_can attributes.
        db_safety (int): Database safety value as defined in mods.ModsDb
        ref_names_and_lens (list): List containing two lists of the same length
            representing chromosome 1) names and 2) lengths.
    """
    mods_db = ModsDb(
        mods_info.mods_db_fn, db_safety=mods_info.db_safety, read_only=False,
        init_db_tables=True)
    mods_db.insert_chrms(ref_names_and_lens)
    mods_db.insert_mod_long_names(
        mods_info.mod_long_names, mods_info.mod_base_to_can)
    mods_db.close()


def _get_mods_queue(
        mods_q, mods_conn, mods_info, map_info, ref_out_info, aux_failed_q):
    def write_mod_alignment(
            read_id, mod_seq, mod_quals, chrm, strand, r_st, fp):
        # convert to reference based sequence
        if strand == -1:
            mod_seq = mh.revcomp(mod_seq)
            mod_quals = mod_quals[::-1]
        a = mapping.prepare_mapping(
            read_id, mod_seq, flag=0 if strand == 1 else 16,
            ref_id=fp.get_tid(chrm), ref_st=r_st,
            qual=array.array('B', mod_quals), map_qual=MOD_MAP_MAX_QUAL,
            tags=[('RG', MOD_MAP_RG_ID)])
        fp.write(a)

    def write_alignment_w_tags(
            read_id, ref_seq, chrm, strand, r_st, fp, mods_scores):
        # convert to reference based sequence
        if strand == -1:
            ref_seq = mh.revcomp(ref_seq)
        a = mapping.prepare_mapping(
            read_id, ref_seq, flag=0 if strand == 1 else 16,
            ref_id=fp.get_tid(chrm), ref_st=r_st,
            map_qual=MOD_MAP_MAX_QUAL,
            tags=[('RG', MOD_MAP_RG_ID)], mods_scores=mods_scores)
        fp.write(a)

    def store_mod_call(mod_res, been_warned_timeout, been_warned_other):
        (r_insert_data, all_mods_seq, per_mod_seqs, mod_out_text), (
            read_id, chrm, strand, r_start, ref_seq, read_len, q_st, q_en,
            cigar) = mod_res
        if mods_info.do_output.db:
            try:
                read_dbid = mods_db.insert_uuid(read_id)
                data_commited = False
                while not data_commited:
                    try:
                        mods_db.insert_read_data(r_insert_data, read_dbid)
                        mods_db.commit()
                        data_commited = True
                    except sqlite3.OperationalError as e:
                        if not been_warned_timeout:
                            LOGGER.warning(
                                DB_TIMEOUT_ERR_MSG.format('data', str(e)))
                            been_warned_timeout = True
                        LOGGER.debug('ModDBTimeout {}'.format(str(e)))
                        mods_db.establish_db_conn()
                        mods_db.set_cursor()
            except Exception as e:
                if not been_warned_other:
                    LOGGER.warning(
                        'Error inserting modified base scores into DB. ' +
                        'See log debug output for error details.')
                    been_warned_other = True
                LOGGER.debug('ModDBInsertError {}\n{}'.format(
                    str(e), traceback.format_exc()))

        if mods_info.do_output.text and len(mod_out_text) > 0:
            mods_txt_fp.write(mod_out_text)
        if ref_out_info.do_output.mod_pr_refs and mapping.read_passes_filters(
                ref_out_info.filt_params, read_len, q_st, q_en, cigar):
            pr_refs_fp.write('>{}\n{}\n'.format(read_id, all_mods_seq.mod_seq))
        if mods_info.do_output.mod_map:
            if mods_info.map_emulate_bisulfite:
                for mod_base, _ in mods_info.mod_long_names:
                    mod_seq, mod_qual = per_mod_seqs[mod_base]
                    write_mod_alignment(
                        read_id, mod_seq, mod_qual, chrm, strand, r_start,
                        mod_map_fps[mod_base])
            else:
                write_alignment_w_tags(
                    read_id, ref_seq, chrm, strand, r_start,
                    mod_map_fp, per_mod_seqs)

        return been_warned_timeout, been_warned_other

    try:
        LOGGER.debug('GetterStarting')
        if mods_info.do_output.db:
            mods_db = ModsDb(
                mods_info.mods_db_fn, db_safety=mods_info.db_safety,
                read_only=False, mod_db_timeout=mods_info.mod_db_timeout)
        if mods_info.do_output.text:
            mods_txt_fp = open(mh.get_megalodon_fn(
                mods_info.out_dir, mh.PR_MOD_TXT_NAME), 'w')
            mods_txt_fp.write('\t'.join(ModsDb.text_field_names) + '\n')
        if ref_out_info.do_output.mod_pr_refs:
            pr_refs_fp = open(mh.get_megalodon_fn(
                mods_info.out_dir, mh.PR_REF_NAME), 'w')
        if mods_info.do_output.mod_map:
            try:
                w_mode = mh.MAP_OUT_WRITE_MODES[map_info.map_fmt]
            except KeyError:
                raise mh.MegaError('Invalid mapping output format')
            header = {
                'HD': {'VN': '1.4'},
                'SQ': [{'LN': ref_len, 'SN': ref_name}
                       for ref_name, ref_len in sorted(
                    zip(*map_info.ref_names_and_lens))],
                'RG': [{'ID': MOD_MAP_RG_ID, 'SM': SAMPLE_NAME}, ]}
            mod_map_bn = mh.get_megalodon_fn(
                mods_info.out_dir, mh.MOD_MAP_NAME)
            if mods_info.map_emulate_bisulfite:
                mod_map_fns = [
                    (mod_base, '{}.{}.'.format(mod_map_bn, mln))
                    for mod_base, mln in mods_info.mod_long_names]
                mod_map_fps = dict((
                    (mod_base, pysam.AlignmentFile(
                        mod_map_fn + map_info.map_fmt, w_mode, header=header,
                        reference_filename=map_info.cram_ref_fn))
                    for mod_base, mod_map_fn in mod_map_fns))
            else:
                mod_map_fp = pysam.AlignmentFile(
                    '{}.{}'.format(mod_map_bn, map_info.map_fmt),
                    w_mode, header=header,
                    reference_filename=map_info.cram_ref_fn)
        been_warned_timeout = been_warned_other = False
        workers_active = True
        LOGGER.debug('GetterInitComplete')
    except Exception as e:
        aux_failed_q.put(('ModsInitError', str(e), traceback.format_exc()))
        return

    try:
        while workers_active or not mods_q.empty():
            try:
                mod_res = mods_q.get(timeout=0.1)
                r_val = mh.log_errors(
                    store_mod_call, mod_res, been_warned_timeout,
                    been_warned_other)
                if r_val is not None:
                    been_warned_timeout, been_warned_other = r_val
            except queue.Empty:
                if mods_conn.poll():
                    workers_active = False
        LOGGER.debug('GetterClosing')
    except Exception as e:
        aux_failed_q.put((
            'ModsProcessingError', str(e), traceback.format_exc()))
    finally:
        if mods_info.do_output.text:
            mods_txt_fp.close()
        if ref_out_info.do_output.mod_pr_refs:
            pr_refs_fp.close()
        if mods_info.do_output.mod_map:
            if mods_info.map_emulate_bisulfite:
                for mod_map_fp in mod_map_fps.values():
                    mod_map_fp.close()
            else:
                mod_map_fp.close()
        if mods_info.do_output.db:
            if not mods_info.skip_db_index:
                LOGGER.debug('CreatingIndex')
                mods_db.create_data_covering_index()
            LOGGER.debug('ClosingDB')
            mods_db.close()


if _PROFILE_MODS_QUEUE:
    _get_mods_queue_wrapper = _get_mods_queue

    def _get_mods_queue(*args):
        import cProfile
        cProfile.runctx('_get_mods_queue_wrapper(*args)', globals(), locals(),
                        filename='mods_getter_queue.prof')


############
# Mod Info #
############

class ModInfo:
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
                motif = mh.compile_motif_pat(raw_motif)
                self.all_mod_motifs.append((motif, pos, mod_bases, raw_motif))

            if not self.distinct_motifs():
                raise mh.MegaError(
                    'One provided motif can be found within another motif. ' +
                    'Only distinct sets of motifs are accepted')

    def __init__(
            self, model_info, all_mod_motifs_raw=None, mod_all_paths=False,
            mod_context_bases=None, mods_calib_fn=None,
            mod_output_fmts=[mh.MOD_BEDMETHYL_NAME],
            edge_buffer=mh.DEFAULT_EDGE_BUFFER, agg_info=DEFAULT_AGG_INFO,
            mod_thresh=0.0, do_ann_all_mods=False, map_emulate_bisulfite=False,
            map_base_conv=None, map_min_prob=mh.DEFAULT_MOD_MIN_PROB,
            mod_db_timeout=mh.DEFAULT_MOD_DATABASE_TIMEOUT, db_safety=0,
            out_dir=None, skip_db_index=False, do_output=None):
        # this is pretty hacky, but these attributes are stored here as
        # they are generally needed alongside other modbase info
        # don't want to pass all of these parameters around individually though
        # as this would make function signatures too complicated
        self.mod_all_paths = mod_all_paths
        self.mod_context_bases = mod_context_bases
        self.mod_long_names = model_info.mod_long_names
        self.calib_table = calibration.ModCalibrator(mods_calib_fn)
        self.mod_output_fmts = mod_output_fmts
        self.edge_buffer = edge_buffer
        self.agg_info = agg_info
        self.mod_thresh = mod_thresh
        # TODO move these attributes to do_output
        self.do_ann_all_mods = do_ann_all_mods
        self.map_emulate_bisulfite = map_emulate_bisulfite
        self.map_base_conv_raw = map_base_conv
        self.map_min_prob = map_min_prob
        self.mod_db_timeout = mod_db_timeout
        self.db_safety = db_safety
        self.out_dir = out_dir
        self.skip_db_index = skip_db_index
        self.do_output = do_output

        self.mods_db_arrays_added = False

        self.mods_db_fn = mh.get_megalodon_fn(self.out_dir, mh.PR_MOD_NAME)

        self.alphabet = model_info.can_alphabet
        self.ncan_base = len(self.alphabet)
        try:
            self.alphabet = self.alphabet.decode()
        except AttributeError:
            pass
        LOGGER.info(model_info.get_alphabet_str())

        self.nbase = len(self.alphabet)
        self.n_can_state = (self.ncan_base + self.ncan_base) * (
            self.ncan_base + 1)
        if model_info.is_cat_mod:
            self.nmod_base = model_info.n_mods
            self.can_base_mods = model_info.can_base_mods
            self.mod_base_to_can = model_info.mod_base_to_can
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

    def add_mods_db_arrays(self, mods_db):
        self.mod_to_dbid = mods_db.mod_to_dbid
        self._chrm_len_lookup = mods_db._chrm_len_lookup
        self._chrm_offset_lookup = mods_db._chrm_offset_lookup
        self.mods_db_arrays_added = True

    def get_pos_dbid(self, chrm, strand, pos):
        """ Compute database ID for position
        """
        if not self.mods_db_arrays_added:
            raise mh.MegaError(
                'Must run add_mods_db_arrays to extract pos_dbid.')
        if pos >= self._chrm_len_lookup[chrm]:
            raise mh.MegaError((
                'Attempt to extract position past the end of a chromosome.' +
                ' {}:{}').format(chrm, pos))
        strand_offset = ModsDb._get_strand_offset(strand)
        return self._chrm_offset_lookup[chrm] + (pos * 2) + strand_offset

    def get_mod_base_dbid(self, mod_base):
        if not self.mods_db_arrays_added:
            raise mh.MegaError(
                'Must run add_mods_db_arrays to extract mod_base_dbid.')
        try:
            return self.mod_to_dbid[mod_base]
        except KeyError:
            raise mh.MegaError(
                'Modified base not found in database: {}'.format(mod_base))


#################
# modVCF Writer #
#################

class ModSite:
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


class ModVcfWriter:
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


class ModBedMethylWriter:
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
                BEDMETHYL_TMPLT.format(
                    chrom=mod_site.chrom, pos=mod_site.pos,
                    end=mod_site.pos + 1, strand=mod_site.strand, cov=cov,
                    score=min(int(cov), 1000),
                    perc=np.around(mod_prop * 100, 1)))
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


class ModWigWriter:
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
                        '{} {}'.format(pos + 1, mod_prop)
                        for pos, mod_prop in sorted(cs_mod_sites))) + '\n')


##########################
# Mods Aggregation Class #
##########################

class AggMods(mh.AbstractAggregationClass):
    """ Class to assist in database queries for per-site aggregation of
    modified base calls over reads.
    """

    def __init__(self, mods_db_fn, agg_info=DEFAULT_AGG_INFO,
                 write_mod_lp=False, load_uuid_index_in_memory=False):
        if load_uuid_index_in_memory:
            self.mods_db = ModsDb(mods_db_fn, in_mem_dbid_to_uuid=True)
        else:
            self.mods_db = ModsDb(mods_db_fn)
        self.n_uniq_stats = None
        self.mods_db.check_data_covering_index_exists()
        assert agg_info.method in mh.MOD_AGG_METHOD_NAMES
        self.agg_method = agg_info.method
        self.binary_thresh = agg_info.binary_threshold
        self.expit_k = agg_info.expit_k
        self.expit_x0 = agg_info.expit_x0
        self.write_mod_lp = write_mod_lp
        self._mod_long_names = self.mods_db.get_mod_long_names()

    def get_mod_long_names(self):
        if self._mod_long_names is None:
            self._mod_long_names = self.mods_db.get_mod_long_names()
        return self._mod_long_names

    def num_uniq(self):
        if self.n_uniq_stats is None:
            self.n_uniq_stats = self.mods_db.get_num_uniq_stats()
        return self.n_uniq_stats

    def iter_uniq(self):
        # fill queue with all stats from each position
        for q_val in self.mods_db.iter_pos_scores():
            yield q_val

    def est_binary_thresh(self, pos_scores):
        valid_cov = 0
        mod_types = set(mt for read_mods in pos_scores.values()
                        for mt in read_mods.keys())
        mods_cov = dict((mt, 0) for mt in mod_types)
        for read_pos_lps in pos_scores.values():
            r_mod_types, mt_lps = zip(*read_pos_lps.items())
            mt_lps = np.array(mt_lps)
            with np.errstate(divide='ignore'):
                can_lp = np.log1p(-np.exp(mt_lps).sum())
            if can_lp > mt_lps.max():
                if np.exp(can_lp) > self.binary_thresh:
                    valid_cov += 1
            else:
                if np.exp(mt_lps.max()) > self.binary_thresh:
                    valid_cov += 1
                    mods_cov[r_mod_types[np.argmax(mt_lps)]] += 1

        if valid_cov == 0:
            return mods_cov, valid_cov
        mods_props = OrderedDict(sorted(
            (mod_type, mod_cov / valid_cov)
            for mod_type, mod_cov in mods_cov.items()))
        return mods_props, valid_cov

    def est_expit(self, pos_scores):
        def expit_scaled_estimate(lps):
            if len(lps) == 0:
                return 0
            return np.sum(1 / (1 + np.exp(-self.expit_k * (
                np.exp(lps) - self.expit_x0))))

        mod_types = set(mt for read_mods in pos_scores.values()
                        for mt in read_mods.keys())
        can_lps = []
        mods_lps = dict((mt, []) for mt in mod_types)
        for read_pos_lps in pos_scores.values():
            r_mod_types, mt_lps = zip(*read_pos_lps.items())
            mt_lps = np.array(mt_lps)
            with np.errstate(divide='ignore'):
                can_lp = np.log1p(-np.exp(mt_lps).sum())
            if can_lp > mt_lps.max():
                can_lps.append(can_lp)
            else:
                mods_lps[r_mod_types[np.argmax(mt_lps)]].append(mt_lps.max())
        mods_lsum = dict(
            (mt, expit_scaled_estimate(mt_lps))
            for mt, mt_lps in mods_lps.items())
        can_lsum = expit_scaled_estimate(can_lps)
        tot_lsum = can_lsum + sum(mods_lsum.values())
        mods_props = OrderedDict(sorted(
            (mt, mt_lsum / tot_lsum) for mt, mt_lsum in mods_lsum.items()))
        return mods_props, len(pos_scores)

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

    def compute_mod_stats(
            self, pos_data, agg_method=None, valid_read_dbids=None):
        if agg_method is None:
            agg_method = self.agg_method
        if agg_method not in mh.MOD_AGG_METHOD_NAMES:
            raise NotImplementedError(
                'No modified base proportion estimation method: {}'.format(
                    agg_method))

        pos_dbid, pos_mod_data = pos_data
        chrm, strand, pos = self.mods_db.get_pos(pos_dbid)
        mod_type_stats = defaultdict(dict)
        for read_dbid, mod_dbid, lp in pos_mod_data:
            if valid_read_dbids is not None and \
               read_dbid not in valid_read_dbids:
                continue
            mod_type_stats[read_dbid][self.mods_db.get_mod_base(mod_dbid)] = lp
        total_cov = len(mod_type_stats)
        if total_cov == 0:
            raise mh.MegaError('No valid reads cover modified base location')
        if agg_method == mh.MOD_BIN_THRESH_NAME:
            mod_props, valid_cov = self.est_binary_thresh(mod_type_stats)
        elif agg_method == mh.MOD_EXPIT:
            mod_props, valid_cov = self.est_expit(mod_type_stats)
        elif agg_method == mh.MOD_EM_NAME:
            mod_props, valid_cov = self.est_em_prop(mod_type_stats)

        strand = mh.int_strand_to_str(strand)
        mod_bases = list(mod_props.keys())
        can_base = self.mods_db.mod_to_can[mod_bases[0]]
        mod_site = ModSite(
            chrom=chrm, pos=pos, strand=strand, ref_seq=can_base,
            ref_mod_pos=0, mod_bases=mod_bases, mod_props=mod_props)
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
