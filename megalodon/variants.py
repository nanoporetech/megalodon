import os
import sys
import queue
import sqlite3
import datetime
from time import sleep
from array import array
from operator import itemgetter
from itertools import product, combinations
from collections import defaultdict, namedtuple, OrderedDict

import pysam
import numpy as np
from scipy import stats

from megalodon import (calibration, decode, logging, mapping,
                       megalodon_helper as mh)
from megalodon._version import MEGALODON_VERSION


_DEBUG_PER_READ = False

VARIANT_DATA = namedtuple('VARIANT_DATA', (
    'np_ref', 'np_alts', 'id', 'chrom', 'start', 'stop',
    'ref', 'alts', 'ref_start'))
# set default value of None for ref, alts and ref_start
VARIANT_DATA.__new__.__defaults__ = (None, None, None)

DIPLOID_MODE = 'diploid'
HAPLIOD_MODE = 'haploid'

SAMPLE_NAME = 'SAMPLE'
# specified by sam format spec
WHATSHAP_MAX_QUAL = 40
WHATSHAP_RG_ID = '1'
FIXED_VCF_MI = [
    'phasing=none',
    'INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
    'FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    'FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
    'FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
    ('FORMAT=<ID=GL,Number=G,Type=Float,' +
     'Description="Log10 likelihoods for genotypes">'),
    ('FORMAT=<ID=PL,Number=G,Type=Integer,' +
     'Description="Normalized, Phred-scaled likelihoods for genotypes">')
]
FORMAT_LOG_PROB_MI = (
    'FORMAT=<ID=LOG_PROBS,Number=A,Type=String,' +
    'Description="Per-read log10 likelihoods for alternative ' +
    'alleles (semi-colon separated)">')


#######################
##### Variants DB #####
#######################

class VarsDb(object):
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
        ('loc', OrderedDict((
            ('loc_id', 'INTEGER PRIMARY KEY'),
            ('loc_chrm', 'INTEGER'),
            ('test_start', 'INTEGER'),
            ('test_end', 'INTEGER'),
            ('var_name', 'TEXT'),
            ('pos', 'INTEGER'),
            ('ref_seq', 'TEXT')))),
        ('alt', OrderedDict((
            ('alt_id', 'INTEGER PRIMARY KEY'),
            ('alt_seq', 'TEXT')))),
        ('read', OrderedDict((
            ('read_id', 'INTEGER PRIMARY KEY'),
            ('uuid', 'TEXT'),
            ('strand', 'INTEGER')))),
        ('data', OrderedDict((
            ('score', 'FLOAT'),
            ('score_loc', 'INTEGER'),
            ('score_alt', 'INTEGER'),
            ('score_read', 'INTEGER'))))))

    # namedtuple for returning var info from a single position
    var_data = namedtuple('var_data', [
        'score', 'pos', 'ref_seq', 'var_name', 'read_id', 'chrm', 'alt_seq'])
    text_field_names = (
        'read_id', 'chrm', 'strand', 'pos', 'ref_log_prob', 'alt_log_prob',
        'ref_seq', 'alt_seq', 'var_id')

    def __init__(self, fn, read_only=True, db_safety=1,
                 loc_index_in_memory=False, chrm_index_in_memory=True,
                 alt_index_in_memory=True, uuid_index_in_memory=True,
                 uuid_strand_index_in_memory=False):
        """ Interface to database containing sequence variant statistics.

        Default settings are for optimal read_only performance.
        """
        self.fn = mh.resolve_path(fn)
        self.read_only = read_only
        self.loc_idx_in_mem = loc_index_in_memory
        self.chrm_idx_in_mem = chrm_index_in_memory
        self.alt_idx_in_mem = alt_index_in_memory
        self.uuid_idx_in_mem = uuid_index_in_memory
        self.uuid_strand_idx_in_mem = uuid_strand_index_in_memory

        if read_only:
            if not os.path.exists(fn):
                logger = logging.get_logger('vars')
                logger.error((
                    'Variant per-read database file ({}) does ' +
                    'not exist.').format(fn))
                raise mh.MegaError('Invalid variant DB filename.')
            self.db = sqlite3.connect('file:' + fn + '?mode=ro', uri=True)
        else:
            self.db = sqlite3.connect(fn, timeout=mh.SQLITE_TIMEOUT)

        self.cur = self.db.cursor()
        if self.read_only:
            # use memory mapped file access
            self.cur.execute('PRAGMA mmap_size = {}'.format(
                mh.MEMORY_MAP_LIMIT))
            if self.chrm_idx_in_mem:
                self.load_chrm_read_index()
            if self.loc_idx_in_mem:
                self.load_loc_read_index()
            if self.alt_idx_in_mem:
                self.load_alt_read_index()
            if self.uuid_idx_in_mem:
                self.load_uuid_read_index()
            if self.uuid_strand_idx_in_mem:
                self.load_uuid_strand_read_index()
        else:
            if db_safety < 2:
                # set asynchronous mode to off for max speed
                self.cur.execute('PRAGMA synchronous = OFF')
            if db_safety < 1:
                # set no rollback mode
                self.cur.execute('PRAGMA journal_mode = OFF')

            # create tables
            for tbl_name, tbl in self.db_tables.items():
                try:
                    self.cur.execute("CREATE TABLE {} ({})".format(
                        tbl_name, ','.join((
                            '{} {}'.format(*ft) for ft in tbl.items()))))
                except sqlite3.OperationalError:
                    raise mh.MegaError(
                        'Sequence variants database already exists. Either ' +
                        'provide location for new database or open in ' +
                        'read_only mode.')

            if self.loc_idx_in_mem:
                self.loc_idx = {}
            else:
                self.create_loc_index()
            if self.chrm_idx_in_mem:
                self.chrm_idx = {}
            else:
                self.create_chrm_index()
            if self.alt_idx_in_mem:
                self.alt_idx = {}
            else:
                self.create_alt_index()
            if self.uuid_idx_in_mem:
                self.uuid_idx = {}

        return

    # insert data functions
    def get_chrm_id_or_insert(self, chrm, chrm_len):
        try:
            if self.chrm_idx_in_mem:
                chrm_id = self.chrm_idx[chrm]
            else:
                chrm_id = self.cur.execute(
                    'SELECT chrm_id FROM chrm WHERE chrm=?',
                    (chrm,)).fetchone()[0]
        except (TypeError, KeyError):
            self.cur.execute('INSERT INTO chrm (chrm, chrm_len) VALUES (?,?)',
                             (chrm, chrm_len))
            chrm_id = self.cur.lastrowid
            if self.chrm_idx_in_mem:
                self.chrm_idx[chrm] = chrm_id
        return chrm_id

    def insert_chrms(self, chrm_names_and_lens):
        next_chrm_id = self.get_num_uniq_chrms() + 1
        self.cur.executemany('INSERT INTO chrm (chrm, chrm_len) VALUES (?,?)',
                             zip(*chrm_names_and_lens))
        if self.chrm_idx_in_mem:
            self.chrm_idx.update(zip(
                chrm_names_and_lens[0],
                range(next_chrm_id,
                      next_chrm_id + len(chrm_names_and_lens[0]))))
        return

    def get_loc_id_or_insert(
            self, chrm_id, test_start, test_end, pos, ref_seq, var_name):
        try:
            if self.loc_idx_in_mem:
                loc_id = self.loc_idx[(chrm_id, test_start, test_end)]
            else:
                loc_id = self.cur.execute(
                    'SELECT loc_id FROM loc WHERE loc_chrm=? AND ' +
                    'test_start=? AND test_end=?',
                    (chrm_id, test_start, test_end)).fetchone()[0]
        except (TypeError, KeyError):
            self.cur.execute(
                'INSERT INTO loc (loc_chrm, test_start, test_end, ' +
                'pos, ref_seq, var_name) VALUES (?,?,?,?,?,?)',
                (chrm_id, test_start, test_end, pos, ref_seq, var_name))
            loc_id = self.cur.lastrowid
            if self.loc_idx_in_mem:
                self.loc_idx[(chrm_id, test_start, test_end)] = loc_id
        return loc_id

    def get_loc_ids_or_insert(self, r_var_scores, chrm_id):
        """ Extract all location IDs and add those locations not currently
        found in the database
        """
        if len(r_var_scores) == 0: return []

        r_locs = dict((
            ((chrm_id, test_start, test_end), (pos, ref_seq, var_name))
            for pos, _, ref_seq, _, var_name,
            test_start, test_end in r_var_scores))
        if self.loc_idx_in_mem:
            locs_to_add = tuple(set(r_locs).difference(self.loc_idx))
        else:
            test_starts, test_ends = map(
                set, list(zip(*r_locs.keys()))[1:])
            loc_ids = dict((
                ((chrm_id, test_start, test_end), loc_id)
                for chrm_id, test_start, test_end, loc_id in  self.cur.execute(
                        ('SELECT loc_chrm, test_start, test_end, loc_id ' +
                         'FROM loc WHERE loc_chrm=? AND ' +
                         'test_start in ({0}) AND test_end in ({1})').format(
                             ','.join(['?',] * len(test_starts)),
                             ','.join(['?',] * len(test_ends))),
                        (chrm_id, *test_starts, *test_ends)).fetchall()))
            locs_to_add = tuple(set(r_locs).difference(loc_ids))

        if len(locs_to_add) > 0:
            next_loc_id = self.get_num_uniq_var_loc() + 1
            self.cur.executemany(
                'INSERT INTO loc (loc_chrm, test_start, test_end, ' +
                'pos, ref_seq, var_name) VALUES (?,?,?,?,?,?)',
                ((*loc_key, *r_locs[loc_key]) for loc_key in locs_to_add))

        loc_idx = self.loc_idx if self.loc_idx_in_mem else loc_ids
        if len(locs_to_add) > 0:
            loc_idx.update(zip(
                locs_to_add,
                range(next_loc_id, next_loc_id + len(locs_to_add))))
        r_loc_ids = [
            loc_idx[(chrm_id, test_start, test_end)]
            for _, _, _, _, _, test_start, test_end in r_var_scores]

        return r_loc_ids

    def get_alt_id_or_insert(self, alt_seq):
        try:
            if self.alt_idx_in_mem:
                alt_id = self.alt_idx[alt_seq]
            else:
                alt_id = self.cur.execute(
                    'SELECT alt_id FROM alt WHERE alt_seq=?',
                    (alt_seq,)).fetchone()[0]
        except (TypeError, KeyError):
            self.cur.execute(
                'INSERT INTO alt (alt_seq) VALUES (?)',
                (alt_seq,))
            alt_id = self.cur.lastrowid
            if self.alt_idx_in_mem:
                self.alt_idx[alt_seq] = alt_id
        return alt_id

    def get_alt_ids_or_insert(self, r_var_scores):
        if len(r_var_scores) == 0: return []

        logger = logging.get_logger()
        r_seqs_and_lps = [
            tuple(zip(alt_seqs, alt_lps))
            for _, alt_lps, _, alt_seqs, _, _, _ in r_var_scores]
        r_uniq_seqs = set(seq_lp_i[0] for loc_seqs_lps in r_seqs_and_lps
                          for seq_lp_i in loc_seqs_lps)
        if self.alt_idx_in_mem:
            alts_to_add = tuple(r_uniq_seqs.difference(self.alt_idx))
        else:
            alt_ids = dict((
                (alt_seq, alt_id)
                for alt_seq, alt_id in  self.cur.execute(
                        ('SELECT alt_seq, alt_id ' +
                         'FROM alt WHERE alt_seq in ({})').format(
                             ','.join(['?',] * len(r_uniq_seqs))),
                        r_uniq_seqs).fetchall()))
            alts_to_add = tuple(r_uniq_seqs.difference(alt_ids))

        if len(alts_to_add) > 0:
            next_alt_id = self.get_num_uniq_alt_seqs() + 1
            self.cur.executemany(
                'INSERT INTO alt (alt_seq) VALUES (?)',
                ((alt_seq,) for alt_seq in alts_to_add))

        alt_idx = self.alt_idx if self.alt_idx_in_mem else alt_ids
        if len(alts_to_add) > 0:
            alt_idx.update(zip(
                alts_to_add,
                range(next_alt_id, next_alt_id + len(alts_to_add))))
        r_alt_ids = [
            tuple((alt_idx[alt_seq], alt_lp)
                  for alt_seq, alt_lp in loc_seqs_lps)
            for loc_seqs_lps in r_seqs_and_lps]

        return r_alt_ids

    def get_read_id_or_insert(self, uuid):
        try:
            if self.uuid_idx_in_mem:
                read_id = self.uuid_idx[uuid]
            else:
                read_id = self.cur.execute(
                    'SELECT read_id FROM read WHERE uuid=?',
                    (uuid,)).fetchone()[0]
        except (TypeError, KeyError):
            self.cur.execute('INSERT INTO read (uuid) VALUES (?)', (uuid,))
            read_id = self.cur.lastrowid
            if self.uuid_idx_in_mem:
                self.uuid_idx[uuid] = read_id
        return read_id

    def insert_read_scores(self, r_var_scores, uuid, chrm, strand):
        self.cur.execute('INSERT INTO read (uuid, strand) VALUES (?,?)',
                         (uuid, strand))
        if len(r_var_scores) == 0: return

        read_id = self.cur.lastrowid
        chrm_id = self.get_chrm_id(chrm)
        loc_ids = self.get_loc_ids_or_insert(r_var_scores, chrm_id)
        alt_ids = self.get_alt_ids_or_insert(r_var_scores)

        read_insert_data = ((alt_lp, loc_id, alt_id, read_id)
                            for loc_id, loc_alts in zip(loc_ids, alt_ids)
                            for alt_id, alt_lp in loc_alts)

        self.cur.executemany(
            'INSERT INTO data VALUES (?,?,?,?)', read_insert_data)
        return

    def insert_data(self, score, loc_id, alt_id, read_id):
        self.cur.execute(
            'INSERT INTO data (score, score_loc, score_alt, score_read) ' +
            'VALUES (?,?,?,?)', (score, loc_id, alt_id, read_id))
        return self.cur.lastrowid

    # create and load index functions
    def create_chrm_index(self):
        self.cur.execute('CREATE UNIQUE INDEX chrm_idx ON chrm(chrm)')
        return

    def load_chrm_read_index(self):
        self.chrm_read_idx = dict(self.cur.execute(
            'SELECT chrm_id, chrm FROM chrm').fetchall())
        return

    def load_uuid_read_index(self):
        self.uuid_read_idx = dict(self.cur.execute(
            'SELECT read_id, uuid FROM read').fetchall())
        return

    def load_uuid_strand_read_index(self):
        self.uuid_strand_read_idx = dict(
            (read_id, (uuid, strand)) for read_id, uuid, strand in
            self.cur.execute(
                'SELECT read_id, uuid, strand FROM read').fetchall())
        return

    def create_alt_index(self):
        self.cur.execute('CREATE UNIQUE INDEX alt_idx ON alt(alt_seq)')
        return

    def load_alt_read_index(self):
        self.alt_read_idx = dict(self.cur.execute(
            'SELECT alt_id, alt_seq FROM alt').fetchall())
        return

    def create_loc_index(self):
        self.cur.execute('CREATE UNIQUE INDEX loc_idx ON loc' +
                         '(loc_chrm, test_start, test_end)')
        return

    def load_loc_read_index(self):
        self.loc_read_idx = dict(
            (loc_id, (pos, ref_seq, var_name))
            for loc_id, pos, ref_seq, var_name in self.cur.execute(
                    'SELECT loc_id, pos, ref_seq, var_name ' +
                    'FROM loc').fetchall())
        return

    def create_data_covering_index(self):
        self.cur.execute('CREATE INDEX data_cov_idx ON data(' +
                         'score_loc, score_alt, score_read, score)')
        return

    # reader functions
    def get_chrm_id(self, chrm):
        try:
            if self.chrm_idx_in_mem:
                chrm_id = self.chrm_idx[chrm]
            else:
                chrm_id = self.cur.execute(
                    'SELECT chrm_id FROM chrm WHERE chrm=?',
                    (chrm,)).fetchone()[0]
        except (TypeError, KeyError):
            raise mh.MegaError('Reference record (chromosome) not found in ' +
                               'database.')
        return chrm_id

    def get_chrm(self, chrm_id):
        try:
            if self.chrm_idx_in_mem:
                chrm = self.chrm_read_idx[chrm_id]
            else:
                chrm = self.cur.execute(
                    'SELECT chrm FROM chrm WHERE chrm_id=?',
                    (chrm_id,)).fetchone()[0]
        except (TypeError, KeyError):
            raise mh.MegaError('Reference record (chromosome) not found in ' +
                               'vars database.')
        return chrm

    def get_all_chrm_and_lens(self):
        try:
            return tuple(map(tuple, zip(*self.cur.execute(
                'SELECT chrm, chrm_len FROM chrm').fetchall())))
        except sqlite3.OperationalError:
            raise mh.MegaError(
                'Old megalodon database scheme detected. Please re-run ' +
                'megalodon processing or downgrade megalodon installation.')

    def get_alt_seq(self, alt_id):
        try:
            if self.alt_idx_in_mem:
                alt_seq = self.alt_read_idx[alt_id]
            else:
                alt_seq = self.cur.execute(
                    'SELECT alt_seq FROM alt WHERE alt_id=?',
                    (alt_id,)).fetchone()[0]
        except (TypeError, KeyError):
            raise mh.MegaError('Alt sequence not found in vars database.')
        return alt_seq

    def get_uuid(self, read_id):
        try:
            if self.uuid_idx_in_mem:
                uuid = self.uuid_read_idx[read_id]
            else:
                uuid = self.cur.execute(
                    'SELECT uuid FROM read WHERE read_id=?',
                    (read_id,)).fetchone()[0]
        except (TypeError, KeyError):
            raise mh.MegaError('Read ID not found in vars database.')
        return uuid

    def get_uuid_strand(self, read_id):
        try:
            if self.uuid_strand_idx_in_mem:
                uuid_strand = self.uuid_strand_read_idx[read_id]
            else:
                uuid_strand = self.cur.execute(
                    'SELECT uuid, strand FROM read WHERE read_id=?',
                    (read_id,)).fetchone()
        except (TypeError, KeyError):
            raise mh.MegaError('Read ID not found in vars database.')
        return uuid_strand

    def get_num_uniq_chrms(self):
        num_chrms = self.cur.execute(
            'SELECT MAX(chrm_id) FROM chrm').fetchone()[0]
        if num_chrms is None:
            num_chrms = 0
        return num_chrms

    def get_num_uniq_var_loc(self):
        num_locs = self.cur.execute('SELECT MAX(loc_id) FROM loc').fetchone()[0]
        if num_locs is None:
            num_locs = 0
        return num_locs

    def get_num_uniq_alt_seqs(self):
        num_alts = self.cur.execute('SELECT MAX(alt_id) FROM alt').fetchone()[0]
        if num_alts is None:
            num_alts = 0
        return num_alts

    def iter_locs(self):
        for loc in self.cur.execute(
                'SELECT loc_id, loc_chrm, pos, ref_seq, var_name ' +
                'FROM loc').fetchall():
            yield loc
        return

    def iter_data(self):
        for data in self.cur.execute(
                'SELECT score, uuid, strand, alt_seq, ref_seq, pos, ' +
                'var_name, test_end, test_start, chrm, chrm_len ' +
                'FROM data ' +
                'INNER JOIN read ON data.score_read = read.read_id ' +
                'INNER JOIN alt ON data.score_alt = alt.alt_id ' +
                'INNER JOIN loc ON data.score_loc = loc.loc_id ' +
                'INNER JOIN chrm ON loc.loc_chrm = chrm.chrm_id').fetchall():
            yield data
        return

    def get_loc_stats(self, loc_data, return_uuids=False):
        read_id_conv = self.get_uuid if return_uuids else lambda x: x
        # these attributes are specified in self.iter_locs
        loc_id, chrm_id, pos, ref_seq, var_name = loc_data
        chrm = self.get_chrm(chrm_id)
        return [
            self.var_data(
                score, pos, ref_seq, var_name, read_id_conv(read_id),
                chrm, self.get_alt_seq(alt_id))
            for score, read_id, alt_id, loc_id in self.cur.execute(
                    'SELECT score, score_read, score_alt, score_loc FROM data ' +
                    'WHERE score_loc=?', (loc_id,)).fetchall()]

    def close(self):
        self.db.commit()
        self.db.close()
        return


############################
##### Helper Functions #####
############################

def logsumexp(x):
    x_max = x.max()
    return np.log(np.sum(np.exp(x - x_max))) + x_max


####################################
##### Per-read Variant Scoring #####
####################################

def write_per_read_debug(
        var_ref_pos, var_id, read_ref_pos, np_s_var_ref_seq, np_s_var_alt_seqs,
        np_s_context_seqs, loc_contexts_ref_lps, loc_contexts_alts_lps,
        w_context, logger):
    ref_seq = mh.int_to_seq(np_s_var_ref_seq)
    if read_ref_pos.strand == -1:
        ref_seq = mh.revcomp(ref_seq)
    alts_seq = [mh.int_to_seq(np_alt) for np_alt in np_s_var_alt_seqs]
    if read_ref_pos.strand == -1:
        alts_seq = [mh.revcomp(alt_seq) for alt_seq in alts_seq]
    ','.join(alts_seq)

    context_seqs = []
    for up_context_seq, dn_context_seq in np_s_context_seqs:
        up_seq, dn_seq = (mh.int_to_seq(up_context_seq),
                          mh.int_to_seq(dn_context_seq))
        if read_ref_pos.strand == -1:
            up_seq, dn_seq = mh.revcomp(dn_seq), mh.revcomp(up_seq)
        context_seqs.append((up_seq, dn_seq))
    out_txt = '\n'
    for ref_lp, alt_lps, (up_seq, dn_seq) in zip(
            loc_contexts_ref_lps, zip(*loc_contexts_alts_lps), context_seqs):
        out_txt += ('VARIANT_FULL_DATA: {}\t{}\t{}\t{}\t{}[{}]{}\t{}\t' +
                    '{:.2f}\t{}\t{}\n').format(
                        read_ref_pos.chrm, read_ref_pos.strand, var_ref_pos,
                        var_id, up_seq, ref_seq, dn_seq, ','.join(alts_seq),
                        ref_lp, ','.join(('{:.2f}'.format(alt_lp)
                                          for alt_lp in alt_lps)),
                        'WITH_CONTEXT' if w_context else 'NO_CONTEXT')
    logger.debug(out_txt)
    return

def score_seq(tpost, seq, tpost_start=0, tpost_end=None,
              all_paths=False):
    """Score a section of log transition posteriors against a proposed sequence
    using a global mapping.
    :param tpost: `ndarray` containing log transition posteriors to be scored
    :param seq: `ndarray` containing integers encoding proposed sequence
    :param tpost_start: start position within post (Default: 0)
    :param tpost_end: end position within post (Default: full posterior)
    :param all_paths: boolean to produce the forwards all paths score
        (default Viterbi best path)
    """
    seq = seq.astype(np.uintp)
    if tpost_end is None:
        tpost_end = post.shape[0]
    if seq.shape[0] >= tpost_end - tpost_start:
        raise mh.MegaError('Mapped signal too short for proposed sequence.')

    score = decode.score_seq(tpost, seq, tpost_start, tpost_end, all_paths)
    if np.isnan(score):
        raise mh.MegaError('Score computation error (likely  memory error).')

    return score

def score_variants_with_context(
        vars_iter, ref_to_block, r_post, read_ref_fwd_seq, r_var_calls,
        read_cached_scores, read_ref_pos=None, all_paths=False,
        calibrate_llr=lambda *x: x[0], logger=None):
    for (np_s_var_ref_seq, np_s_var_alt_seqs, np_s_context_seqs,
         s_ref_start, s_ref_end, variant) in vars_iter:
        ref_cntxt_ref_lp, ref_cntxt_alt_lps = read_cached_scores[(
            variant.id, variant.start, variant.stop)]

        blk_start  = ref_to_block[s_ref_start]
        blk_end = ref_to_block[s_ref_end]

        # filter to context sequences that can be evaluated given the
        # assigned context signal/blocks
        max_var_len = max(np_s_var_ref_seq.shape[0], max(
            var_alt_seq.shape[0] for var_alt_seq in np_s_var_alt_seqs))
        filt_np_s_context_seqs = [
            (up_seq, dn_seq) for up_seq, dn_seq in np_s_context_seqs
            if blk_end - blk_start > len(up_seq) + len(dn_seq) + max_var_len]

        # skip first (reference) context seq as this was cached
        ref_context_seqs = (
            np.concatenate([up_context_seq, np_s_var_ref_seq, dn_context_seq])
            for up_context_seq, dn_context_seq in filt_np_s_context_seqs[1:])
        loc_contexts_ref_lps = np.array([ref_cntxt_ref_lp] + [
            score_seq(r_post, ref_seq, blk_start, blk_end, all_paths)
            for ref_seq in ref_context_seqs])
        loc_ref_lp = logsumexp(loc_contexts_ref_lps)

        loc_alt_llrs = []
        if _DEBUG_PER_READ:
            loc_contexts_alts_lps = []
        for np_s_var_alt_seq, var_alt_seq, ref_cntxt_alt_lp in zip(
                np_s_var_alt_seqs, variant.alts, ref_cntxt_alt_lps):
            alt_context_seqs = (
                np.concatenate([
                    up_context_seq, np_s_var_alt_seq, dn_context_seq])
                for up_context_seq, dn_context_seq in filt_np_s_context_seqs[1:])
            loc_contexts_alt_lps = np.array([ref_cntxt_alt_lp,] + [
                score_seq(r_post, alt_seq, blk_start, blk_end, all_paths)
                for alt_seq in alt_context_seqs])
            loc_alt_lp = logsumexp(loc_contexts_alt_lps)

            if _DEBUG_PER_READ:
                loc_contexts_alts_lps.append(loc_contexts_alt_lps)
            # calibrate log probs
            loc_alt_llrs.append(calibrate_llr(
                loc_ref_lp - loc_alt_lp, variant.ref, var_alt_seq))

        # due to calibration mutli-allelic log likelihoods could result in
        # inferred negative reference likelihood, so re-normalize here
        loc_alt_log_ps = calibration.compute_log_probs(np.array(loc_alt_llrs))

        if _DEBUG_PER_READ and logger is not None and read_ref_pos is not None:
            write_per_read_debug(
                variant.start, variant.id, read_ref_pos,
                np_s_var_ref_seq, np_s_var_alt_seqs, filt_np_s_context_seqs,
                loc_contexts_ref_lps, loc_contexts_alts_lps, True, logger)

        r_var_calls.append((
            variant.ref_start, loc_alt_log_ps, variant.ref,
            variant.alts, variant.id, variant.start,
            variant.start + variant.np_ref.shape[0]))

    return r_var_calls

def score_variants_independently(
        vars_iter, ref_to_block, r_post, read_ref_fwd_seq, read_ref_pos=None,
        all_paths=False, calibrate_llr=lambda *x: x[0],
        context_min_alt_prob=1.0, logger=None):
    r_var_calls = []
    read_cached_scores = {}
    filt_read_variants = []
    for (np_s_var_ref_seq, np_s_var_alt_seqs, np_s_context_seqs,
         s_ref_start, s_ref_end, variant) in vars_iter:
        blk_start  = ref_to_block[s_ref_start]
        blk_end = ref_to_block[s_ref_end]
        if blk_end - blk_start <= max(
                len(up_seq) + len(dn_seq)
                for up_seq, dn_seq in np_s_context_seqs) + max(
                        np_s_var_ref_seq.shape[0], max(
                            var_alt_seq.shape[0]
                            for var_alt_seq in np_s_var_alt_seqs)):
            # no valid mapping over large inserted query bases
            # i.e. need as many "events/strides" as bases for valid mapping
            continue

        np_ref_seq = np.concatenate([
            np_s_context_seqs[0][0], np_s_var_ref_seq, np_s_context_seqs[0][1]])
        loc_ref_lp = score_seq(
            r_post, np_ref_seq, blk_start, blk_end, all_paths)

        loc_alt_lps = []
        loc_alt_llrs = []
        if _DEBUG_PER_READ:
            loc_contexts_alts_lps = []
        for np_s_var_alt_seq, var_alt_seq in zip(
                np_s_var_alt_seqs, variant.alts):
            np_alt_seq = np.concatenate([
                np_s_context_seqs[0][0], np_s_var_alt_seq,
                np_s_context_seqs[0][1]])
            loc_alt_lp = score_seq(
                r_post, np_alt_seq, blk_start, blk_end, all_paths)
            loc_alt_lps.append(loc_alt_lp)
            if _DEBUG_PER_READ:
                loc_contexts_alts_lps.append(np.array([loc_alt_lp,]))
            # calibrate log probs
            loc_alt_llrs.append(calibrate_llr(
                loc_ref_lp - loc_alt_lp, variant.ref, var_alt_seq))

        # due to calibration mutli-allelic log likelihoods could result in
        # inferred negative reference likelihood, so re-normalize here
        loc_alt_log_ps = calibration.compute_log_probs(np.array(loc_alt_llrs))

        if _DEBUG_PER_READ and logger is not None and read_ref_pos is not None:
            write_per_read_debug(
                variant.start, variant.id, read_ref_pos,
                np_s_var_ref_seq, np_s_var_alt_seqs, np_s_context_seqs,
                np.array([loc_ref_lp,]), loc_contexts_alts_lps, False, logger)

        if sum(np.exp(loc_alt_log_ps)) >= context_min_alt_prob:
            # if the probability of a variant scored independently is greater
            # than the minimal threshold then save it for future testing
            filt_read_variants.append(variant)
            read_cached_scores[(variant.id, variant.start, variant.stop)] = (
                loc_ref_lp, loc_alt_lps)
        else:
            # if a variant is less probable independently then save this as
            # the final variant score
            r_var_calls.append((
                variant.ref_start, loc_alt_log_ps, variant.ref,
                variant.alts, variant.id, variant.start,
                variant.start + variant.np_ref.shape[0]))

    return r_var_calls, read_cached_scores, filt_read_variants

def call_read_vars(
        vars_data, read_ref_pos, strand_read_np_ref_seq, ref_to_block, r_post):
    if read_ref_pos.end - read_ref_pos.start <= 2 * vars_data.edge_buffer:
        raise mh.MegaError('Mapped region too short for variant calling.')

    # convert to forward strand sequence in order to annotate with variants
    read_ref_fwd_seq = (strand_read_np_ref_seq if read_ref_pos.strand == 1 else
                        mh.revcomp_np(strand_read_np_ref_seq))
    # call all variantss overlapping this read
    logger = logging.get_logger('per_read_vars')
    read_variants = vars_data.fetch_read_variants(
        read_ref_pos, read_ref_fwd_seq)

    # first pass over variants assuming the reference ground truth
    # (not including context variants)
    vars_iter = vars_data.iter_vars(
        read_variants, read_ref_pos, read_ref_fwd_seq, context_max_dist=0)
    (r_var_calls, read_cached_scores,
     filt_read_variants) = score_variants_independently(
         vars_iter, ref_to_block, r_post, read_ref_fwd_seq, read_ref_pos,
         vars_data.all_paths, vars_data.calib_table.calibrate_llr,
         vars_data.context_min_alt_prob, logger)

    # second round for variants with some evidence for alternative alleles
    # process with other potential variants as context
    vars_iter = vars_data.iter_vars(
        filt_read_variants, read_ref_pos, read_ref_fwd_seq)
    r_var_calls = score_variants_with_context(
        vars_iter, ref_to_block, r_post, read_ref_fwd_seq, r_var_calls,
        read_cached_scores, read_ref_pos, vars_data.all_paths,
        vars_data.calib_table.calibrate_llr, logger)

    # re-sort variants after adding context-included computations
    return sorted(r_var_calls, key=lambda x: x[0])


#######################################
##### Propose All Variant Helpers #####
#######################################

def score_all_single_deletions(
        read_np_seq, r_post, ref_to_block, start, end, context_bases):
    vars_iter = []
    for del_pos in range(start, end):
        ctxt_start, ctxt_end = (
            max(0, del_pos - context_bases),
            min(read_np_seq.shape[0], del_pos + 1 + context_bases))
        np_ref = read_np_seq[del_pos:del_pos + 1]
        np_alts = [np.array([], dtype=np.uintp),]
        ctxt_seqs = [[read_np_seq[ctxt_start:del_pos],
                      read_np_seq[del_pos + 1:ctxt_end]],]
        # technically the variant ref_start should be 1-based, but it will be
        # returned so just use 0-based here.
        variant = VARIANT_DATA(
            np_ref=np_ref, np_alts=np_alts, id=None, chrom=None, start=del_pos,
            stop=del_pos + 1, ref=np_ref, alts=np_alts, ref_start=del_pos)
        vars_iter.append((
            np_ref, np_alts, ctxt_seqs, ctxt_start, ctxt_end, variant))

    return [(x[1][0], x[0]) for x in score_variants_independently(
        vars_iter, ref_to_block, r_post, read_np_seq)[0]]

def score_all_single_insertions(
        read_np_seq, r_post, ref_to_block, start, end, context_bases):
    vars_iter = []
    for del_pos in range(start, end):
        ctxt_start, ctxt_end = (
            max(0, del_pos - context_bases),
            min(read_np_seq.shape[0], del_pos + 1 + context_bases))
        np_ref = np.array([], dtype=np.uintp)
        np_alts = [np.array([b,], dtype=np.uintp) for b in range(4)]
        ctxt_seqs = [[read_np_seq[ctxt_start:del_pos],
                      read_np_seq[del_pos + 1:ctxt_end]],]
        # technically the variant ref_start should be 1-based, but it will be
        # returned so just use 0-based here.
        variant = VARIANT_DATA(
            np_ref=np_ref, np_alts=np_alts, id=None, chrom=None, start=del_pos,
            stop=del_pos + 1, ref=np_ref, alts=np_alts, ref_start=del_pos)
        vars_iter.append((
            np_ref, np_alts, ctxt_seqs, ctxt_start, ctxt_end, variant))

    return [(score, pos, alt) for pos, scores, _, alts, _, _, _ in
            score_variants_independently(
                vars_iter, ref_to_block, r_post, read_np_seq)[0]
            for score, alt in zip(scores, alts)]


###################################
##### Per-read Variant Output #####
###################################

def log_prob_to_phred(log_prob):
    with np.errstate(divide='ignore'):
        return -10 * np.log10(1 - np.exp(log_prob))

def simplify_var_seq(ref_seq, alt_seq):
    trim_before = trim_after = 0
    while (len(ref_seq) > 0 and len(alt_seq) > 0 and
           ref_seq[0] == alt_seq[0]):
        trim_before += 1
        ref_seq = ref_seq[1:]
        alt_seq = alt_seq[1:]
    while (len(ref_seq) > 0 and len(alt_seq) > 0 and
           ref_seq[-1] == alt_seq[-1]):
        trim_after += 1
        ref_seq = ref_seq[:-1]
        alt_seq = alt_seq[:-1]

    return ref_seq, alt_seq, trim_before, trim_after

def iter_non_overlapping_variants(r_var_calls):
    def get_max_prob_allele_var(var_grp):
        """ For overlapping variantss return the variant with the highest
        probability single allele as this one will be added to the reference
        sequence.

        More complex chained variantss could be handled, but are not here.
        For example, a 5 base deletion covering 2 single base swap variants
        could validly result in 2 alternative single base swap alleles, but the
        logic here would only allow one of those alternatives since they
        are covered by the same reference deletion. There are certainly many
        more edge cases than this and each one would require specific logic.
        This likely covers the majority of valid cases and limiting to
        50 base indels by default limits the scope of this issue.
        """
        most_prob_var = None
        for var_data in var_grp:
            with np.errstate(divide='ignore'):
                ref_lp = np.log1p(-np.exp(var_data[1]).sum())
            var_max_lp = max(ref_lp, max(var_data[1]))
            if most_prob_var is None or var_max_lp > most_prob_var[0]:
                most_prob_var = (var_max_lp, ref_lp, var_data)

        _, ref_lp, (var_pos, alt_lps, var_ref_seq,
                    var_alt_seqs, _, _, _) = most_prob_var
        return var_pos, alt_lps, var_ref_seq, var_alt_seqs, ref_lp


    if len(r_var_calls) == 0: return
    r_var_calls_iter = iter(r_var_calls)
    # initialize var_grp with first var
    var_data = next(r_var_calls_iter)
    prev_var_end = var_data[0] + len(var_data[2])
    var_grp = [var_data]
    for var_data in sorted(r_var_calls_iter, key=itemgetter(0)):
        if var_data[0] < prev_var_end:
            prev_var_end = max(var_data[0] + len(var_data[2]), prev_var_end)
            var_grp.append(var_data)
        else:
            yield get_max_prob_allele_var(var_grp)
            prev_var_end = var_data[0] + len(var_data[2])
            var_grp = [var_data]

    # yeild last var grp data
    yield get_max_prob_allele_var(var_grp)
    return

def annotate_variants(r_start, ref_seq, r_var_calls, strand):
    """ Annotate reference sequence with called variants.

    Note: Reference sequence is in read orientation and variant calls are in
    genome coordiates.
    """
    var_seqs, var_quals, var_cigar = [], [], []
    prev_pos, curr_match = 0, 0
    # ref_seq is read-centric so flop order to process vars in genomic order
    if strand == -1:
        ref_seq = ref_seq[::-1]
    for (var_pos, alt_lps, var_ref_seq, var_alt_seqs,
         ref_lp) in iter_non_overlapping_variants(r_var_calls):
        prev_len = var_pos - r_start - prev_pos
        # called canonical
        if ref_lp >= max(alt_lps):
            var_seqs.append(
                ref_seq[prev_pos:var_pos - r_start + len(var_ref_seq)])
            var_quals.extend(
                ([WHATSHAP_MAX_QUAL] * prev_len) +
                ([min(log_prob_to_phred(ref_lp), WHATSHAP_MAX_QUAL)] *
                 len(var_ref_seq)))
            curr_match += prev_len + len(var_ref_seq)
        else:
            alt_seq = var_alt_seqs[np.argmax(alt_lps)]
            # complement since ref_seq is complement seq
            # (not reversed; see loop init)
            read_alt_seq = alt_seq if strand == 1 else mh.comp(alt_seq)
            var_seqs.append(ref_seq[prev_pos:var_pos - r_start] + read_alt_seq)
            var_quals.extend(
                ([WHATSHAP_MAX_QUAL] * prev_len) +
                ([min(log_prob_to_phred(max(alt_lps)), WHATSHAP_MAX_QUAL)] *
                 len(alt_seq)))

            # add cigar information for variant
            t_ref_seq, t_alt_seq, t_before, t_after = simplify_var_seq(
                var_ref_seq, alt_seq)
            curr_match += t_before
            var_cigar.append((7, curr_match + prev_len))
            if len(t_alt_seq) == len(t_ref_seq):
                var_cigar.append((8, len(t_alt_seq)))
            elif len(t_alt_seq) > len(t_ref_seq):
                # left justify mismatch bases in complex insertion
                if len(t_ref_seq) != 0:
                    var_cigar.append((8, len(t_ref_seq)))
                var_cigar.append((1, len(t_alt_seq) - len(t_ref_seq)))
            else:
                # left justify mismatch bases in complex deletion
                if len(t_alt_seq) != 0:
                    var_cigar.append((8, len(t_alt_seq)))
                var_cigar.append((2, len(t_ref_seq) - len(t_alt_seq)))
            curr_match = t_after
        prev_pos = var_pos - r_start + len(var_ref_seq)

    var_seqs.append(ref_seq[prev_pos:])
    var_seq = ''.join(var_seqs)
    if strand == -1:
        var_seq = var_seq[::-1]
    len_remain = len(ref_seq) - prev_pos
    var_quals.extend([WHATSHAP_MAX_QUAL] * len_remain)
    if strand == -1:
        var_quals = var_quals[::-1]
    var_quals = list(map(int, var_quals))
    var_cigar.append((7, len_remain + curr_match))
    if strand == -1:
        var_cigar = var_cigar[::-1]

    return var_seq, var_quals, var_cigar

def _get_variants_queue(
        vars_q, vars_conn, vars_db_fn, vars_txt_fn, db_safety, pr_refs_fn,
        pr_ref_filts, whatshap_map_fn, ref_names_and_lens, ref_fn,
        loc_index_in_memory):
    def write_whatshap_alignment(
            read_id, var_seq, var_quals, chrm, strand, r_st, var_cigar):
        a = pysam.AlignedSegment()
        a.query_name = read_id
        a.flag = 0 if strand == 1 else 16
        a.reference_id = whatshap_map_fp.get_tid(chrm)
        a.reference_start = r_st
        a.template_length = len(var_seq)
        a.mapping_quality = WHATSHAP_MAX_QUAL
        a.set_tags([('RG', WHATSHAP_RG_ID)])

        # convert to reference based sequence
        if strand == -1:
            var_seq = mh.revcomp(var_seq)
            var_quals = var_quals[::-1]
            var_cigar = var_cigar[::-1]
        a.query_sequence = var_seq
        a.query_qualities = array('B', var_quals)
        a.cigartuples = var_cigar
        whatshap_map_fp.write(a)

        return

    def store_var_call(
            r_var_calls, read_id, chrm, strand, r_start, ref_seq, read_len,
            q_st, q_en, cigar, been_warned):
        try:
            vars_db.insert_read_scores(r_var_calls, read_id, chrm, strand)
        except Exception as e:
            if not been_warned:
                logger.warning(
                    'Error inserting sequence variant scores into database. ' +
                    'See log debug output for error details.')
                been_warned = True
            import traceback
            logger.debug(
                'Error inserting modified base scores into database: ' +
                str(e) + '\n' + traceback.format_exc())

        if vars_txt_fp is not None and len(r_var_calls) > 0:
            var_out_text = ''
            for (pos, alt_lps, var_ref_seq, var_alt_seqs, var_id,
                 test_start, test_end) in r_var_calls:
                with np.errstate(divide='ignore'):
                    ref_lp = np.log1p(-np.exp(alt_lps).sum())
                var_out_text += '\n'.join((
                    ('\t'.join('{}' for _ in vars_db.text_field_names)).format(
                        read_id, chrm, strand, pos, ref_lp, alt_lp,
                        var_ref_seq, var_alt_seq, var_id)
                    for alt_lp, var_alt_seq in zip(
                            alt_lps, var_alt_seqs))) + '\n'
            vars_txt_fp.write(var_out_text)
        if do_ann_vars:
            if not mapping.read_passes_filters(
                    pr_ref_filts, read_len, q_st, q_en, cigar):
                return
            var_seq, var_quals, var_cigar = annotate_variants(
                r_start, ref_seq, r_var_calls, strand)
            if pr_refs_fn is not None:
                pr_refs_fp.write('>{}\n{}\n'.format(read_id, var_seq))
            if whatshap_map_fn is not None:
                write_whatshap_alignment(
                    read_id, var_seq, var_quals, chrm, strand, r_start,
                    var_cigar)

        return been_warned


    logger = logging.get_logger('vars_getter')
    been_warned = False

    vars_db = VarsDb(vars_db_fn, db_safety=db_safety, read_only=False,
                     loc_index_in_memory=loc_index_in_memory)
    vars_db.insert_chrms(ref_names_and_lens)
    vars_db.create_chrm_index()
    if vars_txt_fn is None:
        vars_txt_fp = None
    else:
        vars_txt_fp = open(vars_txt_fn, 'w')
        vars_txt_fp.write('\t'.join(vars_db.text_field_names) + '\n')

    if pr_refs_fn is not None:
        pr_refs_fp = open(pr_refs_fn, 'w')

    if whatshap_map_fn is not None:
        _, map_fmt = os.path.splitext(whatshap_map_fn)
        if map_fmt == '.bam': w_mode = 'wb'
        elif map_fmt == '.cram': w_mode = 'wc'
        elif map_fmt == '.sam': w_mode = 'w'
        else:
            raise mh.MegaError('Invalid mapping output format')
        header = {
            'HD': {'VN': '1.4'},
            'SQ': [{'LN': ref_len, 'SN': ref_name}
                   for ref_name, ref_len in sorted(zip(*ref_names_and_lens))],
            'RG': [{'ID':WHATSHAP_RG_ID, 'SM':SAMPLE_NAME},]}
        whatshap_map_fp = pysam.AlignmentFile(
            whatshap_map_fn, w_mode, header=header, reference_filename=ref_fn)

    do_ann_vars = whatshap_map_fn is not None or pr_refs_fn is not None

    while True:
        try:
            r_var_calls, (read_id, chrm, strand, r_start, ref_seq, read_len,
                          q_st, q_en, cigar) = vars_q.get(block=False)
        except queue.Empty:
            if vars_conn.poll():
                break
            sleep(0.001)
            continue
        been_warned = store_var_call(
            r_var_calls, read_id, chrm, strand, r_start, ref_seq, read_len,
            q_st, q_en, cigar, been_warned)
    while not vars_q.empty():
        r_var_calls, (read_id, chrm, strand, r_start, ref_seq, read_len,
                      q_st, q_en, cigar) = vars_q.get(block=False)
        been_warned = store_var_call(
            r_var_calls, read_id, chrm, strand, r_start, ref_seq, read_len,
            q_st, q_en, cigar, been_warned)

    if vars_txt_fp is not None: vars_txt_fp.close()
    if pr_refs_fn is not None: pr_refs_fp.close()
    if whatshap_map_fn is not None: whatshap_map_fp.close()
    vars_db.create_alt_index()
    if vars_db.loc_idx_in_mem:
        vars_db.create_loc_index()
    vars_db.create_data_covering_index()
    vars_db.close()

    return


######################
##### VCF Reader #####
######################

class VarData(object):
    def check_vars_match_ref(
            self, vars_idx, contigs, aligner, num_contigs=5,
            num_sites_per_contig=50):
        """ Validate that the reference sequences in the variant file matches
        a reference sequence file.
        """
        for contig in contigs[:num_contigs]:
            for var_data in list(vars_idx.fetch(contig))[:num_sites_per_contig]:
                ref_seq = aligner.seq(contig, var_data.start, var_data.stop)
                if ref_seq != var_data.ref:
                    # variant reference sequence does not match reference
                    logger = logging.get_logger('vars')
                    logger.debug((
                        'Reference sequence does not match variant reference ' +
                        'sequence at {} expected "{}" got "{}"').format(
                            ref_pos, var_data.ref, ref_seq))
                    return False

        return True

    def __init__(
            self, variant_fn, max_indel_size, all_paths,
            write_vars_txt, context_bases, vars_calib_fn=None,
            call_mode=DIPLOID_MODE, do_pr_ref_vars=False, aligner=None,
            keep_var_fp_open=False, do_validate_reference=True,
            edge_buffer=mh.DEFAULT_EDGE_BUFFER,
            context_min_alt_prob=mh.DEFAULT_CONTEXT_MIN_ALT_PROB,
            loc_index_in_memory=True):
        logger = logging.get_logger('vars')
        self.max_indel_size = max_indel_size
        self.all_paths = all_paths
        self.write_vars_txt = write_vars_txt
        self.vars_calib_fn = vars_calib_fn
        self.calib_table = calibration.VarCalibrator(self.vars_calib_fn)
        self.context_bases = context_bases
        if len(self.context_bases) != 2:
            raise mh.MegaError(
                'Must provide 2 context bases values (for single base ' +
                'variants and indels).')
        self.call_mode = call_mode
        self.do_pr_ref_vars = do_pr_ref_vars
        self.edge_buffer = edge_buffer
        self.context_min_alt_prob = context_min_alt_prob
        self.loc_index_in_memory = loc_index_in_memory
        self.variant_fn = variant_fn
        self.variants_idx = None
        if self.variant_fn is None:
            return

        logger.info('Loading variants.')
        vars_idx = pysam.VariantFile(self.variant_fn)
        try:
            contigs = list(vars_idx.header.contigs.keys())
            try:
                vars_idx.fetch(next(iter(contigs)), 0, 0)
            except StopIteration:
                logger.error('Variants file must contain contigs in header.')
                raise
        except ValueError:
            logger.warn(
                'Variants file must be indexed. Performing indexing now.')
            vars_idx.close()
            self.variant_fn = index_variants(self.variant_fn)
            vars_idx = pysam.VariantFile(self.variant_fn)
        if aligner is None:
            raise mh.MegaError(
                'Must provide aligner if variants filename is provided')
        if len(set(aligner.ref_names_and_lens[0]).intersection(contigs)) == 0:
            raise mh.MegaError((
                'Reference and variant files contain no chromosomes/contigs ' +
                'in common.\n\t\tFirst 3 reference contigs:\t{}\n\t\tFirst 3 ' +
                'variant file contigs:\t{}').format(
                    ', '.join(aligner.ref_names_and_lens[0][:3]),
                    ', '.join(contigs[:3])))
        if do_validate_reference and not self.check_vars_match_ref(
                vars_idx, contigs, aligner):
            raise mh.MegaError(
                'Reference sequence file does not match reference sequence ' +
                'from variants file.')

        if keep_var_fp_open:
            self.variants_idx = vars_idx
        else:
            vars_idx.close()
            self.variants_idx = None

        return

    @property
    def substitution_context(self):
        return self.context_bases[0]
    @property
    def indel_context(self):
        return self.context_bases[1]

    def reopen_variant_index(self):
        if self.variant_fn is not None:
            self.variants_idx = pysam.VariantFile(self.variant_fn)
        return

    @staticmethod
    def compute_variant_distance(var1, var2):
        # if the variants overlap return None
        if not (var1.start >= var2.stop or
                var2.start >= var1.stop):
            return None
        if var1.start >= var2.stop:
            return var1.start - var2.stop
        return var2.start - var1.stop

    @staticmethod
    def any_variants_overlap(variants):
        for var1, var2 in combinations(variants, 2):
            if not (var1.start >= var2.stop or
                    var2.start >= var1.stop):
                return True
        return False

    @staticmethod
    def annotate_context_seqs(
            context_vars, up_context_seq, dn_context_seq, context_ref_start,
            context_rel_var_start, context_rel_var_end):
        # annotate upstream sequence
        ann_up_seq = []
        prev_end = 0
        for context_var, np_alt_seq in context_vars:
            if context_var.stop - context_ref_start > context_rel_var_start:
                break
            ann_up_seq.extend((up_context_seq[
                prev_end:context_var.start - context_ref_start], np_alt_seq))
            prev_end = context_var.stop - context_ref_start
        ann_up_seq.append(up_context_seq[prev_end:])
        ann_up_seq = np.concatenate(ann_up_seq)

        # annotate downstream sequence
        ann_dn_seq = []
        prev_end = 0
        for context_var, np_alt_seq in context_vars:
            if context_var.start - context_ref_start < context_rel_var_end:
                continue
            ann_dn_seq.extend((dn_context_seq[
                prev_end:context_var.start - context_ref_start -
                context_rel_var_end], np_alt_seq))
            prev_end = (context_var.stop - context_ref_start -
                        context_rel_var_end)
        ann_dn_seq.append(dn_context_seq[prev_end:])
        ann_dn_seq = np.concatenate(ann_dn_seq)

        return ann_up_seq, ann_dn_seq

    @staticmethod
    def iter_variant_combos_by_distance(variant, context_variants):
        """ Yield combinations of variants ordered by inclusion of variants
        closer to the variant of interest first.
        """
        def iter_alt_variant_seqs(variants):
            """ Single variants can have multiple alternative sequences,
            so iterate over those here
            """
            vars_w_alts = [[(var, np_alt) for np_alt in var.np_alts]
                           for var in variants]
            for vars_w_alt in product(*vars_w_alts):
                yield sorted(vars_w_alt, key=lambda x: x[0].start)
            return

        dist_vars = defaultdict(list)
        for context_var in context_variants:
            var_dist = VarData.compute_variant_distance(variant, context_var)
            if var_dist is not None:
                dist_vars[var_dist].append(context_var)

        used_vars = []
        # sort list by distance to a variant in order to report most
        # proximal variants first
        for var_dist, dist_context_vars in sorted(dist_vars.items()):
            # loop over numbers of current distance variants (1 or more) and
            # previously used variants (0 or more)
            for n_dist_vars, n_used_vars in product(
                    range(1, len(dist_context_vars) + 1),
                    range(len(used_vars) + 1)):
                # loop over actual selection of variants
                for curr_vars, used_vars_i in product(
                        combinations(dist_context_vars, n_dist_vars),
                        combinations(used_vars, n_used_vars)):
                    # including selected alt seq
                    for vars_w_alt in iter_alt_variant_seqs(
                            list(curr_vars) + list(used_vars_i)):
                        yield vars_w_alt
            used_vars.extend(dist_context_vars)

        return

    @staticmethod
    def iter_context_variants(variants, context_max_dist):
        """ Iterate variants as well as variants within context_max_dist
        """
        vars_iter = iter(variants)
        def next_var_or_none():
            try:
                return next(vars_iter)
            except StopIteration:
                return None

        curr_vars = [next_var_or_none(),]
        next_var = next_var_or_none()
        if next_var is None:
            if curr_vars[0] is not None:
                yield curr_vars[0], []
            return
        curr_var_idx = 0

        while next_var is not None:
            if curr_var_idx == len(curr_vars):
                curr_vars.append(next_var)
                next_var = next_var_or_none()
            curr_var = curr_vars[curr_var_idx]
            # add relevant variants
            while (next_var is not None and
                   next_var.start - context_max_dist <= curr_var.stop):
                curr_vars.append(next_var)
                next_var = next_var_or_none()
            # remove variants that end before the variant of interest
            n_vars_removed = sum(
                var.stop + context_max_dist < curr_var.start + 1
                for var in curr_vars)
            curr_vars = [
                var for var in curr_vars
                if var.stop + context_max_dist >= curr_var.start + 1]
            curr_var_idx -= n_vars_removed - 1
            # yeild variants in range of current variant
            yield curr_var, curr_vars

        # yield final vars from the read
        while len(curr_vars) < curr_var_idx:
            curr_var = curr_vars[curr_var_idx]
            # remove variants that end before the variant of interest
            n_vars_removed = sum(
                var.stop + context_max_dist < curr_var.start
                for var in curr_vars)
            curr_vars = [
                var for var in curr_vars
                if var.stop + context_max_dist >= curr_var.start]
            curr_var_idx -= n_vars_removed - 1
            # yeild variants in range of current variant
            yield curr_var, curr_vars

        return

    @staticmethod
    def add_indel_context_base(
            np_ref_seq, np_alt_seqs, var_start, read_ref_fwd_seq,
            read_ref_pos):
        if np_ref_seq.shape[0] == 0 or any(
                np_alt.shape[0] == 0 for np_alt in np_alt_seqs):
            upstrm_base = mh.ALPHABET[read_ref_fwd_seq[
                var_start - read_ref_pos.start - 1]]
            var_start -= 1
        else:
            upstrm_base = ''
        ref_seq = upstrm_base + mh.int_to_seq(np_ref_seq)
        alt_seqs = tuple((upstrm_base + mh.int_to_seq(np_alt)
                          for np_alt in np_alt_seqs))
        return ref_seq, alt_seqs, var_start

    def merge_variants(
            self, grouped_read_vars, read_ref_fwd_seq, read_ref_pos):
        """ Merge atomized variants into multi-allelic sites.
        if this is not done, allele probabilities will not be normalized
        correctly.
        """
        variants = []
        for _, site_vars in sorted(grouped_read_vars.items()):
            if len(site_vars) == 1:
                # add upstream seq to simple indels
                (out_ref_seq, out_alt_seqs,
                 out_start) = self.add_indel_context_base(
                     site_vars[0].np_ref, site_vars[0].np_alts,
                     site_vars[0].start, read_ref_fwd_seq, read_ref_pos)
                site_var = site_vars[0]._replace(
                    ref=out_ref_seq, alts=out_alt_seqs, ref_start=out_start)
                variants.append(site_var)
                continue

            # join all valid ids
            # skip None ids ('.' in VCF)
            site_var_ids = set(
                var_id for var in site_vars
                if var.id is not None and var.id != '.'
                for var_id in var.id.split(';'))
            # if all ids are None leave id as None
            site_var_ids = (';'.join(sorted(site_var_ids))
                            if len(site_var_ids) > 0 else None)
            # assuming atomized variants with single alt
            # np arrays aren't hash-able, so test for equality manually here
            np_alts = []
            for var in site_vars:
                if any((len(var.np_alts[0]) == len(added_np_alt) and
                        np.all(var.np_alts[0] == added_np_alt))
                       for added_np_alt in np_alts):
                    continue
                np_alts.append(var.np_alts[0])

            # add upstream seq to simple indels
            out_ref_seq, out_alt_seqs, out_start = self.add_indel_context_base(
                site_vars[0].np_ref, np_alts, site_vars[0].start,
                read_ref_fwd_seq, read_ref_pos)

            variants.append(VARIANT_DATA(
                np_ref=site_vars[0].np_ref, np_alts=np_alts,
                id=site_var_ids, chrom=site_vars[0].chrom,
                start=site_vars[0].start, stop=site_vars[0].stop,
                ref=out_ref_seq, alts=out_alt_seqs, ref_start=out_start))

        return variants

    @staticmethod
    def expand_ambig_variant(
            np_ref_seq, np_alt_seq, var_start, read_ref_fwd_seq, read_ref_pos):
        # don't try to expand complex variants
        if np_ref_seq.shape[0] != 0 and np_alt_seq.shape[0] != 0:
            return np_ref_seq, np_alt_seq, var_start
        elif np_ref_seq.shape[0] == 0:
            # expand ambiguous insertion sequence
            indel_seq = np_alt_seq
        else:
            # expand ambiguous deletion sequence
            indel_seq = np_ref_seq

        expand_dnstrm_seq = []
        expand_start = expand_pos = (
            var_start - read_ref_pos.start + np_ref_seq.shape[0])
        # if shifting the insertion one base downstream would be
        # equivalent add this position
        while (expand_pos < read_ref_fwd_seq.shape[0] and
               indel_seq[(expand_pos - expand_start) %
                          indel_seq.shape[0]] ==
               read_ref_fwd_seq[expand_pos]):
            expand_dnstrm_seq.append(read_ref_fwd_seq[expand_pos])
            expand_pos += 1
        if expand_pos == read_ref_fwd_seq.shape[0]:
            raise mh.MegaError(
                'Variant is ambiguous up to the end of the read.')

        # mirror for upstream ambiguous sequence
        expand_upstrm_seq = []
        expand_start = expand_pos = var_start - read_ref_pos.start - 1
        # if shifting the insertion one base upstream would be
        # equivalent add this position
        while (expand_pos >= 0 and
               indel_seq[indel_seq.shape[0] - 1 - (
                   (expand_start - expand_pos) % indel_seq.shape[0])] ==
               read_ref_fwd_seq[expand_pos]):
            expand_upstrm_seq.insert(0, read_ref_fwd_seq[expand_pos])
            expand_pos -= 1
        if expand_pos == -1:
            raise mh.MegaError(
                'Variant is ambiguous up to the start of the read.')

        np_ref_seq = np.concatenate([
            expand_upstrm_seq, np_ref_seq, expand_dnstrm_seq]).astype(np.uintp)
        np_alt_seq = np.concatenate([
            expand_upstrm_seq, np_alt_seq, expand_dnstrm_seq]).astype(np.uintp)
        var_start -= len(expand_upstrm_seq)
        return np_ref_seq, np_alt_seq, var_start

    def iter_atomized_variants(
            self, var, np_ref_seq, np_alt_seq, read_ref_fwd_seq, read_ref_pos):
        # substitutions
        if np_alt_seq.shape[0] == np_ref_seq.shape[0]:
            # convert all substitutions into single base substitutions
            for sub_offset, (np_ref_base, np_alt_base) in enumerate(zip(
                    np_ref_seq, np_alt_seq)):
                if np_ref_base == np_alt_base: continue
                yield (
                    (var.start + sub_offset, var.start + sub_offset + 1),
                    VARIANT_DATA(
                        np_ref=np.array([np_ref_base], dtype=np.uintp),
                        np_alts=(np.array([np_alt_base], dtype=np.uintp),),
                        id=var.id, chrom=var.chrom,
                        start=var.start + sub_offset,
                        stop=var.start + sub_offset + 1))
        else:
            # skip large indels
            if np.abs(np_ref_seq.shape[0] -
                      np_alt_seq.shape[0]) > self.max_indel_size:
                return

            # trim context bases from seq
            np_ref_seq, np_alt_seq, start_trim, _ = simplify_var_seq(
                np_ref_seq, np_alt_seq)
            var_start = var.start + start_trim
            try:
                # expand seqs to include ambiguous locations
                np_ref_seq, np_alt_seq, var_start = self.expand_ambig_variant(
                    np_ref_seq, np_alt_seq, var_start, read_ref_fwd_seq,
                    read_ref_pos)
            except mh.MegaError:
                # if variant is ambiguous to the end of the read, then skip it
                return

            yield ((var_start, var_start + np_ref_seq.shape[0]),
                   VARIANT_DATA(
                       np_ref=np_ref_seq, np_alts=(np_alt_seq,), id=var.id,
                       chrom=var.chrom, start=var_start,
                       stop=var_start + np_ref_seq.shape[0]))
        return

    def atomize_variants(self, fetch_res, read_ref_fwd_seq, read_ref_pos):
        grouped_read_vars = defaultdict(list)
        for var in fetch_res:
            # fetch results include any overlap where only inclusive overlap
            # are valid here
            if (var.stop + self.edge_buffer > read_ref_pos.end or
                var.start - self.edge_buffer < read_ref_pos.start):
                continue

            np_ref_seq = mh.seq_to_int(var.ref)
            for alt_seq in var.alts:
                np_alt_seq = mh.seq_to_int(alt_seq)
                for start_stop, atom_var in self.iter_atomized_variants(
                        var, np_ref_seq, np_alt_seq, read_ref_fwd_seq,
                        read_ref_pos):
                    grouped_read_vars[start_stop].append(atom_var)
        return grouped_read_vars

    def fetch_read_variants(self, read_ref_pos, read_ref_fwd_seq):
        try:
            fetch_res = self.variants_idx.fetch(
                read_ref_pos.chrm, read_ref_pos.start + self.edge_buffer,
                read_ref_pos.end - self.edge_buffer)
        except ValueError:
            raise mh.MegaError('Mapped location not valid for variants file.')

        grouped_read_vars = self.atomize_variants(
            fetch_res, read_ref_fwd_seq, read_ref_pos)
        read_variants = self.merge_variants(
            grouped_read_vars, read_ref_fwd_seq, read_ref_pos)
        return read_variants

    def iter_vars(
            self, read_variants, read_ref_pos, read_ref_fwd_seq,
            max_contexts=16, context_max_dist=mh.CONTEXT_MAX_DIST):
        """Iterator over variants overlapping the read mapped position.

        Args:
            read_variants: List of variant objects (from fetch_read_variants)
            read_ref_pos: Read reference mapping position
                (`megalodon.mapping.MAP_POS`)
            read_ref_fwd_seq: Mapped reference sequence. Forward strand sequence
                no matter the mapped strand.
            max_contexts: Maximum number of context variant combinations to
                include around each variant.

        Yields:
            var_ref_seq: Reference variant sequence on read strand
            var_alt_seqs: Alternative variant sequences on read strand
            context_seqs: Sequences surrounding the variant on read strand
            context_start: Start of variant context in read coordinates
            context_end: End of variant context in read coordinates
            variant_ref: Reference variant sequence on reference
            variant_alts: Alternate variant sequences on reference
            variant_id: string idnentifier for the variant
            pos: variant position (0-based coordinate)

        Variantss within edge buffer of the end of the mapping will be ignored.

        If more than max_contexts variantss exist within context_basss then only
        the max_contexts most proximal to the variant in question will be
        returned.
        """
        def revcomp_variant(context_seqs, np_var_ref, var_alts):
            rc_context_seqs = [
                (mh.revcomp_np(dn_context_seq), mh.revcomp_np(up_context_seq))
                for up_context_seq, dn_context_seq in context_seqs]
            return rc_context_seqs, mh.revcomp_np(np_var_ref), [
                mh.revcomp_np(np_alt) for np_alt in np_var_alts]

        def extract_variant_contexts(variant, context_vars):
            # compute various relative coordinates (in reference coordinates
            # not on read strand, to make it easier to work with variants)
            # select single base substitution or indel context width
            var_context_bases = self.substitution_context if all(
                variant.np_ref.shape[0] == np_alt.shape[0]
                for np_alt in variant.np_alts) else self.indel_context
            context_ref_start = variant.start - var_context_bases
            if context_ref_start < read_ref_pos.start:
                context_ref_start = read_ref_pos.start
            context_ref_end = variant.stop + var_context_bases
            if context_ref_end > read_ref_pos.end:
                context_ref_end = read_ref_pos.end
            context_read_start = context_ref_start - read_ref_pos.start
            context_read_end = context_ref_end - read_ref_pos.start

            context_rel_var_start = variant.start - context_ref_start
            context_rel_var_end = (
                context_rel_var_start + variant.np_ref.shape[0])

            context_ref_seq = read_ref_fwd_seq[
                context_read_start:context_read_end]
            # first context is always reference sequence
            up_context_seq = context_ref_seq[:context_rel_var_start]
            dn_context_seq = context_ref_seq[context_rel_var_end:]
            context_seqs = [(up_context_seq, dn_context_seq),]

            if max_contexts == 1 or len(context_vars) == 0:
                return (context_ref_start, context_read_start, context_read_end,
                        context_seqs)

            for context_vars_i in self.iter_variant_combos_by_distance(
                    variant, context_vars):
                if self.any_variants_overlap(list(zip(*context_vars_i))[0]):
                    continue
                context_seqs.append(self.annotate_context_seqs(
                    context_vars_i, up_context_seq, dn_context_seq,
                    context_ref_start, context_rel_var_start,
                    context_rel_var_end))
                if len(context_seqs) >= max_contexts:
                    break

            return (context_ref_start, context_read_start, context_read_end,
                    context_seqs)


        logger = logging.get_logger('vars')
        for variant, context_variants in self.iter_context_variants(
                read_variants, context_max_dist):
            (context_ref_start, context_read_start, context_read_end,
             np_context_seqs) = extract_variant_contexts(
                 variant, context_variants)

            # revcomp seqs for strand and convert to numpy arrays
            np_var_ref, np_var_alts = variant.np_ref, variant.np_alts
            if read_ref_pos.strand == -1:
                np_context_seqs, np_var_ref, np_var_alts = revcomp_variant(
                    np_context_seqs, np_var_ref, np_var_alts)
                context_read_start, context_read_end = (
                    read_ref_fwd_seq.shape[0] - context_read_end,
                    read_ref_fwd_seq.shape[0] - context_read_start)
            if np.concatenate([np_var_ref,] + list(np_var_alts) + [
                    seq for cntxt_seqs in np_context_seqs
                    for seq in cntxt_seqs]).max() > len(mh.ALPHABET):
                # some sequence contained invalid characters
                logger.debug(
                    'Invalid sequence encountered for variant ' +
                    '"{}" at {}:{}'.format(
                        variant.id, variant.chrom, variant.start))
                continue

            yield (
                np_var_ref, np_var_alts, np_context_seqs,
                context_read_start, context_read_end, variant)

        return

    def calibrate_llr(self, llr, var_ref_seq, var_alt_seq):
        return self.calib_table.calibrate_llr(
            llr, var_ref_seq, var_alt_seq)


######################
##### VCF Writer #####
######################

class Variant(object):
    """ Variant for entry into VcfWriter.
    Currently only handles a single sample.
    """
    def __init__(
            self, chrom, pos, ref, alts, id='.', qual='.', filter='.',
            info=None, sample_dict=None):
        self.chrom = chrom
        self.pos = int(pos)
        self.ref = ref.upper()
        self.alts = [alt.upper() for alt in alts]
        self.alt = ','.join(self.alts)
        self.id = '' if id is None else str(id)
        self.qual = qual
        self.filter = str(filter)
        if info is None:
            info = {}
        self.info_dict = info
        if sample_dict is None:
            sample_dict = OrderedDict()
        self.sample_dict = sample_dict
        return

    @property
    def _sorted_format_keys(self):
        sorted_keys = sorted(self.sample_dict.keys())
        if 'GT' in sorted_keys:
            sorted_keys.insert(0, sorted_keys.pop(sorted_keys.index('GT')))
        if 'LOG_PROBS' in sorted_keys:
            # move log probs to end of format field for easier human readability
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

    def add_haploid_probs(self, probs, gts):
        # phred scaled likelihoods
        with np.errstate(divide='ignore'):
            gl = np.maximum(mh.MIN_GL_VALUE, np.log10(probs))
        raw_pl = -10 * gl
        # "normalized" PL values stored as decsribed by VCF format
        # abs to remove negative 0 from file
        pl = np.abs(np.minimum(raw_pl - raw_pl.min(), mh.MAX_PL_VALUE))
        s_pl = np.sort(pl)

        # add sample tags
        self.add_sample_field('GT', gts[np.argmax(probs)])
        try:
            qual = int(np.around(np.minimum(raw_pl[0], mh.MAX_PL_VALUE)))
        except ValueError:
            logger = logging.get_logger('vars')
            logger.debug(
                'NAN quality value encountered. gts:{}, probs:{}'.format(
                    str(gts), str(probs)))
            qual = mh.MAX_PL_VALUE
        self.qual = '{:.0f}'.format(qual) if qual > 0 else '.'
        self.add_sample_field('GQ', '{:.0f}'.format(np.around(s_pl[1])))
        self.add_sample_field(
            'GL', ','.join('{:.2f}' for _ in range(probs.shape[0])).format(*gl))
        self.add_sample_field(
            'PL', ','.join('{:.0f}' for _ in range(probs.shape[0])).format(
                *np.around(pl)))
        return

    def add_diploid_probs(self, probs, gts):
        # phred scaled likelihoods
        with np.errstate(divide='ignore'):
            gl = np.maximum(mh.MIN_GL_VALUE, np.log10(probs))
        raw_pl = -10 * gl
        # "normalized" PL values stored as decsribed by VCF format
        pl = np.minimum(raw_pl - raw_pl.min(), mh.MAX_PL_VALUE)
        s_pl = np.sort(pl)

        # add sample tags
        self.add_sample_field('GT', gts[np.argmax(probs)])
        try:
            qual = int(np.minimum(np.around(raw_pl[0]), mh.MAX_PL_VALUE))
        except ValueError:
            logger = logging.get_logger('vars')
            logger.debug(
                'NAN quality value encountered. gts:{}, probs:{}'.format(
                    str(gts), str(probs)))
            qual = mh.MAX_PL_VALUE
        self.qual = '{:.0f}'.format(qual) if qual > 0 else '.'
        self.add_sample_field('GQ', '{:.0f}'.format(np.around(s_pl[1])))
        self.add_sample_field(
            'GL', ','.join('{:.2f}' for _ in range(probs.shape[0])).format(*gl))
        self.add_sample_field(
            'PL', ','.join('{:.0f}' for _ in range(probs.shape[0])).format(
                *np.around(pl)))
        return

    def __eq__(self, var2):
        return (self.chrm, self.pos, self.id) == (var2.chrm, var2.pos, var2.id)
    def __ne__(self, var2):
        return (self.chrm, self.pos, self.id) != (var2.chrm, var2.pos, var2.id)
    def __lt__(self, var2):
        return (self.chrm, self.pos, self.id) < (var2.chrm, var2.pos, var2.id)
    def __le__(self, var2):
        return (self.chrm, self.pos, self.id) <= (var2.chrm, var2.pos, var2.id)
    def __gt__(self, var2):
        return (self.chrm, self.pos, self.id) > (var2.chrm, var2.pos, var2.id)
    def __ge__(self, var2):
        return (self.chrm, self.pos, self.id) >= (var2.chrm, var2.pos, var2.id)


class VcfWriter(object):
    """ VCF writer class
    """
    version_options = set(['4.1',])
    def __init__(
            self, filename, mode='w',
            header=('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
                    'INFO', 'FORMAT', SAMPLE_NAME),
            extra_meta_info=FIXED_VCF_MI, version='4.1', ref_fn=None,
            ref_names_and_lens=None, write_vcf_lp=False):
        self.filename = filename
        self.mode = mode
        self.header = header
        if version not in self.version_options:
            raise ValueError('version must be one of {}'.format(
                self.version_options))
        self.version = version
        contig_mis = [] if ref_names_and_lens is None else [
            mh.CONTIG_MI.format(ref_name, ref_len)
            for ref_name, ref_len in zip(*ref_names_and_lens)]
        self.meta = [
            mh.VCF_VERSION_MI.format(self.version),
            mh.FILE_DATE_MI.format(datetime.date.today().strftime("%Y%m%d")),
            mh.SOURCE_MI.format(MEGALODON_VERSION),
            mh.REF_MI.format(ref_fn)] + contig_mis + extra_meta_info
        if write_vcf_lp:
            self.meta.append(FORMAT_LOG_PROB_MI)

        self.handle = open(self.filename, self.mode, encoding='utf-8')
        self.handle.write('\n'.join('##' + line for line in self.meta) + '\n')
        self.handle.write('#' + '\t'.join(self.header) + '\n')
        return

    def close(self):
        self.handle.close()
        return

    def write_variant(self, variant):
        elements = [getattr(variant, field.lower()) for field in self.header]
        elements = ['.' if e == '' else e for e in elements]
        # VCF POS field is 1-based
        elements[self.header.index('POS')] += 1
        self.handle.write('{}\n'.format('\t'.join(map(str, elements))))
        self.handle.flush()

        return


#####################################
##### Variant Aggregation Class #####
#####################################

class AggVars(mh.AbstractAggregationClass):
    """ Class to assist in database queries for per-site aggregation of
    variant calls over reads.
    """
    def __init__(
            self, vars_db_fn, write_vcf_log_probs=False,
            load_in_mem_indices=True):
        # open as read only database
        if load_in_mem_indices:
            self.vars_db = VarsDb(vars_db_fn)
        else:
            self.vars_db = VarsDb(vars_db_fn, chrm_index_in_memory=False,
                                  alt_index_in_memory=False)
        self.n_uniq_vars = None
        self.write_vcf_log_probs = write_vcf_log_probs
        return

    def num_uniq(self):
        if self.n_uniq_vars is None:
            self.n_uniq_vars = self.vars_db.get_num_uniq_var_loc()
        return self.n_uniq_vars

    def iter_uniq(self):
        for q_val in self.vars_db.iter_locs():
            yield q_val
        return

    def compute_diploid_probs(self, ref_lps, alts_lps, het_factor=1.0):
        def compute_het_lp(a1, a2):
            # order by the inverse log likelihood ratio
            llr_ord = np.argsort(all_lps[a2] - all_lps[a1])
            s_a1_lps = all_lps[a1][llr_ord]
            s_a2_lps = all_lps[a2][llr_ord]
            with np.errstate(divide='ignore'):
                # compute log probability of heterozygous genotype by binomial
                # weighted sum of maximum likelihoods
                het_lp = logsumexp(np.array([
                    np.log(stats.binom.pmf(i, all_lps.shape[1], 0.5)) +
                    np.sum(s_a1_lps[:i]) + np.sum(s_a2_lps[i:])
                    for i in range(all_lps.shape[1] + 1)]))
            return het_lp


        all_lps = np.concatenate([ref_lps.reshape(1, -1), alts_lps], axis=0)
        genotype_lps, het_gts, gts = [], [], []
        # order genotypes as described here under the GL genotype fields section
        # http://samtools.github.io/hts-specs/VCFv4.1.pdf
        for a2 in range(all_lps.shape[0]):
            for a1 in range(a2 + 1):
                gts.append('{}/{}'.format(a1, a2))
                if a1 == a2:
                    genotype_lps.append(np.sum(all_lps[a1]))
                    het_gts.append(False)
                else:
                    genotype_lps.append(compute_het_lp(a1, a2))
                    het_gts.append(True)

        log_prior_weights = np.array([
            0.0 if het_gt else all_lps.shape[1] * np.log(het_factor)
            for het_gt in het_gts])
        log_prior_weights = log_prior_weights - logsumexp(log_prior_weights)
        var_lps = np.array(genotype_lps) + log_prior_weights
        post_var_lps = var_lps - logsumexp(var_lps)
        return np.exp(post_var_lps), gts

    def compute_haploid_probs(self, ref_lps, alts_lps):
        var_lps = np.concatenate([[ref_lps.sum()], alts_lps.sum(axis=1)])
        post_var_lps = var_lps - logsumexp(var_lps)
        return np.exp(post_var_lps), list(map(str, range(var_lps.shape[0])))

    def compute_var_stats(
            self, var_loc, het_factors, call_mode=DIPLOID_MODE,
            valid_read_ids=None):
        assert call_mode in (HAPLIOD_MODE, DIPLOID_MODE), (
            'Invalid variant aggregation ploidy call mode: {}.'.format(
                call_mode))

        pr_var_stats = self.vars_db.get_loc_stats(
            var_loc, valid_read_ids is not None)
        alt_seqs = sorted(set(r_stats.alt_seq for r_stats in pr_var_stats))
        pr_alt_lps = defaultdict(dict)
        for r_stats in pr_var_stats:
            if (valid_read_ids is not None and
                r_stats.read_id not in valid_read_ids):
                continue
            pr_alt_lps[r_stats.read_id][r_stats.alt_seq] = r_stats.score
        if len(pr_alt_lps) == 0:
            raise mh.MegaError('No valid reads cover variant')

        alt_seq_lps = [[] for _ in range(len(alt_seqs))]
        for read_lps in pr_alt_lps.values():
            for i, alt_seq in enumerate(alt_seqs):
                try:
                    alt_seq_lps[i].append(read_lps[alt_seq])
                except KeyError:
                    raise mh.MegaError(
                        'Alternative variant seqence must exist for all reads.')
        alts_lps = np.stack(alt_seq_lps, axis=0)
        with np.errstate(all='ignore'):
            ref_lps = np.log1p(-np.exp(alts_lps).sum(axis=0))

        r0_stats = pr_var_stats[0]
        variant = Variant(
            chrom=r0_stats.chrm, pos=r0_stats.pos, ref=r0_stats.ref_seq,
            alts=alt_seqs, id=r0_stats.var_name)
        variant.add_tag('DP', '{}'.format(ref_lps.shape[0]))
        variant.add_sample_field('DP', '{}'.format(ref_lps.shape[0]))

        if self.write_vcf_log_probs:
            variant.add_sample_field('LOG_PROBS', ','.join(
                ';'.join('{:.2f}'.format(lp) for lp in alt_i_lps)
                for alt_i_lps in alts_lps))

        if call_mode == DIPLOID_MODE:
            het_factor = (
                het_factors[0] if len(r0_stats.ref_seq) == 1 and
                len(r0_stats.alt_seq) == 1 else
                het_factors[1])
            diploid_probs, gts = self.compute_diploid_probs(
                ref_lps, alts_lps, het_factor)
            variant.add_diploid_probs(diploid_probs, gts)
        elif call_mode == HAPLIOD_MODE:
            haploid_probs, gts = self.compute_haploid_probs(ref_lps, alts_lps)
            variant.add_haploid_probs(haploid_probs, gts)

        return variant

    def close(self):
        self.vars_db.close()
        return


#########################################
##### Whatshap Mapping Post-process #####
#########################################

def get_whatshap_command(index_variant_fn, whatshap_sort_fn, phase_fn):
    return ('Run following command to obtain phased variants:\n\t\t' +
            'whatshap phase --indels --distrust-genotypes -o {} {} {}').format(
                phase_fn, index_variant_fn, whatshap_sort_fn)

def sort_variants(in_vcf_fn, out_vcf_fn):
    in_vcf_fp = pysam.VariantFile(in_vcf_fn)
    with pysam.VariantFile(
            out_vcf_fn, 'w', header=in_vcf_fp.header) as out_vcf_fp:
        for rec in sorted(in_vcf_fp.fetch(), key=lambda r: (r.chrom, r.start)):
            out_vcf_fp.write(rec)
    return

def index_variants(variant_fn):
    try:
        return pysam.tabix_index(
            variant_fn, force=True, preset='vcf', keep_original=True)
    except OSError:
        # file likely not sorted
        sort_variant_fn = mh.add_fn_suffix(variant_fn, 'sorted')
        sort_variants(variant_fn, sort_variant_fn)
        return pysam.tabix_index(
            sort_variant_fn, force=True, preset='vcf', keep_original=True)


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
