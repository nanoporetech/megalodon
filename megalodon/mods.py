import re
import sys
import queue
import sqlite3
import datetime
from time import sleep
from collections import defaultdict, namedtuple, OrderedDict

import numpy as np

from megalodon import (
    calibration, decode, logging, mapping, megalodon_helper as mh)
from megalodon._version import MEGALODON_VERSION


BIN_THRESH_NAME = 'binary_threshold'
# EM method is broken from transition to multi-mod sites support. Need to update
# function accordingly
EM_NAME = 'em'
AGG_METHOD_NAMES = set((BIN_THRESH_NAME,))
AGG_INFO = namedtuple('AGG_INFO', ('method', 'binary_threshold'))
DEFAULT_AGG_INFO = AGG_INFO(BIN_THRESH_NAME, 0.75)

FIXED_VCF_MI = [
    'INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
    'INFO=<ID=SN,Number=1,Type=String,Description="Strand">',
    'FORMAT=<ID=VALID_DP,Number=1,Type=Integer,Description="Valid Read Depth">',
]
MOD_MI_TMPLTS = [
    'FORMAT=<ID={0},Number=1,Type=Float,Description=' +
    '"{1} Modified Base Proportion">']
FORMAT_LOG_PROB_MI = (
    'FORMAT=<ID=LOG_PROBS,Number=A,Type=String,' +
    'Description="Per-read log10 likelihoods for modified ' +
    'bases (semi-colon separated)">')

OUT_BUFFER_LIMIT = 10000

# allow 64GB for memory mapped sqlite file access
MEMORY_MAP_LIMIT = 64000000000


###################
##### Mods DB #####
###################

class ModsDb(object):
    # note foreign key constraint is not applied here as this
    # drastically reduces efficiency. Namely the search for pos_id
    # when inserting into the data table forces a scan of a very
    # large table or maintainance of a very large pos table index
    # both of which slow data base speed
    # thus foreign constraint must be handled by the class
    chrm_tbl = OrderedDict((
        ('chrm_id', 'INTEGER PRIMARY KEY'),
        ('chrm', 'TEXT')))
    pos_tbl = OrderedDict((
        ('pos_id', 'INTEGER PRIMARY KEY'),
        ('pos_chrm', 'INTEGER'),
        ('strand', 'INTEGER'),
        ('pos', 'INTEGER')))
    mod_tbl = OrderedDict((
        ('mod_id', 'INTEGER PRIMARY KEY'),
        ('mod_base', 'TEXT'),
        ('motif', 'TEXT'),
        ('motif_pos', 'INTEGER'),
        ('raw_motif', 'TEXT')))
    read_tbl = OrderedDict((
        ('read_id', 'INTEGER PRIMARY KEY'),
        ('uuid', 'TEXT')))
    data_tbl = OrderedDict((
        ('score', 'FLOAT'),
        ('score_pos', 'INTEGER'),
        ('score_mod', 'INTEGER'),
        ('score_read', 'INTEGER')))

    # namedtuple for returning mods from a single position
    mod_data = namedtuple('mod_data', [
        'read_id', 'chrm', 'strand', 'pos', 'score', 'mod_base', 'motif',
        'motif_pos', 'raw_motif'])

    def __init__(self, fn, read_only=True, db_safety=1,
                 pos_index_in_memory=False, mod_chrm_index_in_memory=True):
        """ Interface to database containing modified base statistics.

        Default settings are for optimal read_only performance.
        """
        self.fn = mh.resolve_path(fn)
        self.read_only = read_only
        self.pos_idx_in_mem = pos_index_in_memory
        self.cm_idx_in_mem = mod_chrm_index_in_memory

        if read_only:
            self.db = sqlite3.connect('file:' + fn + '?mode=ro', uri=True)
        else:
            self.db = sqlite3.connect(fn)

        self.cur = self.db.cursor()
        if read_only:
            # use memory mapped file access
            self.db.execute('PRAGMA mmap_size = {}'.format(MEMORY_MAP_LIMIT))
            if self.cm_idx_in_mem:
                self.load_chrm_read_index()
                self.load_mod_read_index()
            if self.pos_idx_in_mem:
                self.load_pos_index()
        else:
            if db_safety < 2:
                # set asynchronous mode to off for max speed
                self.db.execute('PRAGMA synchronous = OFF')
            if db_safety < 1:
                # set no rollback mode
                self.db.execute('PRAGMA journal_mode = OFF')

            # create tables
            for tbl_name, tbl in (
                    ('chrm', self.chrm_tbl), ('pos', self.pos_tbl),
                    ('mod', self.mod_tbl), ('read', self.read_tbl),
                    ('data', self.data_tbl)):
                try:
                    self.db.execute("CREATE TABLE {} ({})".format(
                        tbl_name, ','.join((
                            '{} {}'.format(*ft) for ft in tbl.items()))))
                except sqlite3.OperationalError:
                    raise mh.MegaError(
                        'Modified bases database already exists. Either ' +
                        'provide location for new database or open in ' +
                        'read_only mode.')

            if self.pos_idx_in_mem:
                self.pos_idx = {}
            else:
                self.create_pos_index()
            if self.cm_idx_in_mem:
                self.chrm_idx = {}
                self.mod_idx = {}
            else:
                self.create_mod_index()
                self.create_chrm_index()

        return

    def insert_chrm(self, chrm):
        self.cur.execute('INSERT INTO chrm (chrm) VALUES (?)', (chrm,))
        if self.cm_idx_in_mem:
            self.chrm_idx[chrm] = self.cur.lastrowid
        return self.cur.lastrowid

    def get_chrm_id(self, chrm):
        try:
            if self.cm_idx_in_mem:
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
            chrm = self.cur.execute(
                'SELECT chrm FROM chrm WHERE chrm_id=?',
                (chrm_id,)).fetchone()[0]
        except TypeError:
            raise mh.MegaError('Reference record (chromosome) not found in ' +
                               'mods database.')
        return chrm

    def get_pos_id(self, chrm, strand, pos, chrm_id=None):
        if chrm_id is None:
            chrm_id = self.get_chrm_id(chrm)

        try:
            if self.pos_idx_in_mem:
                pos_id = self.pos_idx[(chrm_id, strand, pos)]
            else:
                pos_id = self.cur.execute(
                    'SELECT pos_id FROM pos WHERE pos_chrm=? AND strand=? ' +
                    'AND pos=?', (chrm_id, strand, pos)).fetchone()[0]
        except (TypeError, KeyError):
            raise mh.MegaError(
                'Reference position not found in database.')

        return pos_id

    def get_pos_id_or_insert(self, chrm, strand, pos, chrm_id=None):
        if chrm_id is None:
            chrm_id = self.get_chrm_id(chrm)
        try:
            pos_id = self.get_pos_id(chrm, strand, pos, chrm_id)
        except mh.MegaError:
            self.cur.execute(
                'INSERT INTO pos (pos_chrm, strand, pos) VALUES (?,?,?)',
                (chrm_id, strand, pos))
            pos_id = self.cur.lastrowid
            if self.pos_idx_in_mem:
                self.pos_idx[(chrm_id, strand, pos)] = pos_id
        return pos_id

    def get_mod_base_id(self, mod_base, motif, motif_pos, raw_motif):
        try:
            if self.cm_idx_in_mem:
                mod_id = self.mod_idx[(mod_base, motif, motif_pos, raw_motif)]
            else:
                mod_id = self.cur.execute(
                    'SELECT mod_id FROM mod WHERE mod_base=? AND motif=? AND ' +
                    'motif_pos=? AND raw_motif=?',
                    (mod_base, motif, motif_pos, raw_motif)).fetchone()[0]
        except (TypeError, KeyError):
            raise mh.MegaError('Modified base not found in mods database.')
        return mod_id

    def get_mod_base_id_or_insert(self, mod_base, motif, motif_pos, raw_motif):
        try:
            mod_base_id = self.get_mod_base_id(
                mod_base, motif, motif_pos, raw_motif)
        except mh.MegaError:
            self.cur.execute(
                'INSERT INTO mod (mod_base, motif, motif_pos, raw_motif) ' +
                'VALUES (?,?,?,?)', (mod_base, motif, motif_pos, raw_motif))
            mod_base_id = self.cur.lastrowid
            if self.cm_idx_in_mem:
                self.mod_idx[(mod_base, motif, motif_pos,
                              raw_motif)] = mod_base_id
        return mod_base_id

    def insert_read_scores(self, r_mod_scores, uuid, chrm, strand):
        self.cur.execute('INSERT INTO read (uuid) VALUES (?)', (uuid,))
        read_id = self.cur.lastrowid
        chrm_id = self.get_chrm_id(chrm)
        read_insert_data = []
        for (pos, mod_lps, mod_bases, ref_motif, rel_pos,
             raw_motif) in r_mod_scores:
            pos_id = self.get_pos_id_or_insert(None, strand, pos, chrm_id)
            for mod_lp, mod_base in zip(mod_lps, mod_bases):
                mod_base_id = self.get_mod_base_id_or_insert(
                    mod_base, ref_motif, rel_pos, raw_motif)
                read_insert_data.append((mod_lp, pos_id, mod_base_id, read_id))

        self.cur.executemany(
            'INSERT INTO data VALUES (?,?,?,?)', read_insert_data)
        return

    def create_chrm_index(self):
        self.cur.execute('CREATE UNIQUE INDEX chrm_idx ON chrm(chrm)')
        return

    def load_chrm_read_index(self):
        self.chrm_read_idx = {}
        self.cur.execute('SELECT chrm_id, chrm FROM chrm')
        for chrm_id, chrm in self.cur:
            self.chrm_read_idx[chrm_id] = chrm
        return

    def create_mod_index(self):
        self.cur.execute('CREATE UNIQUE INDEX mod_idx ON ' +
                         'mod(mod_base, motif, motif_pos, raw_motif)')
        return

    def load_mod_read_index(self):
        self.mod_read_idx = {}
        self.cur.execute(
            'SELECT mod_id, mod_base, motif, motif_pos, raw_motif FROM mod')
        for mod_id, mod_base, motif, motif_pos, raw_motif in self.cur:
            self.mod_read_idx[mod_id] = (mod_base, motif, motif_pos, raw_motif)
        return

    def create_pos_index(self):
        self.cur.execute('CREATE UNIQUE INDEX pos_idx ON pos' +
                         '(pos_chrm, strand, pos)')
        return

    def load_pos_index(self):
        self.pos_read_idx = {}
        self.cur.execute('SELECT pos_id, pos_chrm, strand, pos FROM pos')
        for pos_id, chrm_id, strand, pos in self.cur:
            self.pos_read_idx[pos_id] = (chrm_id, strand, pos)
        return

    def create_data_covering_index(self):
        self.cur.execute('CREATE INDEX data_cov_idx ON data(' +
                         'score_pos, score_mod, score_read, score)')
        return

    def close(self):
        self.db.commit()
        self.db.close()
        return

    def get_num_uniq_mod_pos(self):
        return self.cur.execute('SELECT MAX(pos_id) FROM pos').fetchone()[0]

    def iter_pos_id(self):
        self.cur.execute('SELECT pos_id FROM pos')
        for pos_id in self.cur:
            yield pos_id[0]

        return

    def iter_pos(self):
        self.cur.execute('SELECT pos_id, pos_chrm, strand, pos FROM pos')
        for pos in self.cur:
            yield pos

        return

    def iter_pos_id_ordered(self):
        self.cur.execute('SELECT pos_id FROM pos ORDER BY pos_id')
        for pos_id in self.cur:
            yield pos_id[0]

        return

    def iter_pos_ordered(self):
        self.cur.execute('SELECT pos_id, pos_chrm, strand, pos FROM pos ' +
                         'ORDER BY pos_id')
        for pos in self.cur:
            yield pos

        return

    def get_pos_stats(self, pos_data, return_uuids=False):
        pos_id, chrm_id, strand, pos = pos_data
        if self.cm_idx_in_mem:
            if return_uuids:
                self.cur.execute(
                    'SELECT read.uuid, data.score, data.score_mod ' +
                    'FROM data ' +
                    'INNER JOIN read ON read.read_id = data.score_read ' +
                    'WHERE score_pos=?', (pos_id, ))
            else:
                # simplest query using covering index
                self.cur.execute(
                    'SELECT score_read, score, score_mod ' +
                    'FROM data ' +
                    'WHERE score_pos=?', (pos_id, ))
            return [
                self.mod_data(read_id, self.chrm_read_idx[chrm_id], strand,
                              pos, score, *self.mod_read_idx[mod_id])
                for read_id, score, mod_id in self.cur]

        # perform full query from on-disk database
        chrm = self.cur.execute(
            'SELECT chrm FROM chrm WHERE chrm_id=?', (chrm_id,)).fetchone()[0]
        if return_uuids:
            self.cur.execute(
                'SELECT read.uuid, data.score, mod.mod_base, mod.motif, ' +
                'mod.motif_pos, mod.raw_motif ' +
                'FROM data ' +
                'INNER JOIN read ON read.read_id = data.score_read ' +
                'INNER JOIN mod ON mod.mod_id = data.score_mod ' +
                'WHERE score_pos=?', (pos_id, ))
        else:
            self.cur.execute(
                'SELECT data.score_read, data.score, mod.mod_base, ' +
                'mod.motif, mod.motif_pos, mod.raw_motif ' +
                'FROM data ' +
                'INNER JOIN mod ON mod.mod_id = data.score_mod ' +
                'WHERE score_pos=?', (pos_id, ))
        return [self.mod_data(read_id, chrm, strand, pos, score, mod_base,
                              motif, motif_pos, raw_motif)
                for read_id, score, mod_base, motif, motif_pos, raw_motif in
                self.cur]

    def get_read_id(self, uuid):
        try:
            read_id = self.cur.execute(
                'SELECT read_id FROM read WHERE uuid=?', (uuid,)).fetchone()[0]
        except TypeError:
            raise mh.MegaError('Read ID not found in mods data base.')
        return read_id

    def get_read_id_or_insert(self, uuid):
        try:
            read_id = self.get_read_id(uuid)
        except mh.MegaError:
            self.cur.execute('INSERT INTO read (uuid) VALUES (?)', (uuid,))
            read_id = self.cur.lastrowid
        return read_id

    def create_data_read_index(self):
        self.cur.execute('CREATE INDEX data_read_idx ON data(score_read)')
        return

    def get_read_stats(self, uuid):
        # TODO optimize this as with position query
        self.cur.execute(
            'SELECT uuid, chrm.chrm, pos.strand, pos.pos, data.score, ' +
            ' mod.mod_base, mod.motif, mod.motif_pos, mod.raw_motif ' +
            'FROM read ' +
            'INNER JOIN data ON data.score_read = read.uuid ' +
            'INNER JOIN pos ON pos.pos_id = data.score_pos ' +
            'INNER JOIN mod ON mod.mod_id = data.score_mod ' +
            'INNER JOIN chrm ON chrm.chrm_id = pos.pos_chrm ' +
            'WHERE uuid=?', (uuid,))
        return [self.mod_data(*pos_data_i) for pos_data_i in self.cur]


################################
##### Per-read Mod Scoring #####
################################

def score_mod_seq(
        tpost, seq, mod_cats, can_mods_offsets,
        tpost_start=0, tpost_end=None, all_paths=False):
    """Score a section of log transition posteriors against a proposed sequence
    using a global mapping.
    :param tpost: `ndarray` containing log transition posteriors to be scored
    :param seq: `ndarray` containing integers encoding proposed sequence
    :param mod_cats: `ndarray` containing integers encoding proposed modified
        base labels
    :param can_mods_offsets: `ndarray` containing integers encoding proposed
        modified base labels
    :param tpost_start: start position within post (Default: 0)
    :param tpost_end: end position within post (Default: full posterior)
    :param all_paths: boolean to produce the forwards all paths score
        (default Viterbi best path)
    """
    seq = seq.astype(np.uintp)
    if tpost_end is None:
        tpost_end = tpost.shape[0]

    return decode.score_mod_seq(
        tpost, seq, mod_cats, can_mods_offsets, tpost_start, tpost_end,
        all_paths)

def call_read_mods(
        r_ref_pos, r_ref_seq, rl_cumsum, r_to_q_poss, r_post,
        post_mapped_start, mods_info):
    def iter_motif_sites(r_ref_seq):
        max_pos = len(r_ref_seq) - mods_info.edge_buffer
        for motif, rel_pos, mod_bases, raw_motif in mods_info.all_mod_motifs:
            for motif_match in motif.finditer(r_ref_seq):
                m_pos = motif_match.start() + rel_pos
                if m_pos < mods_info.edge_buffer: continue
                if m_pos > max_pos: break
                yield m_pos, mod_bases, motif_match.group(), rel_pos, raw_motif
        return


    logger = logging.get_logger('mods')
    # call all mods overlapping this read
    r_mod_scores = []
    for (pos, mod_bases, ref_motif, rel_pos,
         raw_motif) in iter_motif_sites(r_ref_seq):
        pos_bb, pos_ab = min(mods_info.mod_context_bases, pos), min(
            mods_info.mod_context_bases, len(r_ref_seq) - pos - 1)
        try:
            pos_ref_seq = mh.seq_to_int(
                r_ref_seq[pos - pos_bb:pos + pos_ab + 1])
        except mh.MegaError:
            ref_pos = r_ref_pos.start + pos if r_ref_pos.strand == 1 else \
                      r_ref_pos.start + len(r_ref_seq) - pos - 1
            logger.debug(
                'Invalid sequence encountered calling modified base ' +
                'at {}:{}'.format(r_ref_pos.chrm, ref_pos))
            continue
        pos_can_mods = np.zeros_like(pos_ref_seq)

        blk_start, blk_end = (rl_cumsum[r_to_q_poss[pos - pos_bb]],
                              rl_cumsum[r_to_q_poss[pos + pos_ab]])
        if blk_end - blk_start < (mods_info.mod_context_bases * 2) + 1:
            # no valid mapping over large inserted query bases
            # i.e. need as many "events/strides" as bases for valid mapping
            continue

        loc_can_score = score_mod_seq(
            r_post, pos_ref_seq, pos_can_mods, mods_info.can_mods_offsets,
            post_mapped_start + blk_start, post_mapped_start + blk_end,
            mods_info.mod_all_paths)
        if loc_can_score is None:
            raise mh.MegaError('Score computation error (memory error)')

        calib_llrs = []
        for mod_base in mod_bases:
            pos_mod_mods = pos_can_mods.copy()
            pos_mod_mods[pos_bb] = mods_info.str_to_int_mod_labels[mod_base]
            loc_mod_score = score_mod_seq(
                r_post, pos_ref_seq, pos_mod_mods, mods_info.can_mods_offsets,
                post_mapped_start + blk_start, post_mapped_start + blk_end,
                mods_info.mod_all_paths)
            if loc_mod_score is None:
                raise mh.MegaError('Score computation error (memory error)')

            # calibrate llr scores
            calib_llrs.append(mods_info.calibrate_llr(
                loc_can_score - loc_mod_score, mod_base))

        # due to calibration mutli-mod log likelihoods could result in
        # inferred negative reference likelihood, so re-normalize here
        loc_mod_lps = calibration.compute_log_probs(np.array(calib_llrs))

        m_ref_pos = (pos + r_ref_pos.start if r_ref_pos.strand == 1 else
                     r_ref_pos.end - pos - 1)
        r_mod_scores.append((
            m_ref_pos, loc_mod_lps, mod_bases, ref_motif, rel_pos, raw_motif))

    return r_mod_scores


###############################
##### Per-read Mod Output #####
###############################

def annotate_mods(r_start, ref_seq, r_mod_scores, strand):
    """ Annotate reference sequence with called modified bases.

    Note: Reference sequence is in read orientation and mod calls are in
    genome coordiates.
    """
    mod_seqs = []
    prev_pos = 0
    if strand == -1:
        ref_seq = ref_seq[::-1]
    for mod_pos, mod_lps, mod_bases, _, _, _ in sorted(r_mod_scores):
        can_lp = np.log1p(-np.exp(mod_lps).sum())
        # called canonical
        if can_lp >= mod_lps.max(): continue
        most_prob_mod = np.argmax(mod_lps)
        mod_seqs.append(ref_seq[prev_pos:mod_pos - r_start] +
                        mod_bases[most_prob_mod])
        prev_pos = mod_pos - r_start + 1
    mod_seqs.append(ref_seq[prev_pos:])
    mod_seq = ''.join(mod_seqs)
    if strand == -1:
        mod_seq = mod_seq[::-1]

    return mod_seq

def _get_mods_queue(
        mods_q, mods_conn, mods_db_fn, db_safety, ref_names_and_lens,
        mods_txt_fn, pr_refs_fn, pr_ref_filts):
    def store_mod_call(
            r_mod_scores,  read_id, chrm, strand, r_start, ref_seq,
            read_len, q_st, q_en, cigar, been_warned):
        try:
            mods_db.insert_read_scores(r_mod_scores, read_id, chrm, strand)
        except Exception as e:
            if not been_warned:
                logger.warning(
                    'Error inserting modified base scores into database. See ' +
                    'log debug output for error details.')
                been_warned = True
            import traceback
            var = traceback.format_exc()
            logger.debug(
                'Error inserting modified base scores into database: ' +
                str(e) + '\n' + var)

        if mods_txt_fp is not None and len(r_mod_scores) > 0:
            # would involve batching and creating several conversion tables
            # for var strings (read_id and chrms).
            mods_txt_fp.write('\n'.join((
                ('\t'.join('{}' for _ in field_names)).format(
                    read_id, chrm, strand, pos, mod_lp,
                    np.log1p(-np.exp(mod_lps).sum()), mod_base,
                    '{}:{}'.format(raw_motif, rel_pos))
                for pos, mod_lps, mod_bases, ref_motif, rel_pos, raw_motif
                in r_mod_scores
                for mod_lp, mod_base in zip(mod_lps, mod_bases))) + '\n')
        if pr_refs_fn is not None:
            if not mapping.read_passes_filters(
                    pr_ref_filts, read_len, q_st, q_en, cigar):
                return

            pr_refs_fp.write('>{}\n{}\n'.format(read_id, annotate_mods(
                r_start, ref_seq, r_mod_scores, strand)))

        return been_warned


    logger = logging.get_logger('mods')
    been_warned = False

    mods_db = ModsDb(mods_db_fn, db_safety=db_safety, read_only=False,
                     pos_index_in_memory=True)
    for ref_name in ref_names_and_lens[0]:
        mods_db.insert_chrm(ref_name)
    mods_db.create_chrm_index()

    if mods_txt_fn is None:
        mods_txt_fp = None
    else:
        mods_txt_fp = open(mods_txt_fn, 'w')
        field_names = (
            'read_id', 'chrm', 'strand', 'pos', 'mod_log_prob',
            'can_log_prob', 'mod_base', 'motif')
        mods_txt_fp.write('\t'.join(field_names) + '\n')

    if pr_refs_fn is not None:
        pr_refs_fp = open(pr_refs_fn, 'w')

    while True:
        try:
            # note strand is +1 for fwd or -1 for rev
            r_mod_scores, (
                read_id, chrm, strand, r_start, ref_seq, read_len, q_st, q_en,
                cigar) = mods_q.get(block=False)
        except queue.Empty:
            if mods_conn.poll():
                break
            sleep(0.001)
            continue
        try:
            been_warned = store_mod_call(
                r_mod_scores,  read_id, chrm, strand, r_start, ref_seq,
                read_len, q_st, q_en, cigar, been_warned)
        except Exception as e:
            logger.debug('Error processing mods output for read: ' +
                         '{}\nError type: {}'.format(read_id, str(e)))

    while not mods_q.empty():
        r_mod_scores, (
            read_id, chrm, strand, r_start, ref_seq, read_len, q_st, q_en,
            cigar) = mods_q.get(block=False)
        try:
            been_warned = store_mod_call(
                r_mod_scores,  read_id, chrm, strand, r_start, ref_seq,
                read_len, q_st, q_en, cigar, been_warned)
        except Exception as e:
            logger.debug('Error processing mods output for read: ' +
                         '{}\nError type: {}'.format(read_id, str(e)))

    if mods_txt_fp is not None: mods_txt_fp.close()
    if pr_refs_fn is not None: pr_refs_fp.close()
    mods_db.create_mod_index()
    mods_db.create_data_covering_index()
    mods_db.close()

    return


####################
##### Mod Info #####
####################

class ModInfo(object):
    single_letter_code = {
        'A':'A', 'C':'C', 'G':'G', 'T':'T', 'B':'CGT',
        'D':'AGT', 'H':'ACT', 'K':'GT', 'M':'AC',
        'N':'ACGT', 'R':'AG', 'S':'CG', 'V':'ACG',
        'W':'AT', 'Y':'CT'}

    def distinct_bases(self, b1, b2):
        return len(set(self.single_letter_code[b1]).intersection(
            self.single_letter_code[b2])) == 0

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
                mods_bases = list(mod_bases)
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
                    '[{}]'.format(self.single_letter_code[letter])
                    for letter in raw_motif))
                self.all_mod_motifs.append((motif, pos, mod_bases, raw_motif))

            if not self.distinct_motifs():
                raise mh.MegaError(
                    'One provided motif can be found within another motif. ' +
                    'Only distinct sets of motifs are accepted')

        return

    def __init__(
            self, model_info, all_mod_motifs_raw=None, mod_all_paths=False,
            write_mods_txt=None, mod_context_bases=None,
            do_output_mods=False, do_pr_ref_mods=False, mods_calib_fn=None,
            mod_output_fmts=[mh.MOD_BEDMETHYL_NAME],
            edge_buffer=mh.DEFAULT_EDGE_BUFFER):
        logger = logging.get_logger()
        # this is pretty hacky, but these attributes are stored here as
        # they are generally needed alongside other alphabet info
        # don't want to pass all of these parameters around individually though
        # as this would make function signatures too complicated
        self.mod_all_paths = mod_all_paths
        self.write_mods_txt = write_mods_txt
        self.mod_context_bases = mod_context_bases
        self.do_output_mods = do_output_mods
        self.do_pr_ref_mods = do_pr_ref_mods
        self.mod_long_names = model_info.mod_long_names
        self.calib_table = calibration.ModCalibrator(mods_calib_fn)
        self.mod_output_fmts = mod_output_fmts
        self.edge_buffer = edge_buffer

        self.alphabet = model_info.can_alphabet
        self.ncan_base = len(self.alphabet)
        try:
            self.alphabet = self.alphabet.decode()
        except:
            pass
        if model_info.is_cat_mod:
            # TODO also output "(alt to C)" for each mod
            logger.info(
                'Using canonical alphabet {} and modified bases {}.'.format(
                    self.alphabet, ' '.join(
                        '{}={}'.format(*mod_b)
                        for mod_b in model_info.mod_long_names)))
        else:
            logger.info(
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

        return

    def calibrate_llr(self, score, mod_base):
        return self.calib_table.calibrate_llr(score, mod_base)


#########################
##### modVCF Writer #####
#########################

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
        if 'STRD' not in info:
            info['STRD'] = strand
        self.info_dict = info
        if sample_dict is None:
            sample_dict = OrderedDict()
        self.sample_dict = sample_dict

        if self.mod_props is not None:
            self.add_mod_props(self.mod_props)

        return

    @property
    def _sorted_format_keys(self):
        sorted_keys = sorted(self.sample_dict.keys())
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
            mod2.chrm, mod2.pos, mod2.strand)
    def __lt__(self, var2):
        return (self.chrm, self.pos, self.strand) < (
            mod2.chrm, mod2.pos, mod2.strand)
    def __le__(self, var2):
        return (self.chrm, self.pos, self.strand) <= (
            mod2.chrm, mod2.pos, mod2.strand)
    def __gt__(self, var2):
        return (self.chrm, self.pos, self.strand) > (
            mod2.chrm, mod2.pos, mod2.strand)
    def __ge__(self, var2):
        return (self.chrm, self.pos, self.strand) >= (
            mod2.chrm, mod2.pos, mod2.strand)

class ModVcfWriter(object):
    """ modVCF writer class
    """
    version_options = set(['4.2',])
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
            for ref_name, ref_len in zip(*ref_names_and_lens)]
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
        return

    def write_mod_site(self, mod_site):
        elements = [getattr(mod_site, field.lower()) for field in self.header]
        elements = ['.' if e == '' else e for e in elements]
        # VCF POS field is 1-based
        elements[self.header.index('POS')] += 1
        self.buffer.append('\t'.join(map(str, elements)))
        if len(self.buffer) > self.buffer_limit:
            self.handle.write('\n'.join(self.buffer) + '\n')
            self.buffer = []

        return

    def close(self):
        if len(self.buffer) > 0:
            self.handle.write('\n'.join(self.buffer) + '\n')
        self.handle.close()
        return

class ModBedMethylWriter(object):
    """ bedMethyl writer class

    Note that the bedMethyl format cannot store more than one modification
    type, so multiple file handles will be opened.
    """
    def __init__(self, basename, mods, mode='w', buffer_limit=OUT_BUFFER_LIMIT):
        self.basename = basename
        self.mods = mods
        self.mod_short_names, self.mod_long_names = zip(*self.mods)
        self.mode = mode
        self.buffer_limit = buffer_limit
        self.buffers = dict(
            (mod_short_name, []) for mod_short_name, _ in self.mods)
        self.handles = dict(
            (mod_short_name,
             open('{}.{}.{}'.format(self.basename, mod_long_name,
                                    mh.MOD_OUTPUT_EXTNS[mh.MOD_BEDMETHYL_NAME]),
                  self.mode, encoding='utf-8'))
            for mod_short_name, mod_long_name in self.mods)
        return

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

        return

    def close(self):
        for mod_base, handle in self.handles.items():
            if len(self.buffers[mod_base]) > 0:
                handle.write(
                    '\n'.join(self.buffers[mod_base]) + '\n')
            handle.close()
        return

class ModWigWriter(object):
    """ Modified base wiggle variableStep writer class

    Note that the wiggle/bedgraph format cannot store more than one modification
    type or multiple strands, so multiple file handles will be opened.
    """
    def __init__(
            self, basename, mods, mode='w',
            strands={'+':'fwd_strand', '-':'rev_strand'}):
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

        return

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

        return

    def close(self):
        # write all data on close since all data is required to write
        # wiggle format
        for (mod_base, strand), all_cs_mod_sites in self.mod_sites_data.items():
            with open('{}.{}.{}.{}'.format(
                    self.basename, self.mods_lookup[mod_base],
                    self.strands[strand], mh.MOD_OUTPUT_EXTNS[mh.MOD_WIG_NAME]),
                      self.mode, encoding='utf-8') as wig_fp:
                # write header
                track_name = ('Modified Base {} Proportion Modified ' +
                              '({})').format(
                                  self.mods_lookup[mod_base],
                                  self.strands[strand])
                wig_fp.write(
                    'track type=wiggle_0 name="{0}" description="{0}"\n'.format(
                        track_name))
                for chrom, cs_mod_sites in all_cs_mod_sites.items():
                    wig_fp.write('variableStep chrom={} span=1\n'.format(chrom))
                    wig_fp.write(
                        '\n'.join((
                            '{} {}'.format(pos, mod_prop)
                            for pos, mod_prop in sorted(cs_mod_sites))) + '\n')

        return


##################################
##### Mods Aggregation Class #####
##################################

class AggMods(mh.AbstractAggregationClass):
    """ Class to assist in database queries for per-site aggregation of
    modified base calls over reads.

    Warning, setting pos_index_in_memory for a large database will drastically
    increase the startup time and memory usage (~8GB for a human genome at
    CpG sites).
    """
    def __init__(self, mods_db_fn, agg_info=DEFAULT_AGG_INFO,
                 write_mod_lp=False, pos_index_in_memory=False):
        # open as read only database (default)
        self.mods_db = ModsDb(
            mods_db_fn, pos_index_in_memory=pos_index_in_memory)
        self.n_uniq_mods = None
        assert agg_info.method in AGG_METHOD_NAMES
        self.agg_method = agg_info.method
        self.binary_thresh = agg_info.binary_threshold
        self.write_mod_lp = write_mod_lp
        return

    def num_uniq(self):
        if self.n_uniq_mods is None:
            self.n_uniq_mods = self.mods_db.get_num_uniq_mod_pos()
        return self.n_uniq_mods

    def iter_uniq(self):
        # fill queue with only pos_ids for faster queue filling
        # and let workers extract pos info from an order query
        #for q_val in self.mods_db.iter_pos_id_ordered():
        # fill queue with full position information to make
        # workers avoid the ordered pos data extraction
        for q_val in self.mods_db.iter_pos():
            yield q_val
        return

    def est_binary_thresh(self, pos_scores):
        mod_cov = len(pos_scores)
        valid_cov = 0
        mod_types = set(mt for read_mods in pos_scores.values()
                        for mt in read_mods.keys())
        mods_cov = dict((mt, 0) for mt in mod_types)
        for read_pos_lps in pos_scores.values():
            mt_lps = np.array(list(read_pos_lps.values()))
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
            self, pos_scores, max_iters=5, conv_tol=0.005,
            init_thresh=0, min_prop=0.01, max_prop=0.99):
        """ [DEFUNCT]

        Estimate proportion of modified bases at a position via EM
        computation

        TODO adapt this function to handle multiple modified base types
        """
        curr_mix_prop = np.clip(np.mean(pos_scores < init_thresh),
                                min_prop, max_prop)
        for _ in range(max_iters):
            prev_mix_prop = curr_mix_prop
            if prev_mix_prop < min_prop:
                return 0.0, pos_scores.shape[0]
            if prev_mix_prop > max_prop:
                return 1.0, pos_scores.shape[0]
            curr_mix_prop = np.mean(1 / (1 + (
                np.exp(pos_scores) * (1 - prev_mix_prop) / prev_mix_prop)))
            if np.abs(curr_mix_prop - prev_mix_prop) < conv_tol:
                break

        return curr_mix_prop, pos_scores.shape[0]

    def emp_em(self, pos_scores, max_iters):
        """ Emperical EM estimation

        This method does not actually work yet. Need to compute
        emperical proportion array
        """
        # pre-computed array from ground truth
        # rows are mixture rate, cols are per-read scores
        mod_prob_tbl = np.array()

        curr_mix_prop = np.mean(pos_scores < 0)
        for _ in range(max_iters):
            curr_mix_prop = np.mean(mod_prob_tbl[curr_mix_prop, pos_scores])

        return curr_mix_prop, pos_scores.shape[0]

    def compute_mod_stats(self, mod_loc, agg_method=None, valid_read_ids=None):
        if agg_method is None:
            agg_method = self.agg_method
        if agg_method not in AGG_METHOD_NAMES:
            raise NotImplementedError(
                'No modified base proportion estimation method: {}'.format(
                    agg_method))

        pr_mod_stats = self.mods_db.get_pos_stats(
            mod_loc, return_uuids=valid_read_ids is not None)
        mod_type_stats = defaultdict(dict)
        for r_stats in pr_mod_stats:
            if (valid_read_ids is not None and
                r_stats.read_id not in valid_read_ids):
                continue
            mod_type_stats[r_stats.read_id][r_stats.mod_base] = r_stats.score
        total_cov = len(mod_type_stats)
        if total_cov == 0:
            raise mh.MegaError('No valid reads cover modified base location')
        if agg_method == BIN_THRESH_NAME:
            mod_props, valid_cov = self.est_binary_thresh(mod_type_stats)

        r0_stats = pr_mod_stats[0]
        strand = '+' if r0_stats.strand == 1 else '-'
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
                    for mod_i, mod_lp in enumerate([read_mod_scores[mod_type]
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
        return


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
