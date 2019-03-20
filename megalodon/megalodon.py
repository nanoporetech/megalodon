#!/usr/bin/env python3
import os
# set blas library environment variables (without these the cblas calls
# can completely halt processing)
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'


import re
import sys
import queue
import shutil
import argparse
import threading
import traceback
from time import sleep
import multiprocessing as mp
from collections import defaultdict

import numpy as np
from tqdm import tqdm
from tqdm._utils import _term_move_up

from megalodon import (
    decode, backends, snps, mods, mapping, megalodon_helper as mh)


DEFAULT_CONTEXT_BASES = 10
DEFAULT_EDGE_BUFFER = 5

SINGLE_LETTER_CODE = {
    'A':'A', 'C':'C', 'G':'G', 'T':'T', 'B':'[CGT]',
    'D':'[AGT]', 'H':'[ACT]', 'K':'[GT]', 'M':'[AC]',
    'N':'[ACGT]', 'R':'[AG]', 'S':'[CG]', 'V':'[ACG]',
    'W':'[AT]', 'Y':'[CT]'}

_DO_PROFILE = False
_UNEXPECTED_ERROR_CODE = 'Unexpected error'
_UNEXPECTED_ERROR_FN = 'unexpected_snp_calling_errors.{}.err'
_MAX_NUM_UNEXP_ERRORS = 50
_MAX_QUEUE_SIZE = 1000


###################################
##### Multi-processing Helper #####
###################################

def create_getter_q(getter_func, args):
    q = mp.Queue(maxsize=_MAX_QUEUE_SIZE)
    main_conn, conn = mp.Pipe()
    p = mp.Process(target=getter_func, daemon=True, args=(q, conn, *args))
    p.start()
    return q, p, main_conn


#######################################
##### Aggregate SNP and Mod Stats #####
#######################################

def _agg_snps_worker(
        locs_q, snp_stats_q, snp_prog_q, snps_db_fn, write_vcf_llr):
    agg_snps = snps.AggSnps(snps_db_fn, write_vcf_llr)

    while True:
        try:
            snp_id = locs_q.get(block=False)
        except queue.Empty:
            sleep(0.1)
            continue
        if snp_id is None:
            break

        snp_var = agg_snps.compute_snp_stats(snp_id[0])
        snp_stats_q.put(snp_var)
        snp_prog_q.put(1)

    return

def _get_snp_stats_queue(snp_stats_q, snp_conn, out_dir, do_sort=False):
    agg_snp_fn = os.path.join(out_dir, mh.OUTPUT_FNS[mh.SNP_NAME])
    if do_sort:
        all_snp_vars = []
    else:
        agg_snp_fp = snps.VcfWriter(agg_snp_fn, 'w')

    while True:
        try:
            snp_var = snp_stats_q.get(block=False)
            if do_sort:
                all_snp_vars.append(snp_var)
            else:
                agg_snp_fp.write_variant(snp_var)
        except queue.Empty:
            if snp_conn.poll():
                break
            sleep(0.1)
            continue

    while not snp_stats_q.empty():
        snp_var = snp_stats_q.get(block=False)
        if do_sort:
            all_snp_vars.append(snp_var)
        else:
            agg_snp_fp.write_variant(snp_var)

    if do_sort:
        # sort variants and write to file (requires adding __lt__, __gt__
        # methods to variant class
        with snps.VcfWriter(agg_snp_fn, 'w') as agg_snp_fp:
            for snp_var in sorted(all_snp_vars):
                agg_snp_fp.write_variant(snp_var)
    else:
        agg_snp_fp.close()

    return

def _agg_mods_worker(
        locs_q, mod_stats_q, mod_prog_q, mod_db_fn, mod_calib_fn):
    agg_mods = snps.AggMods(mods_db_fn, mods_calib_fn)

    while True:
        try:
            mod_loc = locs_q.get(block=False)
        except queue.Empty:
            sleep(0.1)
            continue
        if mod_loc is None:
            break

        mod_prop = agg_mods.compute_mod_stats(mod_loc)
        mod_stats_q.put((mod_prop, mod_loc))
        mod_prog_q.put(1)

    return

def _get_mod_stats_queue(mod_stats_q, mod_conn, out_dir, do_sort=False):
    agg_mod_fn = os.path.join(out_dir, mh.OUTPUT_FNS[mh.MOD_NAME])
    if do_sort:
        all_mods = []
    else:
        agg_mod_fp = mods.bedMethylWriter(agg_mod_fn, 'w')

    while True:
        try:
            mod_prop, mod_loc = mod_stats_q.get(block=False)
            if do_sort:
                all_mods.append((mod_loc, mod_prop))
            else:
                agg_mod_fp.write_mod(mod_prop, mod_loc)
        except queue.Empty:
            if mod_conn.poll():
                break
            sleep(0.1)
            continue

    while not mod_stats_q.empty():
        mod_prop, mod_loc = mod_stats_q.get(block=False)
        if do_sort:
            all_mods.append((mod_loc, mod_prop))
        else:
            agg_mod_fp.write_mod(mod_prop, mod_loc)
    if do_sort:
        with mods.bedMethylWriter(agg_mod_fn, 'w') as agg_mod_fp:
            for mod_loc, mod_prop in sorted(all_mods):
                agg_mod_fp.write_mod(mod_prop, mod_loc)
    else:
        agg_mod_fp.close()

    return

def _fill_locs_queue(locs_q, db_fn, agg_class, num_ps):
    agg_db = agg_class(db_fn)
    for loc in agg_db.iter_uniq():
        locs_q.put(loc)
    for _ in range(num_ps):
        locs_q.put(None)

    return

def aggregate_stats(outputs, out_dir, num_ps, write_vcf_llr):
    sys.stderr.write('Aggregating SNPs/Mods at sites over reads.\n')
    if mh.SNP_NAME in outputs and mh.MOD_NAME in outputs:
        num_ps = num_ps // 2

    agg_snps_ps, agg_mods_ps = [], []
    num_snps, num_mods = 0, 0
    if mh.SNP_NAME in outputs:
        snps_db_fn = os.path.join(out_dir, mh.OUTPUT_FNS[mh.PR_SNP_NAME][0])
        num_snps = snps.AggSnps(snps_db_fn).num_uniq()
        # create process to collect snp stats from workers
        snp_stats_q, snp_stats_p, main_snp_stats_conn = create_getter_q(
            _get_snp_stats_queue, (out_dir,))
        # create process to fill snp locs queue
        snp_filler_q = mp.Queue(maxsize=_MAX_QUEUE_SIZE)
        snp_filler_p = mp.Process(
            target=_fill_locs_queue,
            args=(snp_filler_q, snps_db_fn, snps.AggSnps, num_ps), daemon=True)
        snp_filler_p.start()
        # create worker processes to aggregate snps
        snp_prog_q = mp.Queue(maxsize=_MAX_QUEUE_SIZE)
        for _ in range(num_ps):
            p = mp.Process(
                target=_agg_snps_worker,
                args=(snp_filler_q, snp_stats_q, snp_prog_q, snps_db_fn,
                      write_vcf_llr), daemon=True)
            p.start()
            agg_snps_ps.append(p)

    if mh.MOD_NAME in outputs:
        mods_db_fn = os.path.join(out_dir, mh.OUTPUT_FNS[mh.PR_MOD_NAME][0])
        num_mods = mods.AggMods(mods_db_fn).num_uniq()
        # create process to collect mods stats from workers
        mod_stats_q, mod_stats_p, main_mod_stats_conn = create_getter_q(
            _get_mod_stats_queue, (out_dir,))
        # create process to fill mod locs queue
        mod_filler_q = mp.Queue(maxsize=_MAX_QUEUE_SIZE)
        mod_filler_p = mp.Process(
            target=_fill_locs_queue,
            args=(mod_filler_q, mods_db_fn, mods.AggMods, num_ps), daemon=True)
        mod_filler_p.start()
        # create worker processes to aggregate mods
        mod_prog_q = mp.Queue(maxsize=_MAX_QUEUE_SIZE)
        for _ in range(num_ps):
            p = mp.Process(
                target=_agg_mods_worker,
                args=(mod_filler_q, mod_stats_q, mod_prog_q, mods_db_fn,
                      mods_calib_fn), daemon=True)
            p.start()
            agg_mods_ps.append(p)

    # TODO create progress process

    # join filler processes first
    if mh.SNP_NAME in outputs:
        snp_filler_p.join()
        for agg_snps_p in agg_snps_ps:
            agg_snps_p.join()
        # send to conn
        if snp_stats_p.is_alive():
            main_snp_stats_conn.send(True)
        snp_stats_p.join()
    if mh.MOD_NAME in outputs:
        for agg_mods_p in agg_mods_ps:
            agg_mods_p.join()

    return


###########################
##### Read Processing #####
###########################

def process_read(
        raw_sig, read_id, model_info, bc_q, caller_conn, snps_to_test,
        snp_all_paths, snp_calib_tbl, snps_q, mods_q, alphabet_info,
        context_bases=DEFAULT_CONTEXT_BASES, edge_buffer=DEFAULT_EDGE_BUFFER):
    if model_info.is_cat_mod:
        bc_weights, mod_weights = model_info.run_model(
            raw_sig, alphabet_info.n_can_state)
    else:
        bc_weights = model_info.run_model(raw_sig)

    r_post = decode.crf_flipflop_trans_post(bc_weights, log=True)
    r_seq, score, runlen = decode.decode_post(
        r_post, alphabet_info.collapse_alphabet)
    if bc_q is not None:
        bc_q.put((read_id, r_seq))

    # if no mapping connection return after basecalls are passed out
    if caller_conn is None: return

    # map read and record mapping from reference to query positions
    r_ref_seq, r_to_q_poss, r_ref_pos = mapping.map_read(
        r_seq, read_id, caller_conn)
    np_ref_seq = np.array([
        mh.ALPHABET.find(b) for b in r_ref_seq], dtype=np.uintp)

    # get mapped start in post and run len to mapped bit of output
    post_mapped_start = sum(runlen[:r_ref_pos.q_trim_start])
    rl_cumsum = np.cumsum(np.concatenate([
        [0], runlen[r_ref_pos.q_trim_start:r_ref_pos.q_trim_end]]))

    if snps_q is not None:
        try:
            snps_q.put((
                snps.call_read_snps(
                    r_ref_pos, snps_to_test, edge_buffer, context_bases,
                    np_ref_seq, rl_cumsum, r_to_q_poss, r_post,
                    post_mapped_start, snp_all_paths, snp_calib_tbl),
                (read_id, r_ref_pos.chrm, r_ref_pos.strand)))
        except mh.MegaError:
            pass
    if mods_q is not None:
        r_post_w_mods = np.concatenate([r_post, mod_weights], axis=1)
        try:
            mods_q.put((
                mods.call_read_mods(
                    r_ref_pos, edge_buffer, context_bases, r_ref_seq,
                    np_ref_seq, rl_cumsum, r_to_q_poss, r_post_w_mods,
                    post_mapped_start, alphabet_info),
                (read_id, r_ref_pos.chrm, r_ref_pos.strand)))
        except mh.MegaError:
            pass

    return


############################
##### Multi-processing #####
############################

def _get_bc_queue(bc_q, bc_conn, out_dir, bc_fmt):
    bc_fp = open(os.path.join(
        out_dir, mh.OUTPUT_FNS[mh.BC_NAME] + '.' + bc_fmt), 'w')

    while True:
        try:
            # TODO add quality output to add fastq option
            read_id, r_seq = bc_q.get(block=False)
            bc_fp.write('>{}\n{}\n'.format(read_id, r_seq))
            bc_fp.flush()
        except queue.Empty:
            if bc_conn.poll():
                break
            sleep(0.1)
            continue

    while not bc_q.empty():
        read_id, r_seq = bc_q.get(block=False)
        bc_fp.write('>{}\n{}\n'.format(read_id, r_seq))
        bc_fp.flush()
    bc_fp.close()

    return

def get_read_files(fast5s_dir):
    """Get all fast5 files recursively below this directory
    """
    fast5_fns = []
    # walk through directory structure searching for fast5 files
    for root, _, fns in os.walk(fast5s_dir, followlinks=True):
        for fn in fns:
            if not fn.endswith('.fast5'): continue
            fast5_fns.append(os.path.join(root, fn))

    return fast5_fns

def _fill_files_queue(fast5_q, fast5_fns, num_ps):
    for fast5_fn in fast5_fns:
        fast5_q.put(fast5_fn)
    for _ in range(num_ps):
        fast5_q.put(None)

    return

def _process_reads_worker(
        fast5_q, bc_q, mo_q, snps_q, failed_reads_q, mods_q, caller_conn,
        model_info, snps_to_test, snp_all_paths, snp_calib_tbl, alphabet_info):
    model_info.prep_model_worker()

    while True:
        try:
            fast5_fn = fast5_q.get(block=False)
        except queue.Empty:
            sleep(0.1)
            continue

        if fast5_fn is None:
            if caller_conn is not None:
                caller_conn.send(True)
            break

        try:
            raw_sig, read_id = mh.extract_read_data(fast5_fn)
            process_read(
                raw_sig, read_id, model_info, bc_q, caller_conn, snps_to_test,
                snp_all_paths, snp_calib_tbl, snps_q, mods_q, alphabet_info)
            failed_reads_q.put((False, None, None, None))
        except KeyboardInterrupt:
            failed_reads_q.put((True, 'Keyboard interrupt', fast5_fn, None))
            return
        except mh.MegaError as e:
            failed_reads_q.put((True, str(e), fast5_fn, None))
        except:
            failed_reads_q.put((
                True, _UNEXPECTED_ERROR_CODE, fast5_fn, traceback.format_exc()))

    return

if _DO_PROFILE:
    _process_reads_wrapper = _process_reads_worker
    def _process_reads_worker(*args):
        import cProfile
        cProfile.runctx('_process_reads_wrapper(*args)', globals(), locals(),
                        filename='read_processing.prof')
        return


##################################
##### Dynamic error updating #####
##################################

def format_fail_summ(header, fail_summ=[], num_reads=0, num_errs=None):
    summ_errs = sorted(fail_summ)[::-1]
    if num_errs is not None:
        summ_errs = summ_errs[:num_errs]
        if len(summ_errs) < num_errs:
            summ_errs.extend([(None, '') for _ in range(
                num_errs - len(summ_errs))])
    errs_str = '\n'.join("{:8.1f}% ({:>7} reads)".format(
        100 * n_fns / float(num_reads), n_fns) + " : " + '{:<80}'.format(err)
        if (n_fns is not None and num_reads > 0) else
        '     -----' for n_fns, err in summ_errs)
    return '\n'.join((header, errs_str))

def prep_errors_bar(num_update_errors, tot_reads, suppress_progress):
    if num_update_errors > 0 and not suppress_progress:
        # add lines for dynamic error messages
        sys.stderr.write(
            '\n'.join(['' for _ in range(num_update_errors + 2)]))
    bar, prog_prefix, bar_header = None, None, None
    if suppress_progress:
        num_update_errors = 0
    else:
        bar = tqdm(total=tot_reads, smoothing=0)
    if num_update_errors > 0:
        prog_prefix = ''.join(
            [_term_move_up(),] * (num_update_errors + 1)) + '\r'
        bar_header = (
            str(num_update_errors) + ' most common unsuccessful read types:')
        # write failed read update header
        bar.write(prog_prefix + format_fail_summ(
            bar_header, num_errs=num_update_errors), file=sys.stderr)

    return bar, prog_prefix, bar_header

def _get_fail_queue(
        failed_reads_q, f_conn, tot_reads, num_update_errors,
        suppress_progress):
    reads_called = 0
    failed_reads = defaultdict(list)
    bar, prog_prefix, bar_header = prep_errors_bar(
        num_update_errors, tot_reads, suppress_progress)
    while True:
        try:
            is_err, err_type, fast5_fn, err_tb = failed_reads_q.get(block=False)
            if is_err:
                failed_reads[err_type].append(fast5_fn)
                if err_type == _UNEXPECTED_ERROR_CODE:
                    if len(failed_reads[err_type]) == 1:
                        unexp_err_fp = open(_UNEXPECTED_ERROR_FN.format(
                            np.random.randint(10000)), 'w')
                    if len(failed_reads[err_type]) >= _MAX_NUM_UNEXP_ERRORS:
                        unexp_err_fp.close()
                    else:
                        unexp_err_fp.write(
                            fast5_fn + '\n:::\n' + err_tb + '\n\n\n')
                        unexp_err_fp.flush()

            if not suppress_progress: bar.update(1)
            reads_called += 1
            if num_update_errors > 0:
                bar.write(prog_prefix + format_fail_summ(
                    bar_header,
                    [(len(fns), err) for err, fns in failed_reads.items()],
                    reads_called, num_update_errors), file=sys.stderr)
        except queue.Empty:
            # check if all reads are done signal was sent from main thread
            if f_conn.poll():
                break
            try:
                sleep(0.1)
            except KeyboardInterrupt:
                # exit gracefully on keyboard inturrupt
                return
            continue

    if not suppress_progress: bar.close()
    if len(failed_reads[_UNEXPECTED_ERROR_CODE]) >= 1:
        sys.stderr.write((
            '******* WARNING *******\n\tUnexpected errors occured. See full ' +
            'error stack traces for first (up to) {0:d} errors in ' +
            '"{1}"\n').format(_MAX_NUM_UNEXP_ERRORS, unexp_err_fp.name))

    return


###############################
##### All read processing #####
###############################

def process_all_reads(
        fast5s_dir, num_reads, model_info, outputs, out_dir, bc_fmt, aligner,
        add_chr_ref, snps_data, num_ps, num_update_errors, suppress_progress,
        alphabet_info, db_safety):
    sys.stderr.write('Searching for reads.\n')
    fast5_fns = get_read_files(fast5s_dir)
    if num_reads is not None:
        fast5_fns = fast5_fns[:num_reads]

    sys.stderr.write('Preparing workers and calling reads.\n')
    # read filename queue filler
    fast5_q = mp.Queue(maxsize=_MAX_QUEUE_SIZE)
    files_p = mp.Process(
        target=_fill_files_queue, args=(fast5_q, fast5_fns, num_ps),
        daemon=True)
    files_p.start()
    # progress and failed reads getter
    failed_reads_q, f_p, main_f_conn = create_getter_q(
            _get_fail_queue, (len(fast5_fns), num_update_errors,
                              suppress_progress))

    # start output type getters/writers
    (bc_q, bc_p, main_bc_conn, mo_q, mo_p, main_mo_conn, snps_q, snps_p,
     main_snps_conn, mods_q, mods_p, main_mods_conn) = [None,] * 12
    if mh.BC_NAME in outputs:
        bc_q, bc_p, main_bc_conn = create_getter_q(
            _get_bc_queue, (out_dir, bc_fmt))
    if mh.MAP_NAME in outputs:
        mo_q, mo_p, main_mo_conn = create_getter_q(
            mapping._get_map_queue, (out_dir, aligner.ref_names_and_lens,
                             aligner.out_fmt, aligner.ref_fn))
    if mh.PR_SNP_NAME in outputs:
        snps_db_fn, snps_txt_fn = mh.OUTPUT_FNS[mh.PR_SNP_NAME]
        snps_txt_fn = (os.path.join(out_dir, snps_txt_fn)
                       if snps_data.write_snps_txt else None)
        snps_q, snps_p, main_snps_conn = create_getter_q(
            snps._get_snps_queue, (
                snps_data.snp_id_tbl, os.path.join(out_dir, snps_db_fn),
                snps_txt_fn, db_safety))
    if mh.PR_MOD_NAME in outputs:
        mods_db_fn, mods_txt_fn = mh.OUTPUT_FNS[mh.PR_MOD_NAME]
        mods_txt_fn = (os.path.join(out_dir, mods_db_fn)
                       if alphabet_info.write_mods_txt else None)
        mods_q, mods_p, main_mods_conn = create_getter_q(
            mods._get_mods_queue, (
                os.path.join(out_dir, mods_db_fn), mods_txt_fn, db_safety))

    proc_reads_ps, map_conns = [], []
    for _ in range(num_ps):
        if aligner is None:
            map_conn, caller_conn = None, None
        else:
            map_conn, caller_conn = mp.Pipe()
        map_conns.append(map_conn)
        p = mp.Process(
            target=_process_reads_worker, args=(
                fast5_q, bc_q, mo_q, snps_q, failed_reads_q, mods_q,
                caller_conn, model_info, snps_data.snps_to_test,
                snps_data.all_paths, snps_data.calib_table, alphabet_info))
        p.daemon = True
        p.start()
        proc_reads_ps.append(p)
    sleep(0.1)

    # perform mapping in threads for mappy shared memory interface
    if aligner is None:
        map_read_ts = None
    else:
        map_read_ts = []
        for map_conn in map_conns:
            t = threading.Thread(
                target=mapping._map_read_worker,
                args=(aligner, map_conn, mo_q, add_chr_ref))
            t.daemon = True
            t.start()
            map_read_ts.append(t)

    try:
        files_p.join()
    except KeyboardInterrupt:
        sys.stderr.write('Exiting due to keyboard interrupt.\n')
        sys.exit(1)
    for proc_reads_p in proc_reads_ps:
        proc_reads_p.join()
    if map_read_ts is not None:
        for map_t in map_read_ts:
            map_t.join()
    # comm to getter processes to return
    if f_p.is_alive():
        main_f_conn.send(True)
        f_p.join()
    for on, p, main_conn in (
            (mh.BC_NAME, bc_p, main_bc_conn),
            (mh.MAP_NAME, mo_p, main_mo_conn),
            (mh.PR_SNP_NAME, snps_p, main_snps_conn),
            (mh.PR_MOD_NAME, mods_p, main_mods_conn)):
        if on in outputs and p.is_alive():
            main_conn.send(True)
            if on in (mh.PR_SNP_NAME, mh.PR_MOD_NAME):
                sys.stderr.write(
                    'Waiting for snps and/or mods process to complete ' +
                    'database indexing.\n')
            p.join()

    return


################################
########## Parse SNPs ##########
################################

class AlphabetInfo(object):
    def _parse_cat_mods(self):
        self.n_can_state = (self.ncan_base + self.ncan_base) * (
            self.ncan_base + 1)
        self.nmod_base = len(set(self.alphabet)) - self.ncan_base
        # record the canonical base group indices for modified base
        # categorical output layer
        self.can_mods_offsets = np.cumsum([0] + [
            self.collapse_alphabet.count(can_b)
            for can_b in self.alphabet[:self.ncan_base]]).astype(np.uintp)

        # create table of string mod labels to ints taken by
        # categorical modifications model
        self.str_to_int_mod_labels = {}
        can_grouped_mods = defaultdict(int)
        for mod_lab, can_lab in zip(
                self.alphabet[self.ncan_base:],
                self.collapse_alphabet[self.ncan_base:]):
            can_grouped_mods[can_lab] += 1
            self.str_to_int_mod_labels[mod_lab] = can_grouped_mods[can_lab]

        return

    def _parse_mod_motifs(self, all_mod_motifs_raw):
        # note only works for mod_refactor models currently
        self.all_mod_motifs = []
        if all_mod_motifs_raw is None:
            for can_base, mod_base in zip(
                    self.collapse_alphabet[self.ncan_base:],
                    self.alphabet[self.ncan_base:]):
                self.all_mod_motifs.append((
                    re.compile(can_base), 0, mod_base, can_base))
        else:
            # parse detection motifs
            for mod_motifs_raw in all_mod_motifs_raw:
                mod_base, motif_raw = mod_motifs_raw.split(':')
                assert mod_base in self.alphabet[self.ncan_base:], (
                    'Modified base label ({}) not found in model ' +
                    'alphabet ({}).').format(mod_base, self.alphabet)
                raw_motif, pos = motif_raw.split('-')
                pos = int(pos)
                mod_base_idx = self.alphabet.find(mod_base)
                assert mod_base_idx != -1, (
                    'Invalid modified base motif. Mod base ({}) not found ' +
                    'in alphabet ({}).').format(mod_base, self.alphabet)
                assert (self.collapse_alphabet[mod_base_idx] ==
                        raw_motif[pos]), (
                            'Invalid modified base motif. Raw motif modified ' +
                            'position ({}) base ({}) does not match ' +
                            'collapsed alphabet value ({}).').format(
                                pos, raw_motif[pos],
                                self.collapse_alphabet[mod_base_idx])
                motif = re.compile(''.join(
                    SINGLE_LETTER_CODE[letter] for letter in raw_motif))
                self.all_mod_motifs.append((motif, pos, mod_base, raw_motif))

        return

    def __init__(
            self, model_info, all_mod_motifs_raw, mod_all_paths,
            override_alphabet, write_mods_txt):
        self.mod_all_paths = mod_all_paths
        self.write_mods_txt = write_mods_txt

        self.alphabet = model_info.alphabet
        self.collapse_alphabet = model_info.collapse_alphabet
        self.ncan_base = len(set(self.collapse_alphabet))
        if override_alphabet is not None:
            self.alphabet, self.collapse_alphabet = override_alphabet
        try:
            self.alphabet = self.alphabet.decode()
            self.collapse_alphabet = self.collapse_alphabet.decode()
        except:
            pass
        sys.stderr.write('Using alphabet {} collapsed to {}.\n'.format(
            self.alphabet, self.collapse_alphabet))

        self.nbase = len(self.alphabet)
        if model_info.is_cat_mod:
            self._parse_cat_mods()
            assert (
                model_info.output_size - self.n_can_state ==
                self.nmod_base + 1), (
                    'Alphabet ({}) and model number of modified bases ({}) ' +
                    'do not agree. See --override-alphabet.'
                ).format(self.alphabet,
                         model_info.output_size - self.n_can_state - 1)

        # parse mod motifs or use "swap" base if no motif provided
        self._parse_mod_motifs(all_mod_motifs_raw)

        return

def mkdir(out_dir, overwrite):
    if os.path.exists(out_dir):
        if not overwrite:
            sys.stderr.write(
                '*' * 100 + '\nERROR: --output-directory exists and ' +
                '--overwrite is not set.\n' + '*' * 100 + '\n')
            sys.exit(1)
        if os.path.isfile(out_dir) or os.path.islink(out_dir):
            os.remove(out_dir)
        else:
            shutil.rmtree(out_dir)
    os.mkdir(out_dir)

    return


##########################
########## Main ##########
##########################

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'fast5s_dir',
        help='Directory containing raw fast5 (will be searched recursively).')

    mdl_grp = parser.add_argument_group('Model Arguments')
    mdl_grp.add_argument(
        '--taiyaki-model-filename',
        help='Taiyaki model checkpoint file.')
    mdl_grp.add_argument(
        '--flappie-model-name',
        help='Flappie model name.')

    out_grp = parser.add_argument_group('Output Arguments')
    out_grp.add_argument(
        '--outputs', nargs='+',
        default=['basecalls',], choices=tuple(mh.OUTPUT_FNS.keys()),
        help='Output type(s) to produce. Default: %(default)s')
    out_grp.add_argument(
        '--output-directory',
        default='megalodon_results',
        help='Directory to store output results. Default: %(default)s')
    out_grp.add_argument(
        '--overwrite', action='store_true',
        help='Overwrite output directory if it exists.')
    out_grp.add_argument(
        '--num-reads', type=int,
        help='Number of reads to process. Default: All reads')

    bc_grp = parser.add_argument_group('Basecall Arguments')
    bc_grp.add_argument(
        '--basecalls-format', choices=mh.BC_OUT_FMTS, default=mh.BC_OUT_FMTS[0],
        help='Basecalls output format. Choices: {}'.format(
            ', '.join(mh.BC_OUT_FMTS)))

    map_grp = parser.add_argument_group('Mapping Arguments')
    map_grp.add_argument(
        '--reference',
        help='Reference FASTA file used for mapping called reads.')
    map_grp.add_argument(
        '--mappings-format', choices=mh.MAP_OUT_FMTS,
        default=mh.MAP_OUT_FMTS[0],
        help='Mappings output format. Choices: {}'.format(
            ', '.join(mh.MAP_OUT_FMTS)))

    snp_grp = parser.add_argument_group('SNP Arguments')
    snp_grp.add_argument(
        '--snp-filename',
        help='SNPs to call for each read in VCF format (required for output).')
    snp_grp.add_argument(
        '--snp-calibration-filename',
        help='File containing emperical calibration for SNP scores. As ' +
        'created by megalodon/scripts/calibrate_snp_scores.py.')
    snp_grp.add_argument(
        '--prepend-chr-vcf', action='store_true',
        help='Prepend "chr" to chromosome names from VCF to match ' +
        'reference names.')
    snp_grp.add_argument(
        '--max-snp-size', type=int, default=5,
        help='Maximum difference in number of reference and alternate bases. ' +
        'Default: %(default)d')
    snp_grp.add_argument(
        '--snp-all-paths', action='store_true',
        help='Compute forwards algorithm all paths score. (Default: Viterbi ' +
        'best-path score)')
    snp_grp.add_argument(
        '--write-snps-text', action='store_true',
        help='Write per-read SNP calls out to a text file. Default: ' +
        'Only ouput to database.')
    snp_grp.add_argument(
        '--write-vcf-llr', action='store_true',
        help='Write log-likelihood ratios out in non-standard VCF field.')

    mod_grp = parser.add_argument_group('Modified Base Arguments')
    mod_grp.add_argument(
        '--mod-motifs', nargs='+',
        help='Restrict modified base calls to specified motifs. Format as ' +
        '"[mod_base]:[motif]-[relative_pos]". For CpG, dcm and dam calling ' +
        'use "Z:CG-0 Z:CCWGG-1 Y:GATC-1".')
    mod_grp.add_argument(
        '--mod-all-paths', action='store_true',
        help='Compute forwards algorithm all paths score for modified base ' +
        'calls. (Default: Viterbi best-path score)')
    mod_grp.add_argument(
        '--write-mods-text', action='store_true',
        help='Write per-read modified bases out to a text file. Default: ' +
        'Only ouput to database.')

    misc_grp = parser.add_argument_group('Miscellaneous Arguments')
    misc_grp.add_argument(
        '--override-alphabet', nargs=2,
        help='Override alphabet and collapse alphabet from model.')
    misc_grp.add_argument(
        '--prepend-chr-ref', action='store_true',
        help='Prepend "chr" to chromosome names from reference to match ' +
        'VCF names.')
    misc_grp.add_argument(
        '--suppress-progress', action='store_true',
        help='Suppress progress bar output.')
    misc_grp.add_argument(
        '--verbose-read-progress', type=int, default=0,
        help='Output verbose output on read progress. Outputs N most ' +
        'common points where reads could not be processed further. ' +
        'Default: %(default)d')
    misc_grp.add_argument(
        '--processes', type=int, default=1,
        help='Number of parallel processes. Default: %(default)d')
    misc_grp.add_argument(
        '--device', default=None, type=int,
        help='CUDA device to use (only valid for taiyaki), or None to use CPU')
    misc_grp.add_argument(
        '--database-safety', type=int, default=1,
        help='Setting for database performance versus corruption protection. ' +
        'Options: 0 (DB corruption on application crash), 1 (DB corruption ' +
        'on system crash), 2 (DB safe mode). Default: %(default)d')

    return parser


def _main():
    args = get_parser().parse_args()
    if _DO_PROFILE:
        if args.processes > 1:
            msg = ('Running profiling with multiple processes is ' +
                   'not allowed. Setting to single process.')
            args.processes = 1
        else:
            msg = 'Running profiling. This may slow processing.'
        sys.stderr.write(
            '*' * 100 + '\nWARNING: ' + msg + '\n' + '*' * 100 + '\n')

    mkdir(args.output_directory, args.overwrite)
    model_info = backends.ModelInfo(
        args.flappie_model_name, args.taiyaki_model_filename, args.device)

    alphabet_info = AlphabetInfo(
        model_info, args.mod_motifs, args.mod_all_paths, args.override_alphabet,
        args.write_mods_text)
    if mh.PR_MOD_NAME in args.outputs and not mh.MOD_NAME in args.outputs:
        args.output.append(mh.MOD_NAME)
    if model_info.is_cat_mod and mh.PR_MOD_NAME not in args.outputs:
        sys.stderr.write(
            '*' * 100 + '\nWARNING: Categorical modifications model ' +
            'provided, but {} not requested '.format(mh.PR_MOD_NAME) +
            '(via --outputs). Modified base output will not be produced.\n' +
            '*' * 100 + '\n')
    if args.mod_motifs is not None and mh.PR_MOD_NAME not in args.outputs:
        sys.stderr.write(
            '*' * 100 + '\nWARNING: --mod-motifs provided, but ' +
            '{} not requested '.format(mh.PR_MOD_NAME) +
            '(via --outputs). Argument will be ignored.\n' + '*' * 100 + '\n')

    if mh.SNP_NAME in args.outputs and not mh.PR_SNP_NAME in args.outputs:
        args.output.append(mh.PR_SNP_NAME)
    if mh.PR_SNP_NAME in args.outputs and args.snp_filename is None:
        sys.stderr.write(
            '*' * 100 + '\nERROR: {} output requested, '.format(
                mh.PR_SNP_NAME) +
            'but --snp-filename provided.\n' + '*' * 100 + '\n')
        sys.exit(1)
    if mh.PR_SNP_NAME in args.outputs and not (
            model_info.is_cat_mod or
            mh.nstate_to_nbase(model_info.output_size) == 4):
        sys.stderr.write(
            '*' * 100 + '\nERROR: SNP calling from standard modified base ' +
            'flip-flop model is not supported.\n' + '*' * 100 + '\n')
        sys.exit(1)
    # snps data object loads with None snp_fn for easier handling downstream
    snps_data = snps.SnpData(
        args.snp_filename, args.prepend_chr_vcf, args.max_snp_size,
        args.snp_all_paths, args.write_snps_text, args.snp_calibration_filename)
    if args.snp_filename is not None and mh.PR_SNP_NAME not in args.outputs:
        sys.stderr.write(
            '*' * 100 + '\nWARNING: --snps-filename provided, but ' +
            'per_read_snp not requested (via --outputs). Argument will be ' +
            'ignored.\n' + '*' * 100 + '\n')

    do_align = len(mh.ALIGN_OUTPUTS.intersection(args.outputs)) > 0
    if do_align:
        if args.reference is None:
            sys.stderr.write(
                '*' * 100 + '\nERROR: Output(s) requiring reference ' +
                'alignment requested, but --reference not provided.\n' +
                '*' * 100 + '\n')
            sys.exit(1)
        sys.stderr.write('Loading reference.\n')
        aligner = mapping.alignerPlus(
            str(args.reference), preset=str('map-ont'), best_n=1)
        setattr(aligner, 'out_fmt', args.mappings_format)
        setattr(aligner, 'ref_fn', args.reference)
        if mh.MAP_NAME in args.outputs:
            aligner.add_ref_names(args.reference)
    else:
        aligner = None
    if args.reference is not None and not do_align:
        sys.stderr.write(
            '*' * 100 + '\nWARNING: --reference provided, but no --outputs ' +
            'requiring alignment was requested. Argument will be ignored.\n' +
            '*' * 100 + '\n')

    process_all_reads(
        args.fast5s_dir, args.num_reads, model_info, args.outputs,
        args.output_directory, args.basecalls_format, aligner,
        args.prepend_chr_ref, snps_data, args.processes,
        args.verbose_read_progress, args.suppress_progress,
        alphabet_info, args.database_safety)

    if mh.SNP_NAME in args.outputs or mh.MOD_NAME in args.outputs:
        aggregate_stats(args.outputs, args.output_directory, args.processes,
                        args.write_vcf_llr)

    return

if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
