#!/usr/bin/env python3
import os
# set blas library environment variables (without these the cblas calls
# can completely halt processing)
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import re
import sys
import h5py
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
    aggregate, backends, decode, fast5_io, logging, mapping, mods, snps,
    megalodon_helper as mh)
from megalodon._version import MEGALODON_VERSION


_DO_PROFILE = False
_UNEXPECTED_ERROR_CODE = 'Unexpected error'
_UNEXPECTED_ERROR_FN = 'unexpected_snp_calling_errors.{}.err'
_MAX_NUM_UNEXP_ERRORS = 50


###########################
##### Read Processing #####
###########################

def handle_errors(func, args, r_vals, out_q, fast5_fn, failed_reads_q):
    try:
        out_q.put((func(*args), r_vals))
    except KeyboardInterrupt:
        failed_reads_q.put(
            (True, False, 'Keyboard interrupt', fast5_fn, None))
        return
    except mh.MegaError as e:
        failed_reads_q.put((True, False, str(e), fast5_fn, None))
    except:
        failed_reads_q.put((
            True, False, _UNEXPECTED_ERROR_CODE, fast5_fn,
            traceback.format_exc()))
    return

def process_read(
        raw_sig, read_id, model_info, bc_q, caller_conn, snps_to_test,
        snp_all_paths, snp_calib_tbl, snp_context_bases, snps_q, mods_q,
        mods_info, fast5_fn, failed_reads_q, edge_buffer):
    """ Workhorse per-read megalodon function (connects all the parts)
    """
    if model_info.is_cat_mod:
        bc_weights, mod_weights = model_info.run_model(
            raw_sig, mods_info.n_can_state)
        can_nmods = model_info.can_nmods
    else:
        mod_weights, can_nmods = None, None
        bc_weights = model_info.run_model(raw_sig)

    r_post = decode.crf_flipflop_trans_post(bc_weights, log=True)
    if mods_q is not None:
        r_post_w_mods = np.concatenate([r_post, mod_weights], axis=1)
    if not mods_info.do_output_mods:
        mod_weights = None
    r_seq, score, rl_cumsum, mods_scores = decode.decode_post(
        r_post, mods_info.alphabet, mod_weights, can_nmods)
    if bc_q is not None:
        bc_q.put((read_id, r_seq, mods_scores))

    # if no mapping connection return after basecalls are passed out
    if caller_conn is None: return

    # map read and record mapping from reference to query positions
    r_ref_seq, r_to_q_poss, r_ref_pos, r_cigar = mapping.map_read(
        r_seq, read_id, caller_conn)
    np_ref_seq = np.array([
        mh.ALPHABET.find(b) for b in r_ref_seq], dtype=np.uintp)

    # get mapped start in post and run len to mapped bit of output
    post_mapped_start = rl_cumsum[r_ref_pos.q_trim_start]
    mapped_rl_cumsum = rl_cumsum[
        r_ref_pos.q_trim_start:r_ref_pos.q_trim_end + 1] - post_mapped_start

    if snps_q is not None:
        handle_errors(
            func=snps.call_read_snps,
            args=(r_ref_pos, snps_to_test, edge_buffer, snp_context_bases,
                  np_ref_seq, mapped_rl_cumsum, r_to_q_poss, r_post,
                  post_mapped_start, snp_all_paths, snp_calib_tbl),
            r_vals=(read_id, r_ref_pos.chrm, r_ref_pos.strand,
                    r_ref_pos.start, r_ref_seq, len(r_seq),
                    r_ref_pos.q_trim_start, r_ref_pos.q_trim_end, r_cigar),
            out_q=snps_q, fast5_fn=fast5_fn, failed_reads_q=failed_reads_q)
    if mods_q is not None:
        handle_errors(
            func=mods.call_read_mods,
            args=(r_ref_pos, edge_buffer, r_ref_seq, np_ref_seq,
                  mapped_rl_cumsum, r_to_q_poss, r_post_w_mods,
                  post_mapped_start, mods_info),
            r_vals=(read_id, r_ref_pos.chrm, r_ref_pos.strand,
                    r_ref_pos.start, r_ref_seq, len(r_seq),
                    r_ref_pos.q_trim_start, r_ref_pos.q_trim_end, r_cigar),
            out_q=mods_q, fast5_fn=fast5_fn, failed_reads_q=failed_reads_q)

    return


############################
##### Multi-processing #####
############################

def _get_bc_queue(
        bc_q, bc_conn, out_dir, bc_fmt, do_output_mods, mod_long_names):
    bc_fp = open(os.path.join(
        out_dir, mh.OUTPUT_FNS[mh.BC_NAME] + '.' + bc_fmt), 'w')
    if do_output_mods:
        mods_fp = h5py.File(os.path.join(
            out_dir, mh.OUTPUT_FNS[mh.BC_MODS_NAME]))
        mods_fp.create_group('Reads')
        mods_fp.create_dataset(
            'mod_long_names', data=np.array(mod_long_names, dtype='S'),
            dtype=h5py.special_dtype(vlen=str))

    while True:
        try:
            # TODO add quality output to add fastq option
            read_id, r_seq, mods_scores = bc_q.get(block=False)
            bc_fp.write('>{}\n{}\n'.format(read_id, r_seq))
            bc_fp.flush()
            if do_output_mods:
                try:
                    mods_fp.create_dataset(
                        'Reads/' + read_id, data=mods_scores,
                        compression="gzip")
                except RuntimeError:
                    # same read_id encountered previously
                    pass
        except queue.Empty:
            if bc_conn.poll():
                break
            sleep(0.1)
            continue

    while not bc_q.empty():
        read_id, r_seq, mods_scores = bc_q.get(block=False)
        bc_fp.write('>{}\n{}\n'.format(read_id, r_seq))
        bc_fp.flush()
        if do_output_mods:
            try:
                mods_fp.create_dataset(
                    'Reads/' + read_id, data=mods_scores,
                    compression="gzip")
            except RuntimeError:
                # same read_id encountered previously
                pass

    bc_fp.close()
    if do_output_mods:
        mods_fp.close()

    return

def _process_reads_worker(
        read_file_q, bc_q, snps_q, failed_reads_q, mods_q, caller_conn,
        model_info, snps_to_test, snp_all_paths, snp_calib_tbl,
        snp_context_bases, mods_info, edge_buffer, device):
    model_info.prep_model_worker(device)

    while True:
        try:
            try:
                fast5_fn, read_id = read_file_q.get(block=False)
            except queue.Empty:
                sleep(0.1)
                continue

            if fast5_fn is None:
                if caller_conn is not None:
                    caller_conn.send(True)
                break

            raw_sig = fast5_io.get_signal(fast5_fn, read_id, scale=True)
            process_read(
                raw_sig, read_id, model_info, bc_q, caller_conn, snps_to_test,
                snp_all_paths, snp_calib_tbl, snp_context_bases, snps_q, mods_q,
                mods_info, fast5_fn, failed_reads_q, edge_buffer)
            failed_reads_q.put((False, True, None, None, None))
        except KeyboardInterrupt:
            failed_reads_q.put((
                True, True, 'Keyboard interrupt', fast5_fn, None))
            return
        except mh.MegaError as e:
            failed_reads_q.put((True, True, str(e), fast5_fn, None))
        except:
            failed_reads_q.put((
                True, True, _UNEXPECTED_ERROR_CODE, fast5_fn,
                traceback.format_exc()))

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

def _fill_files_queue(
        read_file_q, fast5s_dir, num_reads, recursive, num_ps, num_reads_conn):
    logger = logging.get_logger()
    read_ids = set()
    # fill queue with read filename and read id tuples
    for fast5_fn, read_id in fast5_io.iterate_fast5_reads(
            fast5s_dir, num_reads, recursive):
        if read_id in read_ids:
            logger.debug(
                ('Read ID ({}) found in previous read and will not ' +
                 'process from {}.').format(read_id, fast5_fn))
            continue
        read_file_q.put((fast5_fn, read_id))
        read_ids.add(read_id)
    # add None to indicate that read processes should return
    for _ in range(num_ps):
        read_file_q.put((None, None))
    num_reads_conn.send(len(read_ids))

    return

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

def prep_errors_bar(
        num_update_errors, tot_reads, suppress_progress, curr_num_reads=0,
        start_time=None):
    if num_update_errors > 0 and not suppress_progress:
        # add lines for dynamic error messages
        sys.stderr.write(
            '\n'.join(['' for _ in range(num_update_errors + 2)]))
    bar, prog_prefix, bar_header = None, None, None
    if suppress_progress:
        num_update_errors = 0
    else:
        bar = tqdm(total=tot_reads, smoothing=0, initial=curr_num_reads)
        if start_time is not None:
            bar.start_t = start_time
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
        failed_reads_q, f_conn, getter_num_reads_conn, num_update_errors,
        suppress_progress):
    def update_prog(reads_called, unexp_err_fp):
        if is_err:
            failed_reads[err_type].append(fast5_fn)
            if err_type == _UNEXPECTED_ERROR_CODE:
                if len(failed_reads[_UNEXPECTED_ERROR_CODE]) == 1:
                    unexp_err_fp = open(_UNEXPECTED_ERROR_FN.format(
                        np.random.randint(10000)), 'w')
                if len(failed_reads[err_type]) >= _MAX_NUM_UNEXP_ERRORS:
                    unexp_err_fp.close()
                else:
                    unexp_err_fp.write(
                        fast5_fn + '\n:::\n' + err_tb + '\n\n\n')
                    unexp_err_fp.flush()
        if do_update_prog:
            if not suppress_progress: bar.update(1)
            reads_called += 1
        if num_update_errors > 0:
            bar.write(prog_prefix + format_fail_summ(
                bar_header,
                [(len(fns), err) for err, fns in failed_reads.items()],
                reads_called, num_update_errors), file=sys.stderr)

        return reads_called, unexp_err_fp


    logger = logging.get_logger()
    logger.info('Processing reads.')
    reads_called = 0
    unexp_err_fp = None
    failed_reads = defaultdict(list)
    bar, prog_prefix, bar_header = prep_errors_bar(
        num_update_errors, None, suppress_progress)
    while True:
        try:
            try:
                (is_err, do_update_prog, err_type, fast5_fn,
                 err_tb) = failed_reads_q.get(block=False)
                reads_called, unexp_err_fp = update_prog(
                    reads_called, unexp_err_fp)
            except queue.Empty:
                # get total number of reads once all reads are enumerated
                if bar is not None and bar.total is None:
                    if getter_num_reads_conn.poll():
                        bar.total = getter_num_reads_conn.recv()
                else:
                    # if all reads are done signal was sent from main thread
                    if f_conn.poll():
                        break
                sleep(0.1)
                continue
        except KeyboardInterrupt:
            # exit gracefully on keyboard inturrupt
            return
    if not suppress_progress: bar.close()

    if len(failed_reads[_UNEXPECTED_ERROR_CODE]) >= 1:
        logger.warning((
            'Unexpected errors occured. See full ' +
            'error stack traces for first (up to) {0:d} errors in ' +
            '"{1}"').format(_MAX_NUM_UNEXP_ERRORS, unexp_err_fp.name))
    if any(len(fns) > 0 for fns in failed_reads.values()):
        logger.info(
            format_fail_summ(
                'Unsuccessful processing types:',
                [(len(fns), err) for err, fns in failed_reads.items()
                 if len(fns) > 0], reads_called))
        # TODO flag to output failed read names to file

    return


###############################
##### All read processing #####
###############################

def process_all_reads(
        fast5s_dir, recursive, num_reads, model_info, outputs, out_dir, bc_fmt,
        aligner, add_chr_ref, snps_data, num_ps, num_update_errors,
        suppress_progress, mods_info, db_safety, edge_buffer, pr_ref_filts):
    logger = logging.get_logger()
    logger.info('Preparing workers to process reads.')
    # read filename queue filler
    # Note no maxsize for this queue to compute total number of reads while
    # also not delaying read processing
    read_file_q = mp.Queue()
    num_reads_conn, getter_num_reads_conn = mp.Pipe()
    files_p = mp.Process(
        target=_fill_files_queue, args=(
            read_file_q, fast5s_dir, num_reads, recursive, num_ps,
            num_reads_conn),
        daemon=True)
    files_p.start()
    # progress and failed reads getter
    failed_reads_q, f_p, main_f_conn = mh.create_getter_q(
            _get_fail_queue, (getter_num_reads_conn, num_update_errors,
                              suppress_progress))

    # start output type getters/writers
    (bc_q, bc_p, main_bc_conn, mo_q, mo_p, main_mo_conn, snps_q, snps_p,
     main_snps_conn, mods_q, mods_p, main_mods_conn) = [None,] * 12
    if mh.BC_NAME in outputs or mh.BC_MODS_NAME in outputs:
        if mh.BC_NAME not in outputs:
            outputs.append(mh.BC_NAME)
        bc_q, bc_p, main_bc_conn = mh.create_getter_q(
            _get_bc_queue, (out_dir, bc_fmt, mods_info.do_output_mods,
                            mods_info.mod_long_names))
    if mh.MAP_NAME in outputs:
        do_output_pr_refs = (mh.PR_REF_NAME in outputs and
                             not mods_info.do_pr_ref_mods and
                             not snps_data.do_pr_ref_snps)
        mo_q, mo_p, main_mo_conn = mh.create_getter_q(
            mapping._get_map_queue, (
                out_dir, aligner.ref_names_and_lens, aligner.out_fmt,
                aligner.ref_fn, do_output_pr_refs, pr_ref_filts))
    if mh.PR_SNP_NAME in outputs:
        pr_refs_fn = os.path.join(out_dir, mh.PR_REF_FN) if (
            mh.PR_REF_NAME in outputs and
            snps_data.do_pr_ref_snps) else None
        snps_db_fn, snps_txt_fn = mh.OUTPUT_FNS[mh.PR_SNP_NAME]
        snps_txt_fn = (os.path.join(out_dir, snps_txt_fn)
                       if snps_data.write_snps_txt else None)
        snps_q, snps_p, main_snps_conn = mh.create_getter_q(
            snps._get_snps_queue, (
                snps_data.snp_id_tbl, os.path.join(out_dir, snps_db_fn),
                snps_txt_fn, db_safety, pr_refs_fn, pr_ref_filts))
    if mh.PR_MOD_NAME in outputs:
        pr_refs_fn = os.path.join(out_dir, mh.PR_REF_FN) if (
            mh.PR_REF_NAME in outputs and
            mods_info.do_pr_ref_mods) else None
        mods_db_fn, mods_txt_fn = mh.OUTPUT_FNS[mh.PR_MOD_NAME]
        mods_txt_fn = (os.path.join(out_dir, mods_txt_fn)
                       if mods_info.write_mods_txt else None)
        mods_q, mods_p, main_mods_conn = mh.create_getter_q(
            mods._get_mods_queue, (
                os.path.join(out_dir, mods_db_fn), mods_txt_fn, db_safety,
                pr_refs_fn, pr_ref_filts))

    proc_reads_ps, map_conns = [], []
    for device in model_info.process_devices:
        if aligner is None:
            map_conn, caller_conn = None, None
        else:
            map_conn, caller_conn = mp.Pipe()
        map_conns.append(map_conn)
        p = mp.Process(
            target=_process_reads_worker, args=(
                read_file_q, bc_q, snps_q, failed_reads_q, mods_q,
                caller_conn, model_info, snps_data.snps_to_test,
                snps_data.all_paths, snps_data.calib_table,
                snps_data.context_bases, mods_info, edge_buffer, device))
        p.daemon = True
        p.start()
        proc_reads_ps.append(p)
    sleep(0.1)

    # perform mapping in threads for mappy shared memory interface
    # open threads after all processes have started due to python
    # multiprocess combined with threading instability
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
                if on == mh.PR_SNP_NAME:
                    logger.info(
                        'Waiting for snps database to complete indexing.')
                elif on ==  mh.PR_MOD_NAME:
                    logger.info(
                        'Waiting for mods database to complete indexing.')
                p.join()
    except KeyboardInterrupt:
        logger.error('Exiting due to keyboard interrupt.')
        sys.exit(1)

    return


######################################
########## Input validation ##########
######################################

def aligner_validation(args):
    logger = logging.get_logger()
    if len(mh.ALIGN_OUTPUTS.intersection(args.outputs)) > 0:
        if args.reference is None:
            logger.error(
                ('Output(s) requiring reference alignment requested ({}), ' +
                 'but --reference not provided.').format(', '.join(
                    mh.ALIGN_OUTPUTS.intersection(args.outputs))))
            sys.exit(1)
        logger.info('Loading reference.')
        aligner = mapping.alignerPlus(
            str(args.reference), preset=str('map-ont'), best_n=1)
        setattr(aligner, 'out_fmt', args.mappings_format)
        setattr(aligner, 'ref_fn', args.reference)
        aligner.add_ref_names(args.reference)
    else:
        aligner = None
        if args.reference is not None:
            logger.warning(
                '[--reference] provided, but no [--outputs] requiring ' +
                'alignment was requested. Argument will be ignored.')
    return aligner

def snps_validation(args, is_cat_mod, output_size):
    logger = logging.get_logger()
    if mh.SNP_NAME in args.outputs and not mh.PR_SNP_NAME in args.outputs:
        args.outputs.append(mh.PR_SNP_NAME)
    if mh.PR_SNP_NAME in args.outputs and args.snp_filename is None:
        logger.error(
            '{} output requested, '.format(mh.PR_SNP_NAME) +
            'but --snp-filename provided.')
        sys.exit(1)
    if mh.PR_SNP_NAME in args.outputs and not (
            is_cat_mod or mh.nstate_to_nbase(output_size) == 4):
        logger.error(
            'SNP calling from naive modified base flip-flop model is ' +
            'not supported.')
        sys.exit(1)
    snp_calib_fn = mh.get_snp_calibration_fn(
        args.snp_calibration_filename, args.disable_snp_calibration)
    # snps data object loads with None snp_fn for easier handling downstream
    if args.max_snp_size > snps.HARD_MAX_SNP_SIZE:
        logger.warning((
            'Cannot process SNPs larger than {} bases. Re-setting max '+
            'SNP size.').format(snps.HARD_MAX_SNP_SIZE))
        args.max_snp_size = snps.HARD_MAX_SNP_SIZE
    snps_data = snps.SnpData(
        args.snp_filename, args.prepend_chr_vcf, args.max_snp_size,
        args.snp_all_paths, args.write_snps_text, args.snp_context_bases,
        snp_calib_fn,
        snps.HAPLIOD_MODE if args.haploid else snps.DIPLOID_MODE,
        args.refs_include_snps)
    if args.snp_filename is not None and mh.PR_SNP_NAME not in args.outputs:
        logger.warning(
            '--snps-filename provided, but SNP output not requested ' +
            '(via --outputs). Argument will be ignored.')
    return args, snps_data

def mods_validation(args, model_info):
    logger = logging.get_logger()
    if args.refs_include_mods and mh.PR_MOD_NAME not in args.outputs:
        # TODO don't really have to output this data, but have to compute it
        # so sort out how to compute the output but not output it
        args.outputs.append(mh.PR_MOD_NAME)
    if mh.PR_MOD_NAME not in args.outputs and mh.MOD_NAME in args.outputs:
        args.outputs.append(mh.PR_MOD_NAME)
    if mh.PR_MOD_NAME in args.outputs and not model_info.is_cat_mod:
        logger.error(
            '{} output requested, '.format(mh.PR_MOD_NAME) +
            'but model provided is not a categotical modified base model.\n' +
            'Note that modified base calling from naive modified base ' +
            'model is not currently supported.')
        sys.exit(1)
    if (model_info.is_cat_mod and mh.PR_MOD_NAME not in args.outputs and
        mh.BC_MODS_NAME not in args.outputs):
        logger.warning(
            ('Categorical modifications model provided, but neither {} nor ' +
            '{} requested (via --outputs). Modified base output will not be ' +
             'produced.').format( mh.PR_MOD_NAME, mh.BC_MODS_NAME))
    if args.mod_motif is not None and mh.PR_MOD_NAME not in args.outputs:
        logger.warning((
            '--mod-motif provided, but {} not requested (via --outputs). ' +
            'Argument will be ignored.').format(mh.PR_MOD_NAME))
    if args.refs_include_mods and mh.PR_REF_NAME not in args.outputs:
        logger.warning((
            '--refs-include-mods provided, but {} not requested ' +
            '(via --outputs). Argument will be ignored.').format(
                mh.PR_REF_NAME))
    mod_calib_fn = mh.get_mod_calibration_fn(
        args.mod_calibration_filename, args.disable_mod_calibration)
    mods_info = mods.ModInfo(
        model_info, args.mod_motif, args.mod_all_paths,
        args.write_mods_text, args.mod_context_bases,
        mh.BC_MODS_NAME in args.outputs, args.refs_include_mods, mod_calib_fn)
    return args, mods_info

def parse_pr_ref_output(args):
    logger = logging.get_logger()
    if args.output_per_read_references:
        args.outputs.append(mh.PR_REF_NAME)
        if args.refs_include_snps and args.refs_include_mods:
            logger.error('Cannot output both modified base and SNPs in ' +
                         'per-read references (remove one of ' +
                         '--refs-include-snps or --refs-include-mods).')
            sys.exit(1)
        if args.refs_include_snps and mh.PR_SNP_NAME not in args.outputs:
            args.outputs.append(mh.PR_SNP_NAME)
            logger.warning('--refs-include-snps set, so adding ' +
                           'per_read_snps to --outputs.')
        if args.refs_include_mods and mh.PR_MOD_NAME not in args.outputs:
            args.outputs.append(mh.PR_MOD_NAME)
            logger.warning('--refs-include-mods set, so adding ' +
                           'per_read_mods to --outputs.')
    else:
        if args.refs_include_snps:
            logger.warning(
                '--refs-include-snps but not --output-per-read-references ' +
                'set. Ignoring --refs-include-snps.')
        if args.refs_include_mods:
            logger.warning(
                '--refs-include-mods but not --output-per-read-references ' +
                'set. Ignoring --refs-include-mods.')
    min_len, max_len = (args.refs_length_range
                        if args.refs_length_range is not None else
                        (None, None))
    pr_ref_filts = mh.PR_REF_FILTERS(
        pct_idnt=args.refs_percent_identity_threshold,
        pct_cov=args.refs_percent_coverage_threshold,
        min_len=min_len, max_len=max_len)

    return args, pr_ref_filts

def mkdir(out_dir, overwrite):
    logger = logging.get_logger()
    if os.path.exists(out_dir):
        if not overwrite:
            sys.stderr.write(
                'ERROR: --output-directory exists and --overwrite is ' +
                'not set.\n')
            sys.exit(1)
        if os.path.isfile(out_dir) or os.path.islink(out_dir):
            os.remove(out_dir)
        else:
            shutil.rmtree(out_dir)
    os.mkdir(out_dir)

    return

def profile_validation(args):
    logger = logging.get_logger()
    if args.processes > 1:
        msg = ('Running profiling with multiple processes is ' +
               'not allowed. Setting to single process.')
        args.processes = 1
    else:
        msg = 'Running profiling. This may slow processing.'
    logger.warning(msg)
    return args


##########################
########## Main ##########
##########################

def get_parser():
    # hide more complex arguments for standard help output
    show_hidden_args = '--help-long' in sys.argv
    def hidden_help(help_msg):
        if not show_hidden_args:
            return argparse.SUPPRESS
        return help_msg

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
        help=hidden_help('Flappie model name.'))

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
        '--basecalls-format', choices=mh.BC_OUT_FMTS, default=mh.BC_OUT_FMTS[0],
        help=hidden_help('Basecalls output format. Choices: {}'.format(
            ', '.join(mh.BC_OUT_FMTS))))
    out_grp.add_argument(
        '--num-reads', type=int,
        help=hidden_help('Number of reads to process. Default: All reads'))

    map_grp = parser.add_argument_group('Mapping Arguments')
    map_grp.add_argument(
        '--mappings-format', choices=mh.MAP_OUT_FMTS,
        default=mh.MAP_OUT_FMTS[0],
        help='Mappings output format. Choices: {}'.format(
            ', '.join(mh.MAP_OUT_FMTS)))
    map_grp.add_argument(
        '--reference',
        help='Reference FASTA file used for mapping called reads.')

    map_grp.add_argument(
        '--prepend-chr-ref', action='store_true',
        help=hidden_help('Prepend "chr" to chromosome names from reference ' +
                         'to match VCF names.'))

    snp_grp = parser.add_argument_group('SNP Arguments')
    snp_grp.add_argument(
        '--haploid', action='store_true',
        help='Compute SNP aggregation for haploid genotypes. Default: diploid')
    snp_grp.add_argument(
        '--snp-filename',
        help='SNPs to call for each read in VCF format (required for output).')
    snp_grp.add_argument(
        '--write-snps-text', action='store_true',
        help='Write per-read SNP calls out to a text file. Default: ' +
        'Only ouput to database.')

    snp_grp.add_argument(
        '--disable-snp-calibration', action='store_true',
        help=hidden_help('Use raw SNP scores from the network. ' +
                         'Default: Calibrate score with ' +
                         '--snp-calibration-filename'))
    snp_grp.add_argument(
        '--heterozygous-factors', type=float, nargs=2,
        default=[mh.DEFAULT_SNV_HET_FACTOR, mh.DEFAULT_INDEL_HET_FACTOR],
        help=hidden_help('Bayesian prior factor for snv and indel ' +
                         'heterozygous calls (compared to 1.0 for hom ' +
                         'ref/alt). Default: %(default)s'))
    snp_grp.add_argument(
        '--max-snp-size', type=int, default=5,
        help=hidden_help('Maximum difference in number of reference and ' +
                         'alternate bases. Cannot exceed {}. '.format(
                             snps.HARD_MAX_SNP_SIZE) +
                         'Default: %(default)d'))
    snp_grp.add_argument(
        '--prepend-chr-vcf', action='store_true',
        help=hidden_help('Prepend "chr" to chromosome names from VCF to ' +
                         'match reference names.'))
    snp_grp.add_argument(
        '--snp-all-paths', action='store_true',
        help=hidden_help('Compute forwards algorithm all paths score. ' +
                         '(Default: Viterbi best-path score)'))
    snp_grp.add_argument(
        '--snp-calibration-filename',
        help=hidden_help('File containing emperical calibration for ' +
                         'SNP scores. As created by ' +
                         'megalodon/scripts/calibrate_snp_llr_scores.py. ' +
                         'Default: Load default calibration file.'))
    snp_grp.add_argument(
        '--snp-context-bases', type=int, nargs=2, default=[10, 30],
        help=hidden_help('Context bases for single base SNP and indel ' +
                         'calling. Default: %(default)s'))
    snp_grp.add_argument(
        '--write-vcf-llr', action='store_true',
        help=hidden_help('Write log-likelihood ratios out in non-standard ' +
                         'VCF field.'))

    mod_grp = parser.add_argument_group('Modified Base Arguments')
    mod_grp.add_argument(
        '--mod-motif', action="append", nargs=3,
        metavar=('base', 'motif', 'position'),
        help='Restrict modified base calls to specified motif(s). For ' +
        'example to restrict to CpG, dcm and dam motifs use ' +
        '"--mod-motif Z CG 0 --mod-motif Z CCWGG 1 --mod-motif Y GATC 1".')
    mod_grp.add_argument(
        '--write-mods-text', action='store_true',
        help='Write per-read modified bases out to a text file. Default: ' +
        'Only ouput to database.')

    mod_grp.add_argument(
        '--disable-mod-calibration', action='store_true',
        help=hidden_help('Use raw modified base scores from the network. ' +
                         'Default: Calibrate scores as described in ' +
                         '--mod-calibration-filename'))
    mod_grp.add_argument(
        '--mod-binary-threshold', type=float, nargs=2,
        default=mods.DEFAULT_AGG_INFO.binary_threshold,
        help=hidden_help('Thresholds for modified base aggregation. ' +
                         'Default: %(default)s'))
    mod_grp.add_argument(
        '--mod-calibration-filename',
        help=hidden_help('File containing emperical calibration for ' +
                         'modified base scores. As created by ' +
                         'megalodon/scripts/calibrate_mod_llr_scores.py. ' +
                         'Default: Load default calibration file.'))
    mod_grp.add_argument(
        '--mod-context-bases', type=int, default=10,
        help=hidden_help('Context bases for modified base calling. ' +
                         'Default: %(default)d'))
    mod_grp.add_argument(
        '--mod-all-paths', action='store_true',
        help=hidden_help('Compute forwards algorithm all paths score for ' +
                         'modified base calls. (Default: Viterbi ' +
                         'best-path score)'))

    tai_grp = parser.add_argument_group('Taiyaki Signal Chunking Arguments')
    tai_grp.add_argument(
        '--chunk_size', type=int, default=1000,
        help='Chunk length for base calling. Default: %(default)d')
    tai_grp.add_argument(
        '--chunk_overlap', type=int, default=100,
        help='Overlap between chunks to be stitched together. ' +
        'Default: %(default)d')
    tai_grp.add_argument(
        '--devices', type=int, nargs='+',
        help='CUDA GPU devices to use (only valid for taiyaki), default: CPU')
    tai_grp.add_argument(
        '--max_concurrent_chunks', type=int, default=50,
        help='Only process N chunks concurrently per-read (to avoid GPU ' +
        'memory errors). Default: %(default)d')

    refout_grp = parser.add_argument_group('Reference Output Arguments')
    refout_grp.add_argument(
        '--output-per-read-references', action='store_true',
        help=hidden_help('Output per-read references.'))
    refout_grp.add_argument(
        '--refs-include-mods', action='store_true',
        help=hidden_help('Include modified base calls in per-read ' +
                         'reference output.'))
    refout_grp.add_argument(
        '--refs-include-snps', action='store_true',
        help=hidden_help('Include SNP calls in per-read ' +
                         'reference output.'))
    refout_grp.add_argument(
        '--refs-percent-identity-threshold', type=float,
        help=hidden_help('Only include reads with higher percent identity ' +
                         'in per-read reference output.'))
    refout_grp.add_argument(
        '--refs-percent-coverage-threshold', type=float,
        help=hidden_help('Only include reads with higher read alignment ' +
                         'coverage in per-read reference output.'))
    refout_grp.add_argument(
        '--refs-length-range', type=int, nargs=2,
        help=hidden_help('Only include reads with specified read length ' +
                         'in per-read reference output.'))

    misc_grp = parser.add_argument_group('Miscellaneous Arguments')
    misc_grp.add_argument(
        '--help-long', help='Show all options.', action='help')
    misc_grp.add_argument(
        '--processes', type=int, default=1,
        help='Number of parallel processes. Default: %(default)d')
    misc_grp.add_argument(
        '--verbose-read-progress', type=int, default=0,
        help='Output verbose output on read progress. Outputs N most ' +
        'common points where reads could not be processed further. ' +
        'Default: %(default)d')
    misc_grp.add_argument(
        '-v', '--version', action='version',
        version='Megalodon version: {}'.format(MEGALODON_VERSION),
        help='show megalodon version and exit.')

    misc_grp.add_argument(
        '--database-safety', type=int, default=1,
        help=hidden_help('Setting for database performance versus ' +
                         'corruption protection. Options: 0 (DB corruption ' +
                         'on application crash), 1 (DB corruption on system ' +
                         'crash), 2 (DB safe mode). Default: %(default)d'))
    misc_grp.add_argument(
        '--edge-buffer', type=int, default=100,
        help=hidden_help('Ignore SNP or indel calls near edge of read ' +
                         'mapping. Default: %(default)d'))
    misc_grp.add_argument(
        '--not-recursive', action='store_true',
        help=hidden_help('Only search for fast5 read files directly found ' +
                         'within the fast5 directory. Default: search ' +
                         'recursively'))
    misc_grp.add_argument(
        '--suppress-progress', action='store_true',
        help=hidden_help('Suppress progress bar output.'))

    return parser


def _main():
    args = get_parser().parse_args()

    mkdir(args.output_directory, args.overwrite)
    logging.init_logger(args.output_directory)
    logger = logging.get_logger()
    logger.debug('Command: """' + ' '.join(sys.argv) + '"""')

    if _DO_PROFILE:
        args = profile_validation(args)

    args, pr_ref_filts = parse_pr_ref_output(args)
    model_info = backends.ModelInfo(
        args.flappie_model_name, args.taiyaki_model_filename, args.devices,
        args.processes, args.chunk_size, args.chunk_overlap,
        args.max_concurrent_chunks)
    args, mods_info = mods_validation(args, model_info)
    args, snps_data = snps_validation(
        args, model_info.is_cat_mod, model_info.output_size)
    aligner = aligner_validation(args)

    process_all_reads(
        args.fast5s_dir, not args.not_recursive, args.num_reads, model_info,
        args.outputs, args.output_directory, args.basecalls_format, aligner,
        args.prepend_chr_ref, snps_data, args.processes,
        args.verbose_read_progress, args.suppress_progress,
        mods_info, args.database_safety, args.edge_buffer, pr_ref_filts)

    if mh.SNP_NAME in args.outputs or mh.MOD_NAME in args.outputs:
        mod_names = (mods_info.mod_long_names
                     if mh.MOD_NAME in args.outputs else [])
        mod_agg_info = mods.AGG_INFO(
            mods.BIN_THRESH_NAME, args.mod_binary_threshold)
        aggregate.aggregate_stats(
            args.outputs, args.output_directory, args.processes,
            args.write_vcf_llr, args.heterozygous_factors, snps_data.call_mode,
            mod_names, mod_agg_info, args.suppress_progress)

    return

if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
