#!/usr/bin/env python3
import os
# set blas library environment variables (without these the cblas calls
# can completely halt processing)
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'

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
    aggregate, backends, decode, fast5_io, logging, mapping, mods,
    variants, megalodon_helper as mh)
from megalodon._version import MEGALODON_VERSION


_DO_PROFILE = False
_UNEXPECTED_ERROR_CODE = 'Unexpected error'
_UNEXPECTED_ERROR_FN = 'unexpected_megalodon_errors.{}.err'
_MAX_NUM_UNEXP_ERRORS = 50


###########################
##### Read Processing #####
###########################

def handle_errors(func, args, r_vals, out_q, fast5_fn, failed_reads_q):
    try:
        out_q.put((func(*args), r_vals))
    except KeyboardInterrupt:
        failed_reads_q.put(
            (True, False, 'Keyboard interrupt', fast5_fn, None, 0))
        return
    except mh.MegaError as e:
        failed_reads_q.put((True, False, str(e), fast5_fn, None, 0))
    except:
        failed_reads_q.put((
            True, False, _UNEXPECTED_ERROR_CODE, fast5_fn,
            traceback.format_exc(), 0))
    return

def process_read(
        sig_info, model_info, bc_q, caller_conn, sig_map_q, sig_map_info,
        vars_data, vars_q, mods_q, mods_info, failed_reads_q):
    """ Workhorse per-read megalodon function (connects all the parts)
    """
    r_post, mod_weights, can_nmods = model_info.get_posteriors(sig_info)

    if mods_q is not None:
        r_post_w_mods = np.concatenate([r_post, mod_weights], axis=1)
    if not mods_info.do_output_mods:
        mod_weights = None
    r_seq, score, rl_cumsum, mods_scores = decode.decode_post(
        r_post, mods_info.alphabet, mod_weights, can_nmods)
    if bc_q is not None:
        bc_q.put((sig_info.read_id, r_seq, mods_scores))

    # if no mapping connection return after basecalls are passed out
    if caller_conn is None: return

    # map read and record mapping from reference to query positions
    r_ref_seq, r_to_q_poss, r_ref_pos, r_cigar = mapping.map_read(
        r_seq, sig_info.read_id, caller_conn)
    np_ref_seq = mh.seq_to_int(r_ref_seq)

    sig_map_res = None
    if sig_map_q is not None:
        pass_sig_map_filts = mapping.read_passes_filters(
            sig_map_info, len(r_seq), r_ref_pos.q_trim_start,
            r_ref_pos.q_trim_end, r_cigar)
        sig_map_res = (
            pass_sig_map_filts, sig_info.fast5_fn, sig_info.dacs,
            sig_info.scale_params, r_ref_seq, sig_info.stride,
            sig_map_info.alphabet, sig_info.read_id, r_to_q_poss, rl_cumsum,
            r_ref_pos.q_trim_start)
        if not sig_map_info.annotate_mods and pass_sig_map_filts:
            try:
                sig_map_q.put(signal_mapping.get_remapping(*sig_map_res[1:]))
            except Exception as e:
                logger = logging.get_logger()
                logger.debug((
                    'Read: {} {} failed mapped signal validation with ' +
                    'error: {}').format(
                        sig_info.fast5_fn, sig_info.read_id, str(e)))
                # taiyaki errors can contain newlines so split them here
                failed_reads_q.put((
                    True, False, ' ::: '.join(str(e).strip().split('\n')),
                    sig_info.fast5_fn, None, 0))

    # get mapped start in post and run len to mapped bit of output
    post_mapped_start, post_mapped_end = (rl_cumsum[r_ref_pos.q_trim_start],
                                          rl_cumsum[r_ref_pos.q_trim_end])
    mapped_rl_cumsum = rl_cumsum[
        r_ref_pos.q_trim_start:r_ref_pos.q_trim_end + 1] - post_mapped_start
    ref_to_block = np.array([
        mapped_rl_cumsum[query_pos]
        for ref_pos, query_pos in r_to_q_poss.items()])

    if vars_q is not None:
        mapped_r_post = r_post[post_mapped_start:post_mapped_end]
        handle_errors(
            func=variants.call_read_vars,
            args=(vars_data, r_ref_pos, np_ref_seq, ref_to_block,
                  mapped_r_post),
            r_vals=(sig_info.read_id, r_ref_pos.chrm, r_ref_pos.strand,
                    r_ref_pos.start, r_ref_seq, len(r_seq),
                    r_ref_pos.q_trim_start, r_ref_pos.q_trim_end, r_cigar),
            out_q=vars_q, fast5_fn=sig_info.fast5_fn + ':::' + sig_info.read_id,
            failed_reads_q=failed_reads_q)
    if mods_q is not None:
        mapped_r_post_w_mods = r_post_w_mods[post_mapped_start:post_mapped_end]
        mod_sig_map_q = sig_map_q if sig_map_info.annotate_mods else None
        handle_errors(
            func=mods.call_read_mods,
            args=(r_ref_pos, r_ref_seq, ref_to_block, mapped_r_post_w_mods,
                  mods_info, mod_sig_map_q, sig_map_res),
            r_vals=(sig_info.read_id, r_ref_pos.chrm, r_ref_pos.strand,
                    r_ref_pos.start, r_ref_seq, len(r_seq),
                    r_ref_pos.q_trim_start, r_ref_pos.q_trim_end, r_cigar),
            out_q=mods_q, fast5_fn=sig_info.fast5_fn + ':::' + sig_info.read_id,
            failed_reads_q=failed_reads_q)

    return


############################
##### Multi-processing #####
############################

def _get_bc_queue(
        bc_q, bc_conn, out_dir, bc_fmt, do_output_mods, mod_long_names):
    bc_fp = open(mh.get_megalodon_fn(out_dir, mh.BC_NAME) + '.' + bc_fmt, 'w')
    if do_output_mods:
        mods_fp = h5py.File(mh.get_megalodon_fn(out_dir, mh.BC_MODS_NAME))
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
            sleep(0.001)
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
        read_file_q, bc_q, vars_q, failed_reads_q, mods_q, caller_conn,
        sig_map_q, sig_map_info, model_info, vars_data, mods_info, device):
    # wrap process prep in try loop to avoid stalled command
    logger = logging.get_logger('main')
    try:
        model_info.prep_model_worker(device)
        vars_data.reopen_variant_index()
        logger.debug('Starting read worker {}'.format(mp.current_process()))
    except:
        if caller_conn is not None:
            caller_conn.send(True)
        logger.debug(('Read worker {} has failed process preparation.\n' +
                      'Full error traceback:\n{}').format(
                          mp.current_process(), traceback.format_exc()))
        return

    while True:
        try:
            try:
                fast5_fn, read_id = read_file_q.get(block=False)
            except queue.Empty:
                sleep(0.001)
                continue

            if fast5_fn is None:
                if caller_conn is not None:
                    caller_conn.send(True)
                logger.debug('Gracefully exiting read worker {}'.format(
                    mp.current_process()))
                break
            logger.debug('Analyzing read {}'.format(read_id))
            sig_info = model_info.extract_signal_info(
                fast5_fn, read_id, sig_map_q is not None)
            process_read(
                sig_info, model_info, bc_q, caller_conn, sig_map_q,
                sig_map_info, vars_data, vars_q, mods_q, mods_info,
                failed_reads_q)
            failed_reads_q.put((
                False, True, None, None, None, sig_info.raw_len))
            logger.debug('Successfully processed read {}'.format(read_id))
        except KeyboardInterrupt:
            failed_reads_q.put((
                True, True, 'Keyboard interrupt', fast5_fn, None, 0))
            logger.debug('Keyboard interrupt during read {}'.format(read_id))
            return
        except mh.MegaError as e:
            failed_reads_q.put((
                True, True, str(e), fast5_fn + ':::' + read_id, None,
                sig_info.raw_len))
            logger.debug('Incomplete processing for read {} ::: {}'.format(
                read_id, str(e)))
        except:
            failed_reads_q.put((
                True, True, _UNEXPECTED_ERROR_CODE, fast5_fn + ':::' + read_id,
                traceback.format_exc(), 0))
            logger.debug('Unexpected error for read {}'.format(read_id))

    return

if _DO_PROFILE:
    _process_reads_wrapper = _process_reads_worker
    def _process_reads_worker(*args):
        import cProfile
        cProfile.runctx('_process_reads_wrapper(*args)', globals(), locals(),
                        filename='read_processing.prof')
        return


####################################
##### Post Per-read Processing #####
####################################

def post_process_whatshap(out_dir, map_fmt, ref_fn):
    whatshap_map_bn = mh.get_megalodon_fn(out_dir, mh.WHATSHAP_MAP_NAME)
    whatshap_map_fn = whatshap_map_bn + '.' + map_fmt
    whatshap_sort_fn = whatshap_map_bn + '.sorted.bam'
    whatshap_p = mp.Process(
        target=mapping.sort_and_index_mapping,
        args=(whatshap_map_fn, whatshap_sort_fn, ref_fn), daemon=True)
    whatshap_p.start()
    sleep(0.001)

    return whatshap_sort_fn, whatshap_p

def post_process_mapping(out_dir, map_fmt, ref_fn):
    map_bn = mh.get_megalodon_fn(out_dir, mh.MAP_NAME)
    map_fn = map_bn + '.' + map_fmt
    map_sort_fn = map_bn + '.sorted.bam'
    map_p = mp.Process(
        target=mapping.sort_and_index_mapping,
        args=(map_fn, map_sort_fn, ref_fn), daemon=True)
    map_p.start()
    sleep(0.001)

    return map_p

def post_process_aggregate(
        mods_info, outputs, out_dir, num_ps, write_vcf_lp,
        het_factors, vars_data, write_mod_lp, supp_prog):
    mod_names = mods_info.mod_long_names if mh.MOD_NAME in outputs else []
    aggregate.aggregate_stats(
        outputs, out_dir, num_ps, write_vcf_lp, het_factors,
        vars_data.call_mode, mod_names, mods_info.agg_info,
        write_mod_lp, mods_info.mod_output_fmts, supp_prog)
    return



##################################
##### Dynamic error updating #####
##################################

def _fill_files_queue(
        read_file_q, fast5s_dir, num_reads, read_ids_fn, recursive, num_ps,
        num_reads_conn):
    logger = logging.get_logger()
    valid_read_ids = None
    if read_ids_fn is not None:
        with open(read_ids_fn) as read_ids_fp:
            valid_read_ids = set(line.strip() for line in read_ids_fp)
    used_read_ids = set()
    # fill queue with read filename and read id tuples
    for fast5_fn, read_id in fast5_io.iterate_fast5_reads(
            fast5s_dir, num_reads, recursive):
        if valid_read_ids is not None and read_id not in valid_read_ids:
            continue
        if read_id in used_read_ids:
            logger.debug(
                ('Read ID ({}) found in previous read and will not ' +
                 'process from {}.').format(read_id, fast5_fn))
            continue
        if fast5_fn is None or read_id is None:
            continue
        read_file_q.put((fast5_fn, read_id))
        used_read_ids.add(read_id)
    # add None to indicate that read processes should return
    for _ in range(num_ps):
        read_file_q.put((None, None))
    num_reads_conn.send(len(used_read_ids))

    return

def format_fail_summ(header, fail_summ=[], reads_called=0, num_errs=None):
    summ_errs = sorted(fail_summ)[::-1]
    if num_errs is not None:
        summ_errs = summ_errs[:num_errs]
        if len(summ_errs) < num_errs:
            summ_errs.extend([(None, '') for _ in range(
                num_errs - len(summ_errs))])
    errs_str = '\n'.join("{:8.1f}% ({:>7} reads)".format(
        100 * n_fns / float(reads_called), n_fns) + " : " + '{:<80}'.format(err)
        if (n_fns is not None and reads_called > 0) else
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
        bar = tqdm(total=tot_reads, smoothing=0, initial=curr_num_reads,
                   unit='read', dynamic_ncols=True)
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
    def update_prog(reads_called, sig_called, unexp_err_fp):
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
            if not suppress_progress:
                try:
                    bar.set_postfix({
                        'ksamp/s':(sig_called / 1000) /
                        bar.format_dict['elapsed']})
                except AttributeError:
                    # sometimes get no format_dict error
                    # so don't include ksample/s if so
                    pass
                bar.update(1)
                if num_update_errors > 0:
                    bar.write(prog_prefix + format_fail_summ(
                        bar_header,
                        [(len(fns), err) for err, fns in failed_reads.items()],
                        reads_called, num_update_errors), file=sys.stderr)
            reads_called += 1

        return reads_called, unexp_err_fp


    logger = logging.get_logger()
    logger.info('Processing reads.')
    reads_called, sig_called = 0, 0
    unexp_err_fp = None
    failed_reads = defaultdict(list)
    bar, prog_prefix, bar_header = prep_errors_bar(
        num_update_errors, None, suppress_progress)
    while True:
        try:
            try:
                (is_err, do_update_prog, err_type, fast5_fn,
                 err_tb, n_sig) = failed_reads_q.get(block=False)
                sig_called += n_sig
                reads_called, unexp_err_fp = update_prog(
                    reads_called, sig_called, unexp_err_fp)
            except queue.Empty:
                # get total number of reads once all reads are enumerated
                if bar is not None and bar.total is None:
                    if getter_num_reads_conn.poll():
                        bar.total = getter_num_reads_conn.recv()
                else:
                    # if all reads are done signal was sent from main thread
                    if f_conn.poll():
                        break
                sleep(0.001)
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
    else:
        logger.info('All reads processed successfully.')

    return


###############################
##### All read processing #####
###############################

def process_all_reads(
        fast5s_dir, recursive, num_reads, read_ids_fn, model_info, outputs,
        out_dir, bc_fmt, aligner, vars_data, num_ps, num_update_errors,
        suppress_progress, mods_info, db_safety, pr_ref_filts, sig_map_info):
    logger = logging.get_logger()
    logger.info('Preparing workers to process reads.')
    # read filename queue filler
    # Note no maxsize for this queue to compute total number of reads while
    # also not delaying read processing
    read_file_q = mp.Queue()
    num_reads_conn, getter_num_reads_conn = mp.Pipe()
    files_p = mp.Process(
        target=_fill_files_queue, args=(
            read_file_q, fast5s_dir, num_reads, read_ids_fn, recursive,
            num_ps, num_reads_conn),
        daemon=True)
    files_p.start()
    # progress and failed reads getter (no limit on failed reads queue
    # in case error occurs there, don't halt run
    failed_reads_q, f_p, main_f_conn = mh.create_getter_q(
            _get_fail_queue, (getter_num_reads_conn, num_update_errors,
                              suppress_progress), max_size=None)

    # start output type getters/writers
    (bc_q, bc_p, main_bc_conn, mo_q, mo_p, main_mo_conn, vars_q, vars_p,
     main_vars_conn, mods_q, mods_p, main_mods_conn, sig_map_q, sig_map_p,
     sig_map_conn) = [None,] * 15
    if mh.BC_NAME in outputs or mh.BC_MODS_NAME in outputs:
        if mh.BC_NAME not in outputs:
            outputs.append(mh.BC_NAME)
        bc_q, bc_p, main_bc_conn = mh.create_getter_q(
            _get_bc_queue, (out_dir, bc_fmt, mods_info.do_output_mods,
                            mods_info.mod_long_names))
    if mh.MAP_NAME in outputs:
        do_output_pr_refs = (mh.PR_REF_NAME in outputs and
                             not mods_info.do_pr_ref_mods and
                             not vars_data.do_pr_ref_vars)
        mo_q, mo_p, main_mo_conn = mh.create_getter_q(
            mapping._get_map_queue, (
                out_dir, aligner.ref_names_and_lens, aligner.out_fmt,
                aligner.ref_fn, do_output_pr_refs, pr_ref_filts))
    if mh.PR_VAR_NAME in outputs:
        pr_refs_fn = mh.get_megalodon_fn(out_dir, mh.PR_REF_NAME) if (
            mh.PR_REF_NAME in outputs and vars_data.do_pr_ref_vars) else None
        whatshap_map_fn = (
            mh.get_megalodon_fn(out_dir, mh.WHATSHAP_MAP_NAME) + '.' +
            aligner.out_fmt) if mh.WHATSHAP_MAP_NAME in outputs else None
        vars_txt_fn = (mh.get_megalodon_fn(out_dir, mh.PR_VAR_TXT_NAME)
                       if vars_data.write_vars_txt else None)
        vars_q, vars_p, main_vars_conn = mh.create_getter_q(
            variants._get_variants_queue, (
                mh.get_megalodon_fn(out_dir, mh.PR_VAR_NAME),
                vars_txt_fn, db_safety, pr_refs_fn, pr_ref_filts,
                whatshap_map_fn, aligner.ref_names_and_lens, aligner.ref_fn,
                vars_data.loc_index_in_memory))
    if mh.PR_MOD_NAME in outputs:
        pr_refs_fn = mh.get_megalodon_fn(out_dir, mh.PR_REF_NAME) if (
            mh.PR_REF_NAME in outputs and mods_info.do_pr_ref_mods) else None
        mods_txt_fn = (mh.get_megalodon_fn(out_dir, mh.PR_MOD_TXT_NAME)
                       if mods_info.write_mods_txt else None)
        mods_q, mods_p, main_mods_conn = mh.create_getter_q(
            mods._get_mods_queue, (
                mh.get_megalodon_fn(out_dir, mh.PR_MOD_NAME), db_safety,
                aligner.ref_names_and_lens, mods_txt_fn,
                pr_refs_fn, pr_ref_filts, mods_info.pos_index_in_memory))
    if mh.SIG_MAP_NAME in outputs:
        alphabet_info = signal_mapping.get_alphabet_info(model_info)
        sig_map_fn = mh.get_megalodon_fn(out_dir, mh.SIG_MAP_NAME)
        sig_map_q, sig_map_p, sig_map_conn = mh.create_getter_q(
            signal_mapping.write_signal_mappings,
            (sig_map_fn, alphabet_info))

    proc_reads_ps, map_conns = [], []
    for device in model_info.process_devices:
        if aligner is None:
            map_conn, caller_conn = None, None
        else:
            map_conn, caller_conn = mp.Pipe()
        map_conns.append(map_conn)
        p = mp.Process(
            target=_process_reads_worker, args=(
                read_file_q, bc_q, vars_q, failed_reads_q, mods_q, caller_conn,
                sig_map_q, sig_map_info, model_info, vars_data, mods_info,
                device))
        p.daemon = True
        p.start()
        proc_reads_ps.append(p)
    # ensure process all start up before initializing mapping threads
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
                args=(aligner, map_conn, mo_q))
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
                (mh.SIG_MAP_NAME, sig_map_p, sig_map_conn),
                (mh.PR_VAR_NAME, vars_p, main_vars_conn),
                (mh.PR_MOD_NAME, mods_p, main_mods_conn)):
            if on in outputs and p.is_alive():
                main_conn.send(True)
                if on == mh.PR_VAR_NAME:
                    logger.info(
                        'Waiting for variants database to complete indexing.')
                elif on ==  mh.PR_MOD_NAME:
                    logger.info(
                        'Waiting for mods database to complete indexing.')
                p.join()
    except KeyboardInterrupt:
        logger.error('Exiting due to keyboard interrupt.')
        sys.exit(1)

    return


############################
##### Input validation #####
############################

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
        if not (os.path.exists(args.reference) and
                os.path.isfile(args.reference)):
            logger.error('Provided reference file does not exist or is ' +
                         'not a file.')
            sys.exit(1)
        aligner = mapping.alignerPlus(
            str(args.reference), preset=str('map-ont'), best_n=1)
        setattr(aligner, 'out_fmt', args.mappings_format)
        setattr(aligner, 'ref_fn', mh.resolve_path(args.reference))
        aligner.add_ref_lens()
        mapping.test_open_alignment_out_file(
            args.output_directory, aligner.out_fmt,
            aligner.ref_names_and_lens, aligner.ref_fn)
    else:
        aligner = None
        if args.reference is not None:
            logger.warning(
                '[--reference] provided, but no [--outputs] requiring ' +
                'alignment was requested. Argument will be ignored.')
    return aligner

def vars_validation(args, is_cat_mod, output_size, aligner):
    logger = logging.get_logger()
    if mh.WHATSHAP_MAP_NAME in args.outputs and not mh.VAR_NAME in args.outputs:
        args.outputs.append(mh.VAR_NAME)
    if mh.VAR_NAME in args.outputs and not mh.PR_VAR_NAME in args.outputs:
        args.outputs.append(mh.PR_VAR_NAME)
    if mh.PR_VAR_NAME in args.outputs and args.variant_filename is None:
        logger.error(
            '{} output requested, '.format(mh.PR_VAR_NAME) +
            'but --variant-filename not provided.')
        sys.exit(1)
    if mh.PR_VAR_NAME in args.outputs and not (
            is_cat_mod or mh.nstate_to_nbase(output_size) == 4):
        logger.error(
            'Variant calling from naive modified base flip-flop model is ' +
            'not supported.')
        sys.exit(1)
    var_calib_fn = mh.get_var_calibration_fn(
        args.variant_calibration_filename, args.disable_variant_calibration)
    try:
        vars_data = variants.VarData(
            args.variant_filename, args.max_indel_size,
            args.variant_all_paths, args.write_variants_text,
            args.variant_context_bases, var_calib_fn,
            variants.HAPLIOD_MODE if args.haploid else variants.DIPLOID_MODE,
            args.refs_include_variants, aligner, edge_buffer=args.edge_buffer,
            context_min_alt_prob=args.context_min_alt_prob,
            loc_index_in_memory=not args.variant_locations_on_disk)
    except mh.MegaError as e:
        logger.error(str(e))
        sys.exit(1)
    if args.variant_filename is not None and mh.PR_VAR_NAME not in args.outputs:
        logger.warning(
            '--variants-filename provided, but variants output not requested ' +
            '(via --outputs). Argument will be ignored.')
    return args, vars_data

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
        args.mod_motif = None
    if args.refs_include_mods and mh.PR_REF_NAME not in args.outputs:
        logger.warning((
            '--refs-include-mods provided, but {} not requested ' +
            '(via --outputs). Argument will be ignored.').format(
                mh.PR_REF_NAME))
    mod_calib_fn = mh.get_mod_calibration_fn(
        args.mod_calibration_filename, args.disable_mod_calibration)
    if args.mod_aggregate_method == mods.EM_NAME:
        agg_info = mods.AGG_INFO(mods.EM_NAME, None)
    elif args.mod_aggregate_method == mods.BIN_THRESH_NAME:
        agg_info = mods.AGG_INFO(
            mods.BIN_THRESH_NAME, args.mod_binary_threshold)
    mods_info = mods.ModInfo(
        model_info, args.mod_motif, args.mod_all_paths,
        args.write_mods_text, args.mod_context_bases,
        mh.BC_MODS_NAME in args.outputs, args.refs_include_mods, mod_calib_fn,
        args.mod_output_formats, args.edge_buffer,
        not args.mod_positions_on_disk, agg_info)
    return args, mods_info

def parse_pr_ref_output(args):
    logger = logging.get_logger()
    if args.output_per_read_references:
        args.outputs.append(mh.PR_REF_NAME)
        if args.refs_include_vars and args.refs_include_mods:
            logger.error('Cannot output both modified base and variants in ' +
                         'per-read references (remove one of ' +
                         '--refs-include-variants or --refs-include-mods).')
            sys.exit(1)
        if args.refs_include_variants and mh.PR_VAR_NAME not in args.outputs:
            args.outputs.append(mh.PR_VAR_NAME)
            logger.warning('--refs-include-variants set, so adding ' +
                           'per_read_variants to --outputs.')
        if args.refs_include_mods and mh.PR_MOD_NAME not in args.outputs:
            args.outputs.append(mh.PR_MOD_NAME)
            logger.warning('--refs-include-mods set, so adding ' +
                           'per_read_mods to --outputs.')
    else:
        if args.refs_include_variants:
            logger.warning(
                '--refs-include-variantss but not ' +
                '--output-per-read-references set. Ignoring ' +
                '--refs-include-variants.')
        if args.refs_include_mods:
            logger.warning(
                '--refs-include-mods but not --output-per-read-references ' +
                'set. Ignoring --refs-include-mods.')
    min_len, max_len = (args.refs_length_range
                        if args.refs_length_range is not None else
                        (None, None))
    pr_ref_filts = mh.PR_REF_INFO(
        pct_idnt=args.refs_percent_identity_threshold,
        pct_cov=args.refs_percent_coverage_threshold,
        min_len=min_len, max_len=max_len)

    return args, pr_ref_filts

def parse_sig_map_output(args, model_info):
    logger = logging.get_logger()
    if args.output_signal_mappings:
        from megalodon import signal_mapping
        global signal_mapping
        sig_map_alphabet = signal_mapping.get_alphabet_info(model_info).alphabet
        args.outputs.append(mh.SIG_MAP_NAME)
        if args.signal_map_include_mods and mh.PR_MOD_NAME not in args.outputs:
            args.outputs.append(mh.PR_MOD_NAME)
            logger.warning('--signal-map-include-mods set, so adding ' +
                           '"per_read_mods" to --outputs.')
    else:
        sig_map_alphabet = None
        if args.signal_map_include_mods:
            logger.warning(
                '--signal-map-include-mods but not --output-signal-mappings ' +
                'set. Ignoring --signal-map-include-mods.')
    min_len, max_len = (args.signal_map_length_range
                        if args.signal_map_length_range is not None else
                        (None, None))

    sig_map_info = mh.PR_REF_INFO(
        pct_idnt=args.signal_map_percent_identity_threshold,
        pct_cov=args.signal_map_percent_coverage_threshold,
        min_len=min_len, max_len=max_len, alphabet=sig_map_alphabet,
        annotate_mods=args.signal_map_include_mods)

    return args, sig_map_info

def mkdir(out_dir, overwrite):
    logger = logging.get_logger()
    if os.path.exists(out_dir):
        if not overwrite:
            logger.error(
                '--output-directory exists and --overwrite is not set.')
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

class SelectiveRawFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        # special splitlines command for options output for better readability
        if text.startswith('O|'):
            return text[2:].splitlines()
        # else use standard RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)

def get_parser():
    # hide more complex arguments for standard help output
    show_hidden_args = '--help-long' in sys.argv
    def hidden_help(help_msg):
        if not show_hidden_args:
            return argparse.SUPPRESS
        return help_msg

    parser = argparse.ArgumentParser(formatter_class=SelectiveRawFormatter)
    parser.add_argument(
        'fast5s_dir',
        help='Directory containing raw fast5 (will be searched recursively).')

    mdl_grp = parser.add_argument_group('Model Arguments')
    mdl_grp.add_argument(
        '--taiyaki-model-filename',
        help='Taiyaki basecalling model checkpoint file.')
    mdl_grp.add_argument(
        '--load-default-model', action='store_true',
        help=('Load the default basecalling model included with megalodon ' +
              '({}). Default: Assume guppy --post_out FAST5 files as ' +
              'input').format(mh.MODEL_PRESET_DESC))
    mdl_grp.add_argument(
        '--devices', nargs='+',
        help='GPU devices for taiyaki basecalling backend (--processes will ' +
        'be distributed even over specified --devices).')
    mdl_grp.add_argument(
        '--chunk-size', type=int, default=1000,
        help=hidden_help('Chunk length for base calling. Default: %(default)d'))
    mdl_grp.add_argument(
        '--chunk-overlap', type=int, default=100,
        help=hidden_help('Overlap between chunks to be stitched together. ' +
                         'Default: %(default)d'))
    mdl_grp.add_argument(
        '--max-concurrent-chunks', type=int, default=200,
        help=hidden_help('Only process N chunks concurrently per-read (to ' +
                         'avoid GPU memory errors). Default: %(default)d'))


    out_grp = parser.add_argument_group('Output Arguments')
    out_grp.add_argument(
        '--outputs', nargs='+',
        default=['basecalls',], choices=tuple(mh.OUTPUT_DESCS.keys()),
        # note 'O|' triggers raw formatting for this option alone
        help='O|Desired output(s).\nOptions:\n' +
        '\n'.join(('\t{}: {}'.format(*out_desc)
                   for out_desc in mh.OUTPUT_DESCS.items())) +
        '\nDefault: %(default)s')
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
    out_grp.add_argument(
        '--read-ids-filename',
        help=hidden_help('File containing read ids to process (one per ' +
                         'line). Default: All reads'))

    map_grp = parser.add_argument_group('Mapping Arguments')
    map_grp.add_argument(
        '--mappings-format', choices=mh.MAP_OUT_FMTS,
        default=mh.MAP_OUT_FMTS[0],
        help='Mappings output format. Choices: {}'.format(
            ', '.join(mh.MAP_OUT_FMTS)))
    map_grp.add_argument(
        '--reference',
        help='Reference FASTA or minimap2 index file used for mapping ' +
        'called reads.')

    var_grp = parser.add_argument_group('Sequence Variant Arguments')
    var_grp.add_argument(
        '--haploid', action='store_true',
        help='Compute variant aggregation for haploid genotypes. ' +
        'Default: diploid')
    var_grp.add_argument(
        '--variant-filename',
        help='Sequence variants to call for each read in VCF/BCF format ' +
        '(required for variant output).')
    var_grp.add_argument(
        '--write-variants-text', action='store_true',
        help='Write per-read sequence variant calls out to a text file. ' +
        'Default: Only ouput to database.')

    var_grp.add_argument(
        '--context-min-alt-prob', type=float,
        default=mh.DEFAULT_CONTEXT_MIN_ALT_PROB,
        help=hidden_help('Minimum alternative alleles probability to ' +
                         'include variant in computation of nearby variants. ' +
                         'Default: %(default)f'))
    var_grp.add_argument(
        '--disable-variant-calibration', action='store_true',
        help=hidden_help('Use raw variant scores from the network. ' +
                         'Default: Calibrate score with ' +
                         '--variant-calibration-filename'))
    var_grp.add_argument(
        '--heterozygous-factors', type=float, nargs=2,
        default=[mh.DEFAULT_SNV_HET_FACTOR, mh.DEFAULT_INDEL_HET_FACTOR],
        help=hidden_help('Bayesian prior factor for snv and indel ' +
                         'heterozygous calls (compared to 1.0 for hom ' +
                         'ref/alt). Default: %(default)s'))
    var_grp.add_argument(
        '--max-indel-size', type=int, default=50,
        help=hidden_help('Maximum difference in number of reference and ' +
                         'alternate bases. Default: %(default)d'))
    var_grp.add_argument(
        '--variant-all-paths', action='store_true',
        help=hidden_help('Compute forwards algorithm all paths score. ' +
                         '(Default: Viterbi best-path score)'))
    var_grp.add_argument(
        '--variant-calibration-filename',
        help=hidden_help('File containing emperical calibration for ' +
                         'variant scores. As created by ' +
                         'megalodon/scripts/calibrate_variant_llr_scores.py. ' +
                         'Default: Load default calibration file.'))
    var_grp.add_argument(
        '--variant-locations-on-disk', action='store_true',
        help=hidden_help('Force sequence variant locations to be stored ' +
                         'only within on disk database table. This option ' +
                         'will reduce the RAM memory requirement, but may ' +
                         'drastically slow processing. Default: Store ' +
                         'locations in memory and on disk.'))
    var_grp.add_argument(
        '--variant-context-bases', type=int, nargs=2,
        default=[mh.DEFAULT_SNV_CONTEXT, mh.DEFAULT_INDEL_CONTEXT],
        help=hidden_help('Context bases for single base variant and indel ' +
                         'calling. Default: %(default)s'))
    var_grp.add_argument(
        '--write-vcf-log-probs', action='store_true',
        help=hidden_help('Write per-read alt log probabilities out in ' +
                         'non-standard VCF field.'))

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
        '--mod-all-paths', action='store_true',
        help=hidden_help('Compute forwards algorithm all paths score for ' +
                         'modified base calls. (Default: Viterbi ' +
                         'best-path score)'))
    mod_grp.add_argument(
        '--mod-aggregate-method', choices=list(mods.AGG_METHOD_NAMES),
        default=mods.BIN_THRESH_NAME,
        help=hidden_help('Modified base aggregation method. ' +
                         'Default: %(default)s'))
    mod_grp.add_argument(
        '--mod-binary-threshold', type=float, nargs=1,
        default=mods.DEFAULT_BINARY_THRESH,
        help=hidden_help('Threshold for modified base aggregation ' +
                         '(probability of modified/canonical base). ' +
                         'Only applicable for "--mod-aggregate-method ' +
                         'binary_threshold". Default: %(default)s'))
    mod_grp.add_argument(
        '--mod-calibration-filename',
        help=hidden_help('File containing emperical calibration for ' +
                         'modified base scores. As created by ' +
                         'megalodon/scripts/calibrate_mod_llr_scores.py. ' +
                         'Default: Load default calibration file.'))
    mod_grp.add_argument(
        '--mod-context-bases', type=int, default=mh.DEFAULT_MOD_CONTEXT,
        help=hidden_help('Context bases for modified base calling. ' +
                         'Default: %(default)d'))
    mod_grp.add_argument(
        '--mod-output-formats', nargs='+',
        default=[mh.MOD_BEDMETHYL_NAME,],
        choices=tuple(mh.MOD_OUTPUT_FMTS.keys()),
        help=hidden_help('Modified base aggregated output format(s). ' +
                         'Default: %(default)s'))
    mod_grp.add_argument(
        '--mod-positions-on-disk', action='store_true',
        help=hidden_help('Force modified base positions to be stored ' +
                         'only within on disk database table. This option ' +
                         'will reduce the RAM memory requirement, but may ' +
                         'drastically slow processing. Default: Store ' +
                         'positions in memory and on disk.'))
    mod_grp.add_argument(
        '--write-mod-log-probs', action='store_true',
        help=hidden_help('Write per-read modified base log probabilities ' +
                         'out in non-standard modVCF field.'))

    refout_grp = parser.add_argument_group('Reference Output Arguments')
    refout_grp.add_argument(
        '--output-per-read-references', action='store_true',
        help=hidden_help('Output per-read references.'))
    refout_grp.add_argument(
        '--refs-include-mods', action='store_true',
        help=hidden_help('Include modified base calls in per-read ' +
                         'reference output.'))
    refout_grp.add_argument(
        '--refs-include-variants', action='store_true',
        help=hidden_help('Include variant calls in per-read reference output.'))
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

    sigmap_grp = parser.add_argument_group('Signal Mapping Output Arguments')
    sigmap_grp.add_argument(
        '--output-signal-mappings', action='store_true',
        help=hidden_help('Output signal mapped file (see taiyaki).'))
    sigmap_grp.add_argument(
        '--signal-map-include-mods', action='store_true',
        help=hidden_help('Include modified base calls in signal ' +
                         'mapping output.'))
    sigmap_grp.add_argument(
        '--signal-map-percent-identity-threshold', type=float,
        help=hidden_help('Only include reads with higher percent identity ' +
                         'in signal mapping output.'))
    sigmap_grp.add_argument(
        '--signal-map-percent-coverage-threshold', type=float,
        help=hidden_help('Only include reads with higher read alignment ' +
                         'coverage in signal mapping output.'))
    sigmap_grp.add_argument(
        '--signal-map-length-range', type=int, nargs=2,
        help=hidden_help('Only include reads with specified read length ' +
                         'in signal mapping output.'))

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
        '--edge-buffer', type=int, default=mh.DEFAULT_EDGE_BUFFER,
        help=hidden_help('Do not process sequence variant or modified base ' +
                         'calls near edge of read mapping. ' +
                         'Default: %(default)d'))
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
        logger.warning('Running profiling. This may slow processing.')

    args, pr_ref_filts = parse_pr_ref_output(args)

    tai_model_fn = mh.get_model_fn(
        args.taiyaki_model_filename, args.load_default_model)
    model_info = backends.ModelInfo(
        args.processes, fast5s_dir=args.fast5s_dir,
        taiyaki_model_fn=tai_model_fn, devices=args.devices,
        chunk_size=args.chunk_size, chunk_overlap=args.chunk_overlap,
        max_concur_chunks=args.max_concurrent_chunks)
    args, mods_info = mods_validation(args, model_info)
    aligner = aligner_validation(args)
    args, vars_data = vars_validation(
        args, model_info.is_cat_mod, model_info.output_size, aligner)
    args, sig_map_info = parse_sig_map_output(args, model_info)

    process_all_reads(
        args.fast5s_dir, not args.not_recursive, args.num_reads,
        args.read_ids_filename, model_info, args.outputs,
        args.output_directory, args.basecalls_format, aligner, vars_data,
        args.processes, args.verbose_read_progress, args.suppress_progress,
        mods_info, args.database_safety, pr_ref_filts, sig_map_info)

    if aligner is not None:
        ref_fn = aligner.ref_fn
        map_out_fmt = aligner.out_fmt
        del aligner

    if mh.MAP_NAME in args.outputs:
        logger.info('Spawning process to sort mappings')
        map_p = post_process_mapping(
            args.output_directory, map_out_fmt, ref_fn)

    if mh.WHATSHAP_MAP_NAME in args.outputs:
        logger.info('Spawning process to sort whatshap mappings')
        whatshap_sort_fn, whatshap_p = post_process_whatshap(
            args.output_directory, map_out_fmt, ref_fn)

    if mh.VAR_NAME in args.outputs or mh.MOD_NAME in args.outputs:
        post_process_aggregate(
            mods_info, args.outputs, args.output_directory, args.processes,
            args.write_vcf_log_probs, args.heterozygous_factors, vars_data,
            args.write_mod_log_probs, args.suppress_progress)

    if mh.VAR_NAME in args.outputs:
        logger.info('Sorting output variant file')
        variant_fn = mh.get_megalodon_fn(args.output_directory, mh.VAR_NAME)
        sort_variant_fn = mh.add_fn_suffix(variant_fn, 'sorted')
        variants.sort_variants(variant_fn, sort_variant_fn)
        logger.info('Indexing output variant file')
        index_variant_fn = variants.index_variants(sort_variant_fn)

    if mh.WHATSHAP_MAP_NAME in args.outputs:
        if whatshap_p.is_alive():
            logger.info('Waiting for whatshap mappings sort')
            while whatshap_p.is_alive():
                sleep(0.001)
        logger.info(variants.get_whatshap_command(
            index_variant_fn, whatshap_sort_fn,
            mh.add_fn_suffix(variant_fn, 'phased')))

    if mh.MAP_NAME in args.outputs:
        if map_p.is_alive():
            logger.info('Waiting for mappings sort')
            while map_p.is_alive():
                sleep(0.001)

    return


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
