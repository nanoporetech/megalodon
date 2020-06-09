import os
import sys
import h5py
import queue
import threading
import traceback
from time import sleep
import multiprocessing as mp
from itertools import product
from collections import defaultdict, OrderedDict

import numpy as np
from tqdm import tqdm
from tqdm._utils import _term_move_up

from megalodon import (
    aggregate, backends, fast5_io, logging, mapping, mods,
    variants, megalodon_helper as mh)


# set blas library environment variables (without these the cblas calls
# can completely halt processing)
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'

_DO_PROFILE = False
_UNEXPECTED_ERROR_CODE = 'Unexpected error'
_UNEXPECTED_ERROR_FN = 'unexpected_megalodon_errors.{}.err'
_MAX_NUM_UNEXP_ERRORS = 50
DO_INTERPOLATE_SIG_POS = False
# update error text when 10% more errors are found
ERR_UPDATE_PROP = 0.1
LOGGER = logging.get_logger()


###################
# Read Processing #
###################

def handle_errors(func, args, r_vals, out_q, fast5_fn, failed_reads_q):
    try:
        out_q.put((func(*args), r_vals))
    except KeyboardInterrupt:
        failed_reads_q.put(
            (True, False, 'Keyboard interrupt', fast5_fn, None, 0))
        return
    except mh.MegaError as e:
        failed_reads_q.put((True, False, str(e), fast5_fn, None, 0))
    except Exception:
        failed_reads_q.put((
            True, False, _UNEXPECTED_ERROR_CODE, fast5_fn,
            traceback.format_exc(), 0))


def interpolate_sig_pos(r_to_q_poss, mapped_rl_cumsum):
    # TODO Need to test and optimize this function
    # interpolate signal positions for consecutive reference bases assigned
    # to the same query base
    ref_to_block = np.empty(r_to_q_poss.shape[0], dtype=np.int32)
    prev_query_pos = -1
    curr_stay_bases = 0
    for ref_pos, query_pos in enumerate(r_to_q_poss):
        # backsteps shouldn't be possible, but handled here
        if query_pos <= prev_query_pos:
            curr_stay_bases += 1
            continue
        ref_to_block[ref_pos - curr_stay_bases:ref_pos + 1] = np.around(
            np.linspace(
                start=ref_to_block[ref_pos - curr_stay_bases - 1],
                stop=mapped_rl_cumsum[query_pos], num=curr_stay_bases + 2,
                endpoint=True)[1:]).astype(np.int32)
        curr_stay_bases = 0
        prev_query_pos = query_pos
    # for stay at end of read there is no signal point to interpolate to
    # so copy last value
    if curr_stay_bases > 0:
        ref_to_block[
            ref_to_block.shape[0] - curr_stay_bases:] = ref_to_block[
                ref_to_block.shape[0] - curr_stay_bases - 1]
    return ref_to_block


def process_read(
        sig_info, model_info, bc_q, caller_conn, sig_map_q, ref_out_info,
        vars_data, vars_q, mods_q, mod_pos_conn, mods_info, failed_reads_q,
        signal_reversed):
    """ Workhorse per-read megalodon function (connects all the parts)
    """
    # perform basecalling using loaded backend
    (r_seq, r_qual, rl_cumsum, can_post, sig_info, post_w_mods,
     mods_scores) = model_info.basecall_read(
         sig_info, return_post_w_mods=mods_q is not None,
         return_mod_scores=mods_info.do_output_mods,
         update_sig_info=sig_map_q is not None,
         signal_reversed=signal_reversed)
    if bc_q is not None:
        if signal_reversed:
            if mods_scores is not None:
                mods_scores = mods_scores[::-1]
            bc_q.put((
                sig_info.read_id, r_seq[::-1], r_qual[::-1], mods_scores))
        else:
            bc_q.put((sig_info.read_id, r_seq, r_qual, mods_scores))

    # if no mapping connection return after basecalls are passed out
    if caller_conn is None:
        return

    # map read and record mapping from reference to query positions
    r_ref_seq, r_to_q_poss, r_ref_pos, r_cigar = mapping.map_read(
        r_seq, sig_info.read_id, caller_conn, signal_reversed)
    np_ref_seq = mh.seq_to_int(r_ref_seq, error_on_invalid=False)

    sig_map_res = None
    if sig_map_q is not None:
        pass_sig_map_filts = mapping.read_passes_filters(
            ref_out_info, len(r_seq), r_ref_pos.q_trim_start,
            r_ref_pos.q_trim_end, r_cigar)
        sig_map_res = signal_mapping.SIG_MAP_RESULT(
            pass_sig_map_filts, sig_info.fast5_fn, sig_info.dacs,
            sig_info.scale_params, r_ref_seq, sig_info.stride,
            sig_info.read_id, r_to_q_poss, rl_cumsum, r_ref_pos, ref_out_info)
        if not ref_out_info.annotate_mods and pass_sig_map_filts:
            try:
                sig_map_q.put(signal_mapping.get_remapping(*sig_map_res[1:]))
            except Exception as e:
                LOGGER.debug((
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
    if DO_INTERPOLATE_SIG_POS:
        ref_to_block = interpolate_sig_pos(r_to_q_poss, mapped_rl_cumsum)
    else:
        ref_to_block = mapped_rl_cumsum[r_to_q_poss]

    if vars_q is not None:
        assert not signal_reversed, (
            'Reversed raw signal (RNA) not compatible with sequence ' +
            'variant detection.')
        mapped_can_post = can_post[post_mapped_start:post_mapped_end]
        handle_errors(
            func=variants.call_read_vars,
            args=(vars_data, r_ref_pos, np_ref_seq, ref_to_block,
                  mapped_can_post),
            r_vals=(sig_info.read_id, r_ref_pos.chrm, r_ref_pos.strand,
                    r_ref_pos.start, r_ref_seq, len(r_seq),
                    r_ref_pos.q_trim_start, r_ref_pos.q_trim_end, r_cigar),
            out_q=vars_q,
            fast5_fn=sig_info.fast5_fn + ':::' + sig_info.read_id,
            failed_reads_q=failed_reads_q)
    if mods_q is not None:
        mapped_post_w_mods = post_w_mods[post_mapped_start:post_mapped_end]
        mod_sig_map_q = sig_map_q if ref_out_info.annotate_mods else None
        handle_errors(
            func=mods.call_read_mods,
            args=(r_ref_pos, r_ref_seq, ref_to_block, mapped_post_w_mods,
                  mods_info, mod_pos_conn, mod_sig_map_q, sig_map_res,
                  signal_reversed, sig_info.read_id),
            r_vals=(sig_info.read_id, r_ref_pos.chrm, r_ref_pos.strand,
                    r_ref_pos.start, r_ref_seq, len(r_seq),
                    r_ref_pos.q_trim_start, r_ref_pos.q_trim_end, r_cigar),
            out_q=mods_q,
            fast5_fn=sig_info.fast5_fn + ':::' + sig_info.read_id,
            failed_reads_q=failed_reads_q)


####################
# Multi-processing #
####################

def _get_bc_queue(
        bc_q, bc_conn, out_dir, bc_fmt, do_output_mods, mod_long_names):
    def write_read(read_id, r_seq, r_qual, mods_scores):
        if write_fastq:
            if r_qual is None:
                r_qual = '!' * len(r_seq)
            bc_fp.write('@{}\n{}\n+\n{}\n'.format(read_id, r_seq, r_qual))
        else:
            bc_fp.write('>{}\n{}\n'.format(read_id, r_seq))

        if do_output_mods:
            try:
                mods_fp.create_dataset(
                    'Reads/' + read_id, data=mods_scores,
                    compression="gzip")
            except RuntimeError:
                # same read_id encountered previously
                pass

    bc_fp = open(mh.get_megalodon_fn(out_dir, mh.BC_NAME) + '.' + bc_fmt, 'w')
    write_fastq = bc_fmt == 'fastq'
    # TODO convert this to writing un-aligned sam with htseq recommended format
    if do_output_mods:
        mods_fp = h5py.File(mh.get_megalodon_fn(out_dir, mh.BC_MODS_NAME), 'w')
        mods_fp.create_group('Reads')
        mods_fp.create_dataset(
            'mod_long_names', data=np.array(mod_long_names, dtype='S'),
            dtype=h5py.special_dtype(vlen=str))

    while True:
        try:
            read_id, r_seq, r_qual, mods_scores = bc_q.get(block=False)
            write_read(read_id, r_seq, r_qual, mods_scores)
        except queue.Empty:
            if bc_conn.poll():
                break
            sleep(0.001)
            continue

    while not bc_q.empty():
        read_id, r_seq, r_qual, mods_scores = bc_q.get(block=False)
        write_read(read_id, r_seq, r_qual, mods_scores)

    bc_fp.close()
    if do_output_mods:
        mods_fp.close()


def _process_reads_worker(
        read_file_q, bc_q, vars_q, failed_reads_q, mods_q, mod_pos_conn,
        caller_conn, sig_map_q, ref_out_info, model_info, vars_data, mods_info,
        device, signal_reversed):
    # wrap process prep in try loop to avoid stalled command
    try:
        model_info.prep_model_worker(device)
        vars_data.reopen_variant_index()
        LOGGER.debug('Starting read worker {}'.format(mp.current_process()))
        sig_info = None
    except Exception:
        if caller_conn is not None:
            caller_conn.send(True)
        LOGGER.debug(('Read worker {} has failed process preparation.\n' +
                      'Full error traceback:\n{}').format(
                          mp.current_process(), traceback.format_exc()))

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
                LOGGER.debug('Gracefully exiting read worker {}'.format(
                    mp.current_process()))
                break
            LOGGER.debug('Analyzing read {}'.format(read_id))
            sig_info = model_info.extract_signal_info(
                fast5_fn, read_id, sig_map_q is not None)
            process_read(
                sig_info, model_info, bc_q, caller_conn, sig_map_q,
                ref_out_info, vars_data, vars_q, mods_q, mod_pos_conn,
                mods_info, failed_reads_q, signal_reversed)
            failed_reads_q.put((
                False, True, None, None, None, sig_info.raw_len))
            LOGGER.debug('Successfully processed read {}'.format(read_id))
        except KeyboardInterrupt:
            failed_reads_q.put((
                True, True, 'Keyboard interrupt', fast5_fn, None, 0))
            LOGGER.debug('Keyboard interrupt during read {}'.format(read_id))
            return
        except mh.MegaError as e:
            raw_len = sig_info.raw_len if hasattr(sig_info, 'raw_len') else 0
            failed_reads_q.put((
                True, True, str(e), fast5_fn + ':::' + read_id, None, raw_len))
            LOGGER.debug('Incomplete processing for read {} ::: {}'.format(
                read_id, str(e)))
        except Exception:
            failed_reads_q.put((
                True, True, _UNEXPECTED_ERROR_CODE, fast5_fn + ':::' + read_id,
                traceback.format_exc(), 0))
            LOGGER.debug('Unexpected error for read {}'.format(read_id))

    if mod_pos_conn is not None:
        mod_pos_conn.close()


if _DO_PROFILE:
    _process_reads_wrapper = _process_reads_worker

    def _process_reads_worker(*args):
        import cProfile
        cProfile.runctx('_process_reads_wrapper(*args)', globals(), locals(),
                        filename='read_processing.prof')


############################
# Post Per-read Processing #
############################

def post_process_mapping(map_bn, map_fmt, ref_fn):
    map_fn = map_bn + '.' + map_fmt
    map_sort_fn = map_bn + '.sorted.bam'
    map_p = mp.Process(
        target=mapping.sort_and_index_mapping,
        args=(map_fn, map_sort_fn, ref_fn, True), daemon=True)
    map_p.start()
    sleep(0.001)

    return map_p, map_sort_fn


def start_sort_mapping_procs(outputs, out_dir, map_out_fmt, ref_fn, mlns):
    map_p = mod_map_ps = var_map_p = var_sort_fn = None
    if mh.MAP_NAME in outputs:
        LOGGER.info('Spawning process to sort mappings')
        map_p, _ = post_process_mapping(
            mh.get_megalodon_fn(out_dir, mh.MAP_NAME), map_out_fmt, ref_fn)
    if mh.MOD_MAP_NAME in outputs:
        LOGGER.info('Spawning process to sort modified base mappings')
        mod_map_ps = [post_process_mapping(
            '{}.{}'.format(mh.get_megalodon_fn(out_dir, mh.MOD_MAP_NAME), mln),
            map_out_fmt, ref_fn)[0] for _, mln in mlns]
    if mh.VAR_MAP_NAME in outputs:
        LOGGER.info('Spawning process to sort variant mappings')
        var_map_p, var_sort_fn = post_process_mapping(
            mh.get_megalodon_fn(out_dir, mh.VAR_MAP_NAME), map_out_fmt, ref_fn)
    return map_p, mod_map_ps, var_map_p, var_sort_fn


def get_map_procs(
        map_p, mod_map_ps, var_map_p, var_sort_fn, index_variant_fn,
        variant_fn):
    if var_map_p is not None:
        if var_map_p.is_alive():
            LOGGER.info('Waiting for variant mappings sort')
            while var_map_p.is_alive():
                sleep(0.001)
        if index_variant_fn is not None and var_sort_fn is not None:
            LOGGER.info(variants.get_whatshap_command(
                index_variant_fn, var_sort_fn,
                mh.add_fn_suffix(variant_fn, 'phased')))
    if mod_map_ps is not None:
        if any(mod_map_p.is_alive() for mod_map_p in mod_map_ps):
            LOGGER.info('Waiting for modified base mappings sort')
            while any(mod_map_p.is_alive() for mod_map_p in mod_map_ps):
                sleep(0.001)
    if map_p is not None:
        if map_p.is_alive():
            LOGGER.info('Waiting for mappings sort')
            while map_p.is_alive():
                sleep(0.001)


##########################
# Dynamic error updating #
##########################

def _fill_files_queue(
        read_file_q, fast5s_dir, num_reads, read_ids_fn, recursive, num_ps,
        num_reads_conn):
    valid_read_ids = None
    if read_ids_fn is not None:
        with open(read_ids_fn) as read_ids_fp:
            valid_read_ids = set(line.strip() for line in read_ids_fp)
    used_read_ids = set()
    # fill queue with read filename and read id tuples
    for fast5_fn, read_id in fast5_io.iterate_fast5_reads(
            fast5s_dir, recursive=recursive):
        if valid_read_ids is not None and read_id not in valid_read_ids:
            continue
        if read_id in used_read_ids:
            LOGGER.debug(
                ('Read ID ({}) found in previous read and will not ' +
                 'process from {}.').format(read_id, fast5_fn))
            continue
        if fast5_fn is None or read_id is None:
            continue
        read_file_q.put((fast5_fn, read_id))
        used_read_ids.add(read_id)
        if num_reads is not None and len(used_read_ids) >= num_reads:
            break
    # add None to indicate that read processes should return
    for _ in range(num_ps):
        read_file_q.put((None, None))
    num_reads_conn.send(len(used_read_ids))


def format_fail_summ(header, fail_summ=[], reads_called=0, num_errs=None):
    summ_errs = sorted(fail_summ)[::-1]
    if num_errs is not None:
        summ_errs = summ_errs[:num_errs]
        if len(summ_errs) < num_errs:
            summ_errs.extend([(None, '') for _ in range(
                num_errs - len(summ_errs))])
    errs_str = '\n'.join(
        "{:8.1f}% ({:>7} reads)".format(
            100 * n_fns / float(reads_called), n_fns) +
        " : " + '{:<80}'.format(err)
        if (n_fns is not None and reads_called > 0) else
        '     -----' for n_fns, err in summ_errs)
    return '\n'.join((header, errs_str))


def prep_errors_bar(
        num_update_errors, tot_reads, suppress_progress, do_show_qs, getter_qs,
        curr_num_reads=0, start_time=None):
    num_qs = 0
    if do_show_qs:
        valid_q_names = [q_name for q_name, q_vals in getter_qs.items()
                         if q_vals.queue is not None]
        num_qs = len(valid_q_names)
    if num_update_errors > 0 and not suppress_progress:
        # add lines for dynamic error messages
        # note 2 extra lines for header and bar
        sys.stderr.write(
            '\n'.join(['' for _ in range(num_update_errors + 2)]))
    bar = prog_prefix = bar_header = q_bars = None
    if suppress_progress:
        num_update_errors = 0
    else:
        bar = tqdm(total=tot_reads, smoothing=0, initial=curr_num_reads,
                   unit=' read(s)', dynamic_ncols=True, position=0,
                   desc='Read Processing')
        if start_time is not None:
            bar.start_t = start_time
        if num_qs > 0:
            q_bars = OrderedDict((q_name, tqdm(
                desc=q_name, total=mh._MAX_QUEUE_SIZE, smoothing=0,
                dynamic_ncols=True, position=q_num + 1,
                bar_format='output queue capacity {desc: <20}: ' +
                '{percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt}'))
                                 for q_num, q_name in enumerate(valid_q_names))
    if num_update_errors > 0:
        prog_prefix = ''.join(
            [_term_move_up(), ] * (num_update_errors + 1)) + '\r'
        if num_qs > 0:
            bar_header = (
                '{} most common unsuccessful processing stages (full ' +
                'output queues indicate I/O bottleneck):').format(
                    num_update_errors)
        else:
            bar_header = (
                '{} most common unsuccessful processing stages:').format(
                    num_update_errors)
        # write failed read update header
        bar.write(prog_prefix + format_fail_summ(
            bar_header, num_errs=num_update_errors), file=sys.stderr)

    return bar, q_bars, prog_prefix, bar_header


def _get_fail_queue(
        failed_reads_q, f_conn, getter_num_reads_conn, num_update_errors,
        suppress_progress, do_show_qs, getter_qs):
    def update_prog(reads_called, sig_called, unexp_err_fp, last_err_write,
                    read_called=True):
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
                        'ksamp/s': (sig_called / 1000) /
                        bar.format_dict['elapsed']}, refresh=False)
                except AttributeError:
                    # sometimes get no format_dict error
                    # so don't include ksample/s if so
                    pass
                if q_bars is not None:
                    for q_name, q_bar in q_bars.items():
                        q_bar.n = getter_qs[q_name].queue.qsize()
                        # trigger display refresh
                        q_bar.update(0)
                if read_called:
                    bar.update()
                if num_update_errors > 0:
                    err_types = [
                        (len(fns), err) for err, fns in failed_reads.items()]
                    num_errs = sum((x[0] for x in err_types))
                    if num_errs > 0 and (
                            last_err_write == 0 or
                            num_errs / last_err_write > 1 + ERR_UPDATE_PROP):
                        last_err_write = num_errs
                        bar.write(prog_prefix + format_fail_summ(
                            bar_header, err_types, reads_called,
                            num_update_errors), file=sys.stderr)

        return unexp_err_fp, last_err_write

    LOGGER.info('Processing reads.')
    reads_called = sig_called = last_err_write = 0
    unexp_err_fp = None
    failed_reads = defaultdict(list)
    bar, q_bars, prog_prefix, bar_header = prep_errors_bar(
        num_update_errors, None, suppress_progress, do_show_qs, getter_qs)
    while True:
        try:
            try:
                (is_err, do_update_prog, err_type, fast5_fn,
                 err_tb, n_sig) = failed_reads_q.get(block=False)
                sig_called += n_sig
                reads_called += 1
                unexp_err_fp, last_err_write = update_prog(
                    reads_called, sig_called, unexp_err_fp, last_err_write)
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
    if not suppress_progress:
        if q_bars is not None:
            while any(not getter_qs[q_name].queue.empty()
                      for q_name in q_bars.keys()):
                unexp_err_fp, last_err_write = update_prog(
                    reads_called, 0, unexp_err_fp, last_err_write, False)
        bar.close()
        if q_bars is not None:
            for q_bar in q_bars.values():
                q_bar.close()

    if len(failed_reads[_UNEXPECTED_ERROR_CODE]) >= 1:
        LOGGER.warning((
            'Unexpected errors occured. See full ' +
            'error stack traces for first (up to) {0:d} errors in ' +
            '"{1}"').format(_MAX_NUM_UNEXP_ERRORS, unexp_err_fp.name))
    if any(len(fns) > 0 for fns in failed_reads.values()):
        LOGGER.info(
            format_fail_summ(
                'Unsuccessful processing types:',
                [(len(fns), err) for err, fns in failed_reads.items()
                 if len(fns) > 0], reads_called))
        # TODO flag to output failed read names to file
    else:
        LOGGER.info('All reads processed successfully.')


#######################
# All read processing #
#######################

def process_all_reads(
        fast5s_dir, recursive, num_reads, read_ids_fn, model_info, outputs,
        out_dir, bc_fmt, aligner, vars_data, num_ps, num_update_errors,
        suppress_progress, mods_info, db_safety, ref_out_info, signal_reversed,
        do_show_qs):
    LOGGER.info('Preparing workers to process reads.')
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

    # start output type getters/writers
    getter_qs = OrderedDict(
        (out_name, mh.GETTER_PROC(None, None, None)) for out_name in (
            mh.BC_NAME, mh.MAP_NAME, mh.SIG_MAP_NAME, mh.PR_VAR_NAME,
            mh.PR_MOD_NAME))
    if mh.BC_NAME in outputs or mh.BC_MODS_NAME in outputs:
        if mh.BC_NAME not in outputs:
            outputs.append(mh.BC_NAME)
        getter_qs[mh.BC_NAME] = mh.create_getter_q(
            _get_bc_queue, (out_dir, bc_fmt, mods_info.do_output_mods,
                            mods_info.mod_long_names))
    if mh.MAP_NAME in outputs:
        do_output_pr_refs = (mh.PR_REF_NAME in outputs and
                             not ref_out_info.annotate_mods and
                             not ref_out_info.annotate_vars)
        getter_qs[mh.MAP_NAME] = mh.create_getter_q(
            mapping._get_map_queue, (
                out_dir, aligner.ref_names_and_lens, aligner.out_fmt,
                aligner.ref_fn, do_output_pr_refs, ref_out_info))
    if mh.PR_VAR_NAME in outputs:
        pr_refs_fn = mh.get_megalodon_fn(out_dir, mh.PR_REF_NAME) if (
            mh.PR_REF_NAME in outputs and ref_out_info.annotate_vars) else None
        var_map_fn = (
            mh.get_megalodon_fn(out_dir, mh.VAR_MAP_NAME) + '.' +
            aligner.out_fmt) if mh.VAR_MAP_NAME in outputs else None
        vars_txt_fn = (mh.get_megalodon_fn(out_dir, mh.PR_VAR_TXT_NAME)
                       if vars_data.write_vars_txt else None)
        getter_qs[mh.PR_VAR_NAME] = mh.create_getter_q(
            variants._get_variants_queue, (
                mh.get_megalodon_fn(out_dir, mh.PR_VAR_NAME),
                vars_txt_fn, db_safety, pr_refs_fn, ref_out_info,
                var_map_fn, aligner.ref_names_and_lens, aligner.ref_fn,
                vars_data.loc_index_in_memory))
    if mh.PR_MOD_NAME in outputs:
        pr_refs_fn = mh.get_megalodon_fn(out_dir, mh.PR_REF_NAME) if (
            mh.PR_REF_NAME in outputs and ref_out_info.annotate_mods) else None
        mod_map_fns = None
        if mh.MOD_MAP_NAME in outputs:
            mod_map_fns = [(mod_base, '{}.{}.'.format(
                mh.get_megalodon_fn(out_dir, mh.MOD_MAP_NAME), mln))
                           for mod_base, mln in mods_info.mod_long_names]
        mods_txt_fn = (mh.get_megalodon_fn(out_dir, mh.PR_MOD_TXT_NAME)
                       if mods_info.write_mods_txt else None)
        mods_db_fn = mh.get_megalodon_fn(out_dir, mh.PR_MOD_NAME)
        mods.init_mods_db(
            mods_db_fn, db_safety, aligner.ref_names_and_lens, mods_info)
        # TODO handle on-disk index creation
        getter_qs[mh.PR_MOD_NAME] = mh.create_getter_q(
            mods._get_mods_queue, (
                mods_db_fn, db_safety, aligner.ref_names_and_lens,
                aligner.ref_fn, mods_txt_fn, pr_refs_fn, ref_out_info,
                mod_map_fns, aligner.out_fmt, mods_info))
    if mh.SIG_MAP_NAME in outputs:
        alphabet_info = signal_mapping.get_alphabet_info(
            ref_out_info.alphabet, ref_out_info.collapse_alphabet,
            ref_out_info.mod_long_names)
        sig_map_fn = mh.get_megalodon_fn(out_dir, mh.SIG_MAP_NAME)
        getter_qs[mh.SIG_MAP_NAME] = mh.create_getter_q(
            signal_mapping.write_signal_mappings,
            (sig_map_fn, alphabet_info))
    # progress and failed reads getter (no limit on failed reads queue
    # in case error occurs there, don't halt run
    fr_prog_getter = mh.create_getter_q(
        _get_fail_queue, (getter_num_reads_conn, num_update_errors,
                          suppress_progress, do_show_qs, getter_qs),
        max_size=None)

    proc_reads_ps, map_conns, mod_pos_conns = [], [], []
    for device in model_info.process_devices:
        if aligner is None:
            map_conn, caller_conn = None, None
        else:
            map_conn, caller_conn = mp.Pipe()
        map_conns.append(map_conn)
        if mh.PR_MOD_NAME in outputs:
            mod_pos_conn, mod_pos_db_conn = mp.Pipe()
        else:
            mod_pos_conn, mod_pos_db_conn = None, None
        mod_pos_conns.append(mod_pos_db_conn)
        p = mp.Process(
            target=_process_reads_worker, args=(
                read_file_q, getter_qs[mh.BC_NAME].queue,
                getter_qs[mh.PR_VAR_NAME].queue, fr_prog_getter.queue,
                getter_qs[mh.PR_MOD_NAME].queue, mod_pos_conn, caller_conn,
                getter_qs[mh.SIG_MAP_NAME].queue, ref_out_info, model_info,
                vars_data, mods_info, device, signal_reversed))
        p.daemon = True
        p.start()
        proc_reads_ps.append(p)
        mod_pos_conn.close()

    # extract and enter modified base position and mod_base database ids in
    # separate processes in order to get around bottleneck in
    # single threaded mod base database data entry
    mod_pos_p = None
    if mh.PR_MOD_NAME in outputs:
        mod_pos_p = mp.Process(
            target=mods._mod_aux_table_inserts, args=(
                mods_db_fn, db_safety, mods_info.pos_index_in_memory,
                mod_pos_conns), daemon=True)
        mod_pos_p.start()

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
                args=(aligner, map_conn, getter_qs[mh.MAP_NAME].queue))
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
        if fr_prog_getter.proc.is_alive():
            fr_prog_getter.conn.send(True)
            fr_prog_getter.proc.join()
        mod_pos_p.join()
        for out_name, getter_q in getter_qs.items():
            if out_name in outputs and getter_q.proc.is_alive():
                getter_q.conn.send(True)
                if out_name == mh.PR_VAR_NAME:
                    LOGGER.info(
                        'Waiting for variants database to complete indexing.')
                elif out_name == mh.PR_MOD_NAME:
                    LOGGER.info(
                        'Waiting for mods database to complete indexing.')
                getter_q.proc.join()
    except KeyboardInterrupt:
        LOGGER.error('Exiting due to keyboard interrupt.')
        sys.exit(1)


####################
# Input validation #
####################

def parse_aligner_args(args):
    if len(mh.ALIGN_OUTPUTS.intersection(args.outputs)) > 0:
        if args.reference is None:
            LOGGER.error(
                ('Output(s) requiring reference alignment requested ({}), ' +
                 'but --reference not provided.').format(', '.join(
                    mh.ALIGN_OUTPUTS.intersection(args.outputs))))
            sys.exit(1)
        LOGGER.info('Loading reference.')
        if not (os.path.exists(args.reference) and
                os.path.isfile(args.reference)):
            LOGGER.error('Provided reference file does not exist or is ' +
                         'not a file.')
            sys.exit(1)
        aligner = mapping.alignerPlus(
            str(args.reference), preset=str('map-ont'), best_n=1)
        setattr(aligner, 'out_fmt', args.mappings_format)
        setattr(aligner, 'ref_fn', mh.resolve_path(args.reference))
        aligner.add_ref_lens()
        if mh.MAP_NAME in args.outputs:
            # test that alignment file can be opened
            mapping.test_open_alignment_out_file(
                args.output_directory, aligner.out_fmt,
                aligner.ref_names_and_lens, aligner.ref_fn)
    else:
        aligner = None
        if args.reference is not None:
            LOGGER.warning(
                '[--reference] provided, but no [--outputs] requiring ' +
                'alignment was requested. Argument will be ignored.')
    return aligner


def parse_var_args(args, model_info, aligner):
    if args.ref_include_variants and mh.PR_VAR_NAME not in args.outputs:
        LOGGER.warning('--ref-include-variants set, so adding ' +
                       '"per_read_vars" to --outputs.')
        args.outputs.append(mh.PR_VAR_NAME)
    if mh.VAR_MAP_NAME in args.outputs and \
       mh.PR_VAR_NAME not in args.outputs:
        LOGGER.warning((
            'Adding "{}" to --outputs since "{}" was requested. For full ' +
            'phased variant pipeline add "{}" or run aggregation after run ' +
            'is complete.').format(mh.PR_VAR_NAME, mh.VAR_MAP_NAME,
                                   mh.VAR_NAME))
        args.outputs.append(mh.PR_VAR_NAME)
    if mh.VAR_NAME in args.outputs and mh.PR_VAR_NAME not in args.outputs:
        LOGGER.warning((
            'Adding "{}" to --outputs since "{}" was requested.').format(
                mh.PR_VAR_NAME, mh.VAR_NAME))
        args.outputs.append(mh.PR_VAR_NAME)
    if mh.PR_VAR_NAME in args.outputs and args.variant_filename is None:
        LOGGER.error(
            '{} output requested, '.format(mh.PR_VAR_NAME) +
            'but --variant-filename not provided.')
        sys.exit(1)
    if mh.PR_VAR_NAME in args.outputs and not (
            model_info.is_cat_mod or
            mh.nstate_to_nbase(model_info.output_size) == 4):
        LOGGER.error(
            'Variant calling from naive modified base flip-flop model is ' +
            'not supported.')
        sys.exit(1)
    var_calib_fn = mh.get_var_calibration_fn(
        model_info.params.pyguppy.config, args.variant_calibration_filename,
        args.disable_variant_calibration) \
        if mh.PR_VAR_NAME in args.outputs else None
    try:
        vars_data = variants.VarData(
            args.variant_filename, args.max_indel_size,
            args.variant_all_paths, args.write_variants_text,
            args.variant_context_bases, var_calib_fn,
            variants.HAPLIOD_MODE if args.haploid else variants.DIPLOID_MODE,
            aligner, edge_buffer=args.edge_buffer,
            context_min_alt_prob=args.context_min_alt_prob,
            loc_index_in_memory=not args.variant_locations_on_disk,
            variants_are_atomized=args.variants_are_atomized)
    except mh.MegaError as e:
        LOGGER.error(str(e))
        sys.exit(1)
    if args.variant_filename is not None and \
       mh.PR_VAR_NAME not in args.outputs:
        LOGGER.warning(
            '--variants-filename provided, but variants output not ' +
            'requested (via --outputs). Argument will be ignored.')
    if args.rna and mh.PR_VAR_NAME in args.outputs:
        LOGGER.error(
            'Sequence variant analysis of RNA data is not currently ' +
            'supported.')
        sys.exit(1)
    return args, vars_data


def parse_mod_args(args, model_info):
    if args.ref_include_mods and args.ref_mods_all_motifs is not None:
        LOGGER.warning(
            '--ref-include-mods and --ref-mods-all-motifs are not ' +
            'compatible. Ignoring --ref-include-mods')
        args.ref_include_mods = False
    if args.ref_include_mods and not (mh.SIG_MAP_NAME in args.outputs or
                                      mh.PR_REF_NAME in args.outputs):
        LOGGER.warning((
            '--ref-include-mods specified, but neither {} or {} specified ' +
            'in outputs. Ignoring --ref-include-mods').format(
                mh.SIG_MAP_NAME, mh.PR_REF_NAME))
        args.ref_include_mods = False
    if mh.MOD_MAP_NAME in args.outputs and \
       mh.PR_MOD_NAME not in args.outputs:
        LOGGER.warning((
            'Adding "{}" to --outputs since "{}" was requested.').format(
                mh.PR_MOD_NAME, mh.MOD_MAP_NAME))
        args.outputs.append(mh.PR_MOD_NAME)
    if args.ref_include_mods and mh.PR_MOD_NAME not in args.outputs:
        LOGGER.warning('--ref-include-mods set, so adding ' +
                       '"per_read_mods" to --outputs.')
        args.outputs.append(mh.PR_MOD_NAME)
    if mh.PR_MOD_NAME not in args.outputs and mh.MOD_NAME in args.outputs:
        LOGGER.warning('"mods" output requested, so "per_read_mods" will ' +
                       'be added to outputs.')
        args.outputs.append(mh.PR_MOD_NAME)
    if mh.PR_MOD_NAME in args.outputs and not model_info.is_cat_mod:
        LOGGER.error((
            '{} output requested, but specified model does not support ' +
            'calling modified bases.').format(mh.PR_MOD_NAME))
        sys.exit(1)
    if model_info.is_cat_mod and mh.PR_MOD_NAME not in args.outputs and \
       mh.BC_MODS_NAME not in args.outputs:
        LOGGER.warning(
            ('Model supporting modified base calling specified, but neither ' +
             '{} nor {} requested.').format(mh.PR_MOD_NAME, mh.BC_MODS_NAME))
    if args.mod_motif is not None and mh.PR_MOD_NAME not in args.outputs:
        LOGGER.warning(('--mod-motif provided, but {} not requested. ' +
                        'Ignoring --mod-motif.').format(mh.PR_MOD_NAME))
        args.mod_motif = None

    mod_calib_fn = (mh.get_mod_calibration_fn(
        model_info.params.pyguppy.config, args.mod_calibration_filename,
        args.disable_mod_calibration)
                    if mh.PR_MOD_NAME in args.outputs else None)
    if args.mod_aggregate_method == mh.MOD_EM_NAME:
        agg_info = mods.AGG_INFO(mh.MOD_EM_NAME, None)
    elif args.mod_aggregate_method == mh.MOD_BIN_THRESH_NAME:
        agg_info = mods.AGG_INFO(
            mh.MOD_BIN_THRESH_NAME, args.mod_binary_threshold)
    mods_info = mods.ModInfo(
        model_info=model_info, all_mod_motifs_raw=args.mod_motif,
        mod_all_paths=args.mod_all_paths, write_mods_txt=args.write_mods_text,
        mod_context_bases=args.mod_context_bases,
        do_output_mods=mh.BC_MODS_NAME in args.outputs,
        mods_calib_fn=mod_calib_fn, mod_output_fmts=args.mod_output_formats,
        edge_buffer=args.edge_buffer,
        pos_index_in_memory=not args.mod_positions_on_disk, agg_info=agg_info,
        mod_thresh=args.ref_mod_threshold,
        do_ann_all_mods=args.ref_include_mods,
        do_ann_per_mod=mh.MOD_MAP_NAME in args.outputs,
        map_base_conv=args.mod_map_base_conv)
    return args, mods_info


def parse_ref_mods_all_motifs(ref_mod_motifs_raw, sm_alphabet_info):
    sm_alphabet = sm_alphabet_info.alphabet
    sm_coll_alphabet = sm_alphabet_info.collapse_alphabet
    sm_mlns = sm_alphabet_info.mod_long_names
    ref_mod_motifs_init = []
    for mod_base, mln, motif, rel_pos in ref_mod_motifs_raw:
        rel_pos = int(rel_pos)
        ref_base = motif[rel_pos]
        if mod_base not in sm_alphabet:
            # add mod base to alphabet
            sm_alphabet += mod_base
            sm_coll_alphabet += ref_base
            sm_mlns.append(mln)
        else:
            # check that mod base matches current model
            base_idx = sm_alphabet.index(mod_base)
            if sm_coll_alphabet[base_idx] != ref_base:
                raise mh.MegaError((
                    'Canonical base ({}) specified by --ref-mods-all-motifs ' +
                    'does not match model alphabet base ({}).').format(
                        ref_base, sm_coll_alphabet[base_idx]))
            mod_idx = sm_alphabet_info.mod_bases.index(mod_base)
            if sm_mlns[mod_idx] != mln:
                raise mh.MegaError((
                    'Modified base long name ({}) specified by ' +
                    '--ref-mods-all-motifs does not match model alphabet ' +
                    'long name ({}).').format(
                        mln, sm_mlns[mod_idx]))
        ref_mod_motifs_init.append((mod_base, mln, motif, rel_pos))
    new_sm_alphabet_info = signal_mapping.get_alphabet_info(
        sm_alphabet, sm_coll_alphabet, sm_mlns)
    ref_mod_motifs_w_ints = []
    for mod_base, mln, motif, rel_pos in ref_mod_motifs_init:
        int_mod_base = new_sm_alphabet_info.alphabet.index(mod_base)
        # iterate over all real sequences for the canonical motif
        for motif_bases in product(*(mh.SINGLE_LETTER_CODE[b] for b in motif)):
            int_motif = np.array([new_sm_alphabet_info.alphabet.index(base)
                                  for base in motif_bases])
            ref_mod_motifs_w_ints.append(
                (mod_base, int_mod_base, mln, int_motif, rel_pos))
    return ref_mod_motifs_w_ints, new_sm_alphabet_info


def parse_ref_out_args(args, model_info):
    output_pr_refs = output_sig_maps = False
    if mh.SIG_MAP_NAME in args.outputs or mh.PR_REF_NAME in args.outputs:
        if args.ref_include_variants and args.ref_include_mods:
            LOGGER.error('Cannot output both modified base and variants in ' +
                         'per-read references (remove one of ' +
                         '--refs-include-variants or --refs-include-mods).')
            sys.exit(1)
        if args.ref_include_mods and args.ref_mods_all_motifs is not None:
            LOGGER.warning(
                '--ref-include-mods and --ref-mods-all-motifs are not ' +
                'compatible. Ignoring --ref-include-mods')
            args.ref_include_mods = False

    sm_alphabet_info = None
    if mh.SIG_MAP_NAME in args.outputs:
        output_sig_maps = True
        # import here so that taiyaki is not required unless outputting
        # signal mappings
        from megalodon import signal_mapping
        global signal_mapping
        sm_alphabet_info = signal_mapping.get_alphabet_info_from_model(
            model_info)
        if args.ref_include_mods and mh.PR_MOD_NAME not in args.outputs:
            LOGGER.warning(('--ref-include-mods set, so adding ' +
                            '"{}" to --outputs.').format(mh.PR_MOD_NAME))
            args.outputs.append(mh.PR_MOD_NAME)

    if mh.PR_REF_NAME in args.outputs:
        output_pr_refs = True
        if args.ref_include_variants and mh.PR_VAR_NAME not in args.outputs:
            LOGGER.warning('--refs-include-variants set, so adding ' +
                           'per_read_variants to --outputs.')
            args.outputs.append(mh.PR_VAR_NAME)
    else:
        if args.ref_include_variants:
            LOGGER.warning(
                '{0} set but {1} not requested. Ignoring {0}.'.format(
                    '--ref-include-variants', mh.PR_REF_NAME))
            args.ref_include_variants = None

    if mh.SIG_MAP_NAME not in args.outputs and \
       mh.PR_REF_NAME not in args.outputs:
        for arg_val, arg_str in (
                (args.ref_mods_all_motifs, '--ref-mods-all-motifs'),
                (args.ref_length_range, '--ref-length-range'),
                (args.ref_percent_identity_threshold,
                 '--ref-percent-identity-threshold'),
                (args.ref_percent_coverage_threshold,
                 '--ref-percent-coverage-threshold')):
            if arg_val is not None:
                LOGGER.warning((
                    '{0} set but neither {1} nor {2} requested. Ignoring ' +
                    '{0}.').format(arg_str, mh.SIG_MAP_NAME, mh.PR_REF_NAME))
                arg_val = None
        for arg_flag, arg_str in (
                (args.ref_include_mods, '--ref-include-mods'), ):
            if arg_flag:
                LOGGER.warning((
                    '{0} set but neither {1} nor {2} requested. Ignoring ' +
                    '{0}.').format(arg_str, mh.SIG_MAP_NAME, mh.PR_REF_NAME))
                arg_flag = False

    min_len, max_len = (args.ref_length_range
                        if args.ref_length_range is not None else
                        (None, None))
    ref_mods_all_motifs = None
    if args.ref_mods_all_motifs is not None:
        ref_mods_all_motifs, sm_alphabet_info = parse_ref_mods_all_motifs(
            args.ref_mods_all_motifs, sm_alphabet_info)

    sm_alphabet = sm_coll_alphabet = sm_mlns = None
    if sm_alphabet_info is not None:
        sm_alphabet = sm_alphabet_info.alphabet
        sm_coll_alphabet = sm_alphabet_info.collapse_alphabet
        sm_mlns = sm_alphabet_info.mod_long_names
    ref_out_info = mh.REF_OUT_INFO(
        pct_idnt=args.ref_percent_identity_threshold,
        pct_cov=args.ref_percent_coverage_threshold,
        min_len=min_len, max_len=max_len, alphabet=sm_alphabet,
        collapse_alphabet=sm_coll_alphabet, mod_long_names=sm_mlns,
        annotate_mods=args.ref_include_mods,
        annotate_vars=args.ref_include_variants,
        output_sig_maps=output_sig_maps, output_pr_refs=output_pr_refs,
        ref_mods_all_motifs=ref_mods_all_motifs)

    return args, ref_out_info


########
# Main #
########

def _main(args):
    try:
        mh.mkdir(args.output_directory, args.overwrite)
    except mh.MegaError as e:
        LOGGER.error(str(e))
        sys.exit(1)

    logging.init_logger(args.output_directory)
    LOGGER.debug('Command: """' + ' '.join(sys.argv) + '"""')
    if _DO_PROFILE:
        LOGGER.warning('Running profiling. This may slow processing.')

    backend_params = backends.parse_backend_params(args)
    with backends.ModelInfo(backend_params, args.processes) as model_info:
        # process ref out first as it might add mods or variants to outputs
        args, ref_out_info = parse_ref_out_args(args, model_info)
        args, mods_info = parse_mod_args(args, model_info)
        # aligner can take a while to load, so load as late as possible
        aligner = parse_aligner_args(args)
        args, vars_data = parse_var_args(args, model_info, aligner)

        process_all_reads(
            args.fast5s_dir, not args.not_recursive, args.num_reads,
            args.read_ids_filename, model_info, args.outputs,
            args.output_directory, args.basecalls_format, aligner, vars_data,
            args.processes, args.verbose_read_progress,
            args.suppress_progress_bars, mods_info, args.database_safety,
            ref_out_info, args.rna, not args.suppress_queues_status)

    if aligner is None:
        # all other tasks require aligner
        return

    ref_fn = aligner.ref_fn
    map_out_fmt = aligner.out_fmt
    del aligner

    # start mapping processes before other post-per-read tasks
    map_p, mod_map_ps, var_map_p, var_sort_fn = start_sort_mapping_procs(
        args.outputs, args.output_directory, map_out_fmt, ref_fn,
        mods_info.mod_long_names)

    if mh.VAR_NAME in args.outputs or mh.MOD_NAME in args.outputs:
        aggregate.aggregate_stats(
            args.outputs, args.output_directory, args.processes,
            args.write_vcf_log_probs, args.heterozygous_factors,
            vars_data.call_mode, mods_info.agg_info, args.write_mod_log_probs,
            mods_info.mod_output_fmts, args.suppress_progress_bars)

    variant_fn = index_variant_fn = None
    if mh.VAR_NAME in args.outputs:
        LOGGER.info('Sorting output variant file')
        variant_fn = mh.get_megalodon_fn(args.output_directory, mh.VAR_NAME)
        sort_variant_fn = mh.add_fn_suffix(variant_fn, 'sorted')
        variants.sort_variants(variant_fn, sort_variant_fn)
        LOGGER.info('Indexing output variant file')
        index_variant_fn = variants.index_variants(sort_variant_fn)

    get_map_procs(
        map_p, mod_map_ps, var_map_p, var_sort_fn, index_variant_fn,
        variant_fn)


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
