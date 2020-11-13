import os
import sys
import queue
import threading
import traceback
from time import sleep
import multiprocessing as mp
from itertools import product
from collections import defaultdict, OrderedDict

import mappy
import numpy as np
from tqdm import tqdm

from megalodon._version import MEGALODON_VERSION
from megalodon import (
    aggregate, backends, fast5_io, logging, mapping, mods,
    variants, megalodon_helper as mh, megalodon_multiprocessing as mega_mp)


LOGGER = logging.get_logger()
# set blas library environment variables (without these the cblas calls
# can completely halt processing)
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'

_DO_PROFILE = False
_UNEXPECTED_ERROR_CODE = 'Unexpected error'
_UNEXPECTED_ERROR_FN = 'unexpected_megalodon_errors.{}.err'
_SIG_EXTRACT_GETTER_NAME = 'extract_signal'
_FAILED_READ_GETTER_NAME = 'failed_reads'
_MAX_NUM_UNEXP_ERRORS = 50
DO_INTERPOLATE_SIG_POS = False
# update error text when 10% more errors are found
ERR_UPDATE_PROP = 0.1


############################
# Post Per-read Processing #
############################

def post_process_mapping(map_bn, map_fmt, ref_fn, samtools_exec):
    map_fn = '{}.{}'.format(map_bn, map_fmt)
    map_sort_fn = '{}.sorted.{}'.format(map_bn, map_fmt)
    map_p = mp.Process(
        target=mapping.sort_and_index_mapping,
        args=(samtools_exec, map_fn, map_sort_fn, map_fmt, ref_fn),
        daemon=True)
    map_p.start()

    return map_p, map_sort_fn


def start_sort_mapping_procs(map_info, mods_info, vars_info):
    map_p = mod_map_ps = var_map_p = var_sort_fn = None
    if map_info.do_output_mappings and map_info.do_sort_mappings:
        LOGGER.info('Spawning process to sort mappings')
        map_p, _ = post_process_mapping(
            mh.get_megalodon_fn(map_info.out_dir, mh.MAP_NAME),
            map_info.map_fmt, map_info.cram_ref_fn, map_info.samtools_exec)
    if mods_info.do_output.mod_map and map_info.do_sort_mappings:
        LOGGER.info('Spawning process to sort modified base mappings')
        mod_map_bn = mh.get_megalodon_fn(
            mods_info.out_dir, mh.MOD_MAP_NAME)
        if mods_info.map_emulate_bisulfite:
            mod_map_bns = [
                '{}.{}'.format(mod_map_bn, mln)
                for mod_base, mln in mods_info.mod_long_names]
        else:
            mod_map_bns = [mod_map_bn]
        mod_map_ps = [post_process_mapping(
            mod_map_bn, map_info.map_fmt, map_info.cram_ref_fn,
            map_info.samtools_exec)[0] for mod_map_bn in mod_map_bns]
    if vars_info.do_output.var_map and map_info.do_sort_mappings:
        LOGGER.info('Spawning process to sort variant mappings')
        var_map_p, var_sort_fn = post_process_mapping(
            mh.get_megalodon_fn(vars_info.out_dir, mh.VAR_MAP_NAME),
            map_info.map_fmt, map_info.cram_ref_fn, map_info.samtools_exec)
    return map_p, mod_map_ps, var_map_p, var_sort_fn


def get_map_procs(
        map_p, mod_map_ps, var_map_p, var_sort_fn, index_variant_fn,
        variant_fn):
    if var_map_p is not None:
        if var_map_p.is_alive():
            LOGGER.info('Waiting for variant mappings sort')
            var_map_p.join()
        if index_variant_fn is not None and var_sort_fn is not None:
            LOGGER.info(variants.get_whatshap_command(
                index_variant_fn, var_sort_fn))
    if mod_map_ps is not None:
        if any(mod_map_p.is_alive() for mod_map_p in mod_map_ps):
            LOGGER.info('Waiting for modified base mappings sort')
            for mod_map_p in mod_map_ps:
                mod_map_p.join()
    if map_p is not None:
        if map_p.is_alive():
            LOGGER.info('Waiting for mappings sort')
            map_p.join()


###################
# Read Processing #
###################

def handle_errors(func, args, r_vals, out_q, fast5_fn, failed_reads_q):
    try:
        out_q.put((func(*args), r_vals))
    except KeyboardInterrupt:
        failed_reads_q.put(tuple(mh.READ_STATUS(
            is_err=True, do_update_prog=False, err_type='Keyboard interrupt',
            fast5_fn=fast5_fn)))
        return
    except mh.MegaError as e:
        failed_reads_q.put(tuple(mh.READ_STATUS(
            is_err=True, do_update_prog=False, err_type=str(e),
            fast5_fn=fast5_fn)))
    except Exception:
        failed_reads_q.put(tuple(mh.READ_STATUS(
            is_err=True, do_update_prog=False, err_type=_UNEXPECTED_ERROR_CODE,
            fast5_fn=fast5_fn, err_tb=traceback.format_exc())))


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
        getter_qpcs, caller_conn, bc_res, ref_out_info, vars_info, mods_info,
        bc_info):
    """ Workhorse per-read megalodon function (connects all the parts)
    """
    (sig_info, seq_summ_info, called_read, rl_cumsum, can_post, post_w_mods,
     mods_scores) = bc_res
    if bc_info.do_output.any:
        # convert seq_summ_info to tuple since namedtuples can't be
        # pickled for passing through a queue.
        if bc_info.rev_sig:
            # sequence is stored internally in sequencing direction. Send
            # to basecall output in reference direction.
            getter_qpcs[mh.BC_NAME].queue.put((
                sig_info.read_id, called_read.seq[::-1],
                called_read.qual[::-1], mods_scores, tuple(seq_summ_info)))
        else:
            getter_qpcs[mh.BC_NAME].queue.put((
                sig_info.read_id, called_read.seq, called_read.qual,
                mods_scores, tuple(seq_summ_info)))

    # if no mapping connection return after basecalls are passed out
    if caller_conn is None:
        return

    # map read and record mapping from reference to query positions
    map_q = getter_qpcs[mh.MAP_NAME].queue \
        if mh.MAP_NAME in getter_qpcs else None
    r_ref_seq, r_to_q_poss, r_ref_pos, r_cigar = mapping.map_read(
        caller_conn, called_read, sig_info, map_q, bc_info.rev_sig, rl_cumsum)
    np_ref_seq = mh.seq_to_int(r_ref_seq, error_on_invalid=False)

    failed_reads_q = getter_qpcs[_FAILED_READ_GETTER_NAME].queue
    sig_map_res = None
    if ref_out_info.do_output.sig_maps:
        pass_sig_map_filts = mapping.read_passes_filters(
            ref_out_info.filt_params, len(called_read.seq),
            r_ref_pos.q_trim_start, r_ref_pos.q_trim_end, r_cigar)
        sig_map_res = signal_mapping.SIG_MAP_RESULT(
            pass_sig_map_filts, sig_info.fast5_fn, sig_info.dacs,
            sig_info.scale_params, r_ref_seq, sig_info.stride,
            sig_info.read_id, r_to_q_poss, rl_cumsum, r_ref_pos, ref_out_info)
        if ref_out_info.do_output.can_sig_maps and pass_sig_map_filts:
            try:
                getter_qpcs[mh.SIG_MAP_NAME].queue.put(
                    signal_mapping.get_remapping(*sig_map_res[1:]))
            except Exception as e:
                LOGGER.debug('{} SignalMappingError {}'.format(
                    sig_info.read_id, str(e)))
                # taiyaki errors can contain newlines so split them here
                failed_reads_q.put(tuple(mh.READ_STATUS(
                    is_err=True, do_update_prog=False,
                    err_type=' ::: '.join(str(e).strip().split('\n')),
                    fast5_fn=sig_info.fast5_fn)))

    # get mapped start in post and run len to mapped bit of output
    post_mapped_start, post_mapped_end = (rl_cumsum[r_ref_pos.q_trim_start],
                                          rl_cumsum[r_ref_pos.q_trim_end])
    mapped_rl_cumsum = rl_cumsum[
        r_ref_pos.q_trim_start:r_ref_pos.q_trim_end + 1] - post_mapped_start
    if DO_INTERPOLATE_SIG_POS:
        ref_to_block = interpolate_sig_pos(r_to_q_poss, mapped_rl_cumsum)
    else:
        ref_to_block = mapped_rl_cumsum[r_to_q_poss]

    if vars_info.do_output.db:
        assert not bc_info.rev_sig, (
            'Reversed raw signal (RNA) not compatible with sequence ' +
            'variant detection.')
        mapped_can_post = can_post[post_mapped_start:post_mapped_end]
        handle_errors(
            func=variants.call_read_vars,
            args=(vars_info, r_ref_pos, np_ref_seq, ref_to_block,
                  mapped_can_post),
            r_vals=(sig_info.read_id, r_ref_pos.chrm, r_ref_pos.strand,
                    r_ref_pos.start, r_ref_seq, len(called_read.seq),
                    r_ref_pos.q_trim_start, r_ref_pos.q_trim_end, r_cigar),
            out_q=getter_qpcs[mh.PR_VAR_NAME].queue,
            fast5_fn=sig_info.fast5_fn + ':::' + sig_info.read_id,
            failed_reads_q=failed_reads_q)
    if mods_info.do_output.any:
        mapped_post_w_mods = post_w_mods[post_mapped_start:post_mapped_end]
        mod_sig_map_q = getter_qpcs[mh.SIG_MAP_NAME].queue \
            if ref_out_info.do_output.mod_sig_maps else None
        handle_errors(
            func=mods.call_read_mods,
            args=(r_ref_pos, r_ref_seq, ref_to_block, mapped_post_w_mods,
                  mods_info, mod_sig_map_q, sig_map_res, bc_info.rev_sig,
                  sig_info.read_id, failed_reads_q, sig_info.fast5_fn),
            r_vals=(sig_info.read_id, r_ref_pos.chrm, r_ref_pos.strand,
                    r_ref_pos.start, r_ref_seq, len(called_read.seq),
                    r_ref_pos.q_trim_start, r_ref_pos.q_trim_end, r_cigar),
            out_q=getter_qpcs[mh.PR_MOD_NAME].queue,
            fast5_fn=sig_info.fast5_fn + ':::' + sig_info.read_id,
            failed_reads_q=failed_reads_q)


########################
# Process reads worker #
########################

def _process_reads_worker(
        signal_q, getter_qpcs, caller_conn, ref_out_info, model_info,
        vars_info, mods_info, bc_info, device):
    def create_batch_gen():
        for _ in range(bc_info.reads_per_batch):
            try:
                read_sig_data = signal_q.queue.get(timeout=0.01)
            except queue.Empty:
                continue
            if read_sig_data is None:
                LOGGER.debug('Closing')
                # send signal to end main loop then end this iterator
                gen_conn.send(True)
                break
            sig_info, seq_summ_info = read_sig_data
            # convert tuples back to namedtuples after multiprocessing
            sig_info = backends.SIGNAL_DATA(*sig_info)
            yield sig_info, mh.SEQ_SUMM_INFO(*seq_summ_info)
            LOGGER.debug('{} Processing'.format(sig_info.read_id))

    def iter_bc_res():
        while not full_iter_conn.poll():
            reads_batch_gen = create_batch_gen()
            # perform basecalling using loaded backend
            for bc_res in model_info.iter_basecalled_reads(
                    reads_batch_gen,
                    return_post_w_mods=mods_info.do_output.any,
                    return_mod_scores=bc_info.do_output.mod_basecalls,
                    update_sig_info=ref_out_info.do_output.sig_maps,
                    signal_reversed=bc_info.rev_sig,
                    mod_bc_min_prob=bc_info.mod_bc_min_prob,
                    failed_reads_q=failed_reads_q):
                yield bc_res

    # wrap process prep in try loop to avoid stalled command
    try:
        # prepare connection for communication with batch generator
        full_iter_conn, gen_conn = mp.Pipe(duplex=False)
        failed_reads_q = getter_qpcs[_FAILED_READ_GETTER_NAME].queue
        model_info.prep_model_worker(device)
        vars_info.reopen_variant_index()
        LOGGER.debug('Starting')
    except Exception:
        LOGGER.debug('InitFailed traceback: {}'.format(traceback.format_exc()))
        return

    for bc_res in iter_bc_res():
        sig_info = bc_res[0]
        try:
            process_read(
                getter_qpcs, caller_conn, bc_res, ref_out_info, vars_info,
                mods_info, bc_info)
            failed_reads_q.put(
                tuple(mh.READ_STATUS(n_sig=sig_info.raw_len)))
            LOGGER.debug('{} Success'.format(sig_info.read_id))
        except KeyboardInterrupt:
            failed_reads_q.put(tuple(mh.READ_STATUS(
                is_err=True, do_update_prog=False,
                err_type='Keyboard interrupt', fast5_fn=sig_info.fast5_fn)))
            LOGGER.debug('{} KeyboardInterrupt'.format(sig_info.read_id))
            return
        except mh.MegaError as e:
            raw_len = sig_info.raw_len if hasattr(sig_info, 'raw_len') else 0
            fn_rid = '{}:::{}'.format(sig_info.fast5_fn, sig_info.read_id)
            failed_reads_q.put(tuple(mh.READ_STATUS(
                is_err=True, do_update_prog=True, err_type=str(e),
                fast5_fn=fn_rid, n_sig=raw_len)))
            LOGGER.debug('{} Failed {}'.format(sig_info.read_id, str(e)))
        except Exception as e:
            fn_rid = '{}:::{}'.format(sig_info.fast5_fn, sig_info.read_id)
            failed_reads_q.put(tuple(mh.READ_STATUS(
                is_err=True, do_update_prog=True,
                err_type=_UNEXPECTED_ERROR_CODE, fast5_fn=fn_rid,
                err_tb=traceback.format_exc())))
            LOGGER.debug(
                '{} UnexpectedFail {}'.format(sig_info.read_id, str(e)))


if _DO_PROFILE:
    _process_reads_wrapper = _process_reads_worker

    def _process_reads_worker(*args):
        import cProfile
        cProfile.runctx('_process_reads_wrapper(*args)', globals(), locals(),
                        filename='read_processing.prof')


#################
# Queue getters #
#################

def _get_bc_queue(bc_q, bc_conn, bc_info, aux_failed_q):
    def write_read(bc_res):
        read_id, r_seq, r_qual, mods_scores, seq_summ_info = bc_res
        # convert seq_summ_info back into namedtuple after passing
        # through mp.queue
        seq_summ_info = mh.SEQ_SUMM_INFO(*seq_summ_info)
        if bc_info.do_output.basecalls:
            if bc_info.bc_fmt == mh.BC_FMT_FQ:
                if r_qual is None:
                    r_qual = '!' * len(r_seq)
                bc_fp.write('@{}\n{}\n+\n{}\n'.format(read_id, r_seq, r_qual))
            else:
                bc_fp.write('>{}\n{}\n'.format(read_id, r_seq))
            seq_summ_fp.write('\t'.join(map(str, seq_summ_info)) + '\n')

        if bc_info.do_output.mod_basecalls:
            # 4 indicates unmapped
            mods_fp.write(mapping.prepare_mapping(
                read_id, r_seq, qual=[ord(q) - 33 for q in r_qual],
                mods_scores=mods_scores, flag=4))

    try:
        LOGGER.debug('GetterStarting')
        if bc_info.do_output.basecalls:
            bc_fp = open('{}.{}'.format(mh.get_megalodon_fn(
                bc_info.out_dir, mh.BC_NAME), bc_info.bc_fmt), 'w')
            seq_summ_fp = open(
                mh.get_megalodon_fn(bc_info.out_dir, mh.SEQ_SUMM_NAME), 'w')
            seq_summ_fp.write('\t'.join(mh.SEQ_SUMM_INFO._fields) + '\n')
        if bc_info.do_output.mod_basecalls:
            LOGGER.debug('outputting mod basecalls')
            mods_fp = mapping.open_unaligned_alignment_file(
                mh.get_megalodon_fn(bc_info.out_dir, mh.BC_MODS_NAME),
                bc_info.mod_bc_fmt, bc_info.mod_long_names)
        workers_active = True
        LOGGER.debug('GetterInitComplete')
    except Exception as e:
        aux_failed_q.put((
            'BasecallsInitError', str(e), traceback.format_exc()))
        return

    try:
        while workers_active or not bc_q.empty():
            try:
                bc_res = bc_q.get(timeout=0.1)
                mh.log_errors(write_read, bc_res)
            except queue.Empty:
                if bc_conn.poll():
                    workers_active = False
        LOGGER.debug('GetterClosing')
    except Exception as e:
        aux_failed_q.put((
            'BasecallsProcessingError', str(e), traceback.format_exc()))
    finally:
        if bc_info.do_output.basecalls:
            bc_fp.close()
        if bc_info.do_output.mod_basecalls:
            mods_fp.close()


####################
# Dynamic progress #
####################

def iter_most_common_errs(err_types, reads_called, num_errs=None):
    summ_errs = sorted(err_types)[::-1]
    if num_errs is not None:
        summ_errs = summ_errs[:num_errs]
    if reads_called == 0:
        summ_errs = []
    # skip errors that were checked but not registered
    summ_errs = [(n_fns, err) for n_fns, err in summ_errs if n_fns > 0]
    for n_fns, err in summ_errs:
        yield '{:8.1f}% ({:>7} reads) : {:<80}'.format(
            100 * n_fns / float(reads_called), n_fns, err)
    if num_errs is not None and len(summ_errs) < num_errs:
        for _ in range(num_errs - len(summ_errs)):
            yield '    -----'


def update_err(failed_reads, err_type, fast5_fn, err_tb, unexp_err_fp):
    failed_reads[err_type].append(fast5_fn)
    if err_type == _UNEXPECTED_ERROR_CODE:
        # if this is the first unexpected error open file
        if len(failed_reads[_UNEXPECTED_ERROR_CODE]) == 1:
            unexp_err_fp = open(_UNEXPECTED_ERROR_FN.format(
                np.random.randint(10000)), 'w')
        if len(failed_reads[err_type]) >= _MAX_NUM_UNEXP_ERRORS:
            # if this is the unexpected error limit close the file
            unexp_err_fp.close()
        else:
            # else write the error
            unexp_err_fp.write(
                fast5_fn + '\n:::\n' + err_tb + '\n\n\n')
            unexp_err_fp.flush()
    return unexp_err_fp


def update_prog(
        err_statuses, bar, q_bars, sig_called, getter_qpcs, status_info,
        reads_called, failed_reads, last_err_write, read_called=True):
    try:
        bar.set_postfix(
            {'samples/s': sig_called / bar.format_dict['elapsed']},
            refresh=False)
    except AttributeError:
        # sometimes get no format_dict error, if so don't include samples/s
        pass
    if q_bars is not None:
        for q_name, q_bar in q_bars.items():
            q_bar.n = getter_qpcs[q_name].queue.qsize()
            # trigger display refresh
            q_bar.update(0)
    if read_called:
        bar.update()
    if status_info.num_prog_errs > 0:
        err_types = [
            (len(fns), err) for err, fns in failed_reads.items()]
        num_errs = sum((x[0] for x in err_types))
        if num_errs > 0 and (
                last_err_write == 0 or
                num_errs / last_err_write > 1 + ERR_UPDATE_PROP):
            last_err_write = num_errs
            for err_num, err_status in enumerate(iter_most_common_errs(
                    err_types, reads_called, status_info.num_prog_errs)):
                err_statuses[err_num].set_description_str(err_status)

    return last_err_write


def prep_errors_bar(status_info, getter_qpcs):
    if status_info.suppress_prog_bar:
        return None, None, None

    # prep queue capacity status bars if requested
    q_labs = None
    if not status_info.suppress_queues:
        def io_str(q_name):
            if q_name == _SIG_EXTRACT_GETTER_NAME:
                return ' input queue capacity {}'
            return 'output queue capacity {}'

        sys.stderr.write(
            'Full output or empty input queues indicate I/O bottleneck\n')
        q_labs = [
            (q_num, q_name, io_str(q_name).format(q_name))
            for q_num, q_name in enumerate([
                q_name for q_name in getter_qpcs
                if q_name != _FAILED_READ_GETTER_NAME])]

    err_statuses = q_bars = None
    if status_info.num_prog_errs > 0:
        # write header for progress bars if dynamic error reporting is on
        sys.stderr.write(
            '{} most common unsuccessful processing stages:\n'.format(
                status_info.num_prog_errs))
        err_statuses = [
            tqdm(position=err_num, bar_format='{desc}', desc='    -----')
            for err_num in range(status_info.num_prog_errs)]
    # open main progress bar
    bar = tqdm(total=None, smoothing=0, initial=0,
               unit='reads', dynamic_ncols=True,
               position=status_info.num_prog_errs, desc='Read Processing')
    if q_labs is not None:
        q_bars = OrderedDict((q_name, tqdm(
            desc=q_lab, total=mh._MAX_QUEUE_SIZE, smoothing=0,
            dynamic_ncols=True,
            position=q_num + 1 + status_info.num_prog_errs,
            bar_format='{desc: <42}: {percentage:3.0f}%|{bar}| ' +
            '{n_fmt}/{total_fmt}')) for q_num, q_name, q_lab in q_labs)

    return err_statuses, bar, q_bars


def _get_fail_queue(
        failed_reads_q, f_conn, status_info, gnr_conn, getter_qpcs,
        signal_q):
    LOGGER.info('Processing reads')
    LOGGER.debug('GetterStarting')
    reads_called = sig_called = last_err_write = 0
    unexp_err_fp = None
    failed_reads = defaultdict(list)
    getter_qpcs[_SIG_EXTRACT_GETTER_NAME] = signal_q
    getter_qpcs.move_to_end(_SIG_EXTRACT_GETTER_NAME, last=False)
    err_statuses, bar, q_bars = prep_errors_bar(status_info, getter_qpcs)
    workers_active = True
    LOGGER.debug('GetterInitComplete')
    try:
        while workers_active or not failed_reads_q.empty():
            try:
                read_status = failed_reads_q.get(timeout=0.1)
                read_status = mh.READ_STATUS(*read_status)
                sig_called += read_status.n_sig
                if read_status.do_update_prog:
                    reads_called += 1
                if read_status.is_err:
                    unexp_err_fp = update_err(
                        failed_reads, read_status.err_type,
                        read_status.fast5_fn, read_status.err_tb, unexp_err_fp)
                if read_status.do_update_prog and \
                   not status_info.suppress_prog_bar:
                    last_err_write = update_prog(
                        err_statuses, bar, q_bars, sig_called, getter_qpcs,
                        status_info, reads_called, failed_reads,
                        last_err_write)
                # get total number of reads once all reads are enumerated
                if bar is not None and bar.total is None:
                    if gnr_conn.poll():
                        bar.total = gnr_conn.recv()
            except queue.Empty:
                if f_conn.poll():
                    workers_active = False
        LOGGER.debug('GetterClosing')
    except KeyboardInterrupt:
        # exit gracefully on keyboard inturrupt
        return

    # cleanup progressbars
    if not status_info.suppress_prog_bar:
        # wait for getter queues to flush
        if q_bars is not None:
            while any(not getter_qpcs[q_name].queue.empty()
                      for q_name in q_bars.keys()):
                for q_name, q_bar in q_bars.items():
                    q_bar.n = getter_qpcs[q_name].queue.qsize()
                    q_bar.update(0)
        # close all open progressbars
        if err_statuses is not None:
            for err_status in err_statuses:
                err_status.close()
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
        err_types = [(len(fns), err) for err, fns in failed_reads.items()]
        LOGGER.info('Unsuccessful processing types:\n{}'.format('\n'.join(
            iter_most_common_errs(err_types, reads_called))))
        # TODO flag to output failed read names to file
    else:
        LOGGER.info('All reads processed successfully')


#######################
# All read processing #
#######################

def wait_for_completion(
        files_p, extract_sig_ps, signal_q, proc_reads_ps, map_read_ts,
        getter_qpcs, aux_failed_q, input_info):
    def kill_all_proc():
        for p in [files_p, ] + proc_reads_ps:
            if p.is_alive():
                p.terminate()
                p.join()
        for q in list(getter_qpcs.values()):
            if q.proc.is_alive():
                q.proc.terminate()
                q.proc.join()

    try:
        # wait for file enumeration process to finish first
        while files_p.is_alive():
            try:
                aux_err = aux_failed_q.get(block=False)
                kill_all_proc()
                sleep(0.01)
                for msg in aux_err:
                    LOGGER.error(msg)
                sys.exit(1)
            except queue.Empty:
                # TODO check for failed workers and create mechanism to restart
                sleep(1)
        LOGGER.debug('JoiningMain: FileFiller')
        files_p.join()
        LOGGER.debug('JoinedMain: FileFiller')

        # wait for signal extraction to finish next
        while any(p.is_alive() for p in extract_sig_ps):
            try:
                aux_err = aux_failed_q.get(block=False)
                kill_all_proc()
                sleep(0.01)
                for msg in aux_err:
                    LOGGER.error(msg)
                sys.exit(1)
            except queue.Empty:
                # TODO check for failed workers and create mechanism to restart
                sleep(1)
        LOGGER.debug('JoiningMain: SignalExtractors')
        for extract_sig_p in extract_sig_ps:
            extract_sig_p.join()
        LOGGER.debug('JoinedMain: SignalExtractors')
        # add None to indicate to read workers that signal extraction is
        # complete
        for _ in range(input_info.num_ps):
            signal_q.queue.put(None)

        # wait for signal queue to be exhausted
        while not signal_q.queue.empty():
            try:
                aux_err = aux_failed_q.get(block=False)
                # if an auxiliary process fails exit megalodon
                kill_all_proc()
                sleep(0.01)
                for msg in aux_err:
                    LOGGER.error(msg)
                sys.exit(1)
            except queue.Empty:
                # check that a getter queue has not failed with a segfault
                for g_name, getter_qpc in getter_qpcs.items():
                    if not getter_qpc.proc.is_alive() and \
                       not signal_q.queue.empty():
                        kill_all_proc()
                        sleep(0.01)
                        LOGGER.error((
                            '{} Getter queue has unexpectedly died likely ' +
                            'via a segfault error. Please log this ' +
                            'issue.').format(g_name))
                        sys.exit(1)
                # TODO check for failed workers and create mechanism to restart
                sleep(1)

        LOGGER.debug('JoiningMain: Workers')
        # join worker/getter processes
        for proc_reads_p in proc_reads_ps:
            LOGGER.debug('JoiningMain: {}'.format(proc_reads_p.name))
            proc_reads_p.join()
            LOGGER.debug('JoinedMain: {}'.format(proc_reads_p.name))
        LOGGER.debug('JoiningMain: Mappers')
        if map_read_ts is not None:
            for map_t in map_read_ts:
                LOGGER.debug('JoiningMain: {}'.format(map_t.name))
                map_t.join()
                LOGGER.debug('JoinedMain: {}'.format(map_t.name))
        for out_name, getter_qpc in getter_qpcs.items():
            getter_qpc.conn.send(True)
        LOGGER.debug('JoiningMain: Getters')
        # comm to getter processes to return
        for out_name, getter_qpc in getter_qpcs.items():
            if getter_qpc.proc.is_alive():
                if out_name == mh.PR_VAR_NAME:
                    LOGGER.info(
                        'Waiting for variants database to complete indexing')
                elif out_name == mh.PR_MOD_NAME:
                    LOGGER.info(
                        'Waiting for mods database to complete indexing')
                LOGGER.debug('JoiningMain: {}'.format(out_name))
                getter_qpc.proc.join()
                LOGGER.debug('JoinedMain: {}'.format(out_name))
        LOGGER.debug('JoiningMain: Complete')

    except KeyboardInterrupt:
        LOGGER.error('Exiting due to keyboard interrupt')
        sys.exit(1)


def process_all_reads(
        status_info, input_info, model_info, bc_info, aligner, map_info,
        mods_info, vars_info, ref_out_info):
    # queue to communicate with the main process that an auxiliary/non-worker
    # process has failed unexpectedly to trigger full shutdown
    aux_failed_q = mp.Queue()
    LOGGER.info('Preparing workers to process reads')
    # read filename queue filler
    # Note no maxsize for this queue to compute total number of reads while
    # also not delaying read processing
    fn_read_ids_q = mega_mp.CountingMPQueue()
    nr_conn, gnr_conn = mp.Pipe()
    files_p = mp.Process(
        target=fast5_io._fill_files_queue, daemon=True, name='FileFiller',
        args=(input_info, fn_read_ids_q, nr_conn, aux_failed_q))
    files_p.start()
    # set maxsize to limit memory usage from loaded raw signal
    signal_q = mega_mp.GETTER_QPC(
        queue=mega_mp.CountingMPQueue(maxsize=mh._MAX_QUEUE_SIZE),
        proc=None, conn=None)
    extract_sig_ps = []
    for es_i in range(input_info.num_extract_sig_proc):
        extract_sig_ps.append(mp.Process(
            target=fast5_io._extract_signal, daemon=True,
            name='SignalExtractor{:03d}'.format(es_i),
            args=(fn_read_ids_q, signal_q.queue, aux_failed_q, input_info,
                  model_info, ref_out_info.do_output.sig_maps)))
        extract_sig_ps[-1].start()

    # collate information about the output/getter processes
    getters_info = [
        mh.GETTER_INFO(
            name=mh.BC_NAME, do_output=bc_info.do_output.any,
            func=_get_bc_queue, args=(bc_info, aux_failed_q)),
        mh.GETTER_INFO(
            name=mh.MAP_NAME, do_output=map_info.do_output_mappings,
            func=mapping._get_map_queue,
            args=(map_info, ref_out_info, aux_failed_q)),
        mh.GETTER_INFO(
            name=mh.PR_VAR_NAME, do_output=vars_info.do_output.db,
            func=variants._get_variants_queue,
            args=(vars_info, ref_out_info, map_info, aux_failed_q)),
        mh.GETTER_INFO(
            name=mh.PR_MOD_NAME, do_output=mods_info.do_output.any,
            func=mods._get_mods_queue,
            args=(mods_info, map_info, ref_out_info, aux_failed_q)),
        mh.GETTER_INFO(
            name=mh.SIG_MAP_NAME, do_output=ref_out_info.do_output.sig_maps,
            func=ref_out_info.get_sig_map_func,
            args=(ref_out_info, aux_failed_q))]
    getter_qpcs = OrderedDict(
        (gi.name, mega_mp.create_getter_qpc(gi.func, gi.args, name=gi.name))
        for gi in getters_info if gi.do_output)
    # failed reads queue needs access to other getter queues
    getter_qpcs[_FAILED_READ_GETTER_NAME] = mega_mp.create_getter_qpc(
        _get_fail_queue, (status_info, gnr_conn, getter_qpcs, signal_q),
        max_size=None, name=_FAILED_READ_GETTER_NAME)

    proc_reads_ps, map_conns = [], []
    for pnum, device in enumerate(model_info.process_devices):
        map_conn, caller_conn = (None, None) if aligner is None else mp.Pipe()
        map_conns.append(map_conn)
        proc_reads_ps.append(mp.Process(
            target=_process_reads_worker, args=(
                signal_q, getter_qpcs, caller_conn, ref_out_info,
                model_info, vars_info, mods_info, bc_info, device),
            daemon=True, name='ReadWorker{:03d}'.format(pnum)))
        proc_reads_ps[-1].start()
        if caller_conn is not None:
            caller_conn.close()
        del caller_conn
    # ensure process all start up before initializing mapping threads
    sleep(0.1)
    # perform mapping in threads for mappy shared memory interface
    # open threads after all processes have started due to python
    # multiprocess combined with threading instability
    map_read_ts = None if aligner is None else []
    if aligner is not None:
        for ti, map_conn in enumerate(map_conns):
            map_read_ts.append(threading.Thread(
                target=mapping._map_read_worker, args=(aligner, map_conn),
                daemon=True, name='Mapper{:03d}'.format(ti)))
            map_read_ts[-1].start()

    wait_for_completion(
        files_p, extract_sig_ps, signal_q, proc_reads_ps, map_read_ts,
        getter_qpcs, aux_failed_q, input_info)


############################
# Input validation/parsing #
############################

def parse_aligner_args(args):
    if len(mh.ALIGN_OUTPUTS.intersection(args.outputs)) > 0:
        if args.reference is None:
            LOGGER.error(
                ('Output(s) requiring reference alignment requested ({}), ' +
                 'but --reference not provided.').format(', '.join(
                     mh.ALIGN_OUTPUTS.intersection(args.outputs))))
            sys.exit(1)
        LOGGER.info('Loading reference')
        if not (os.path.exists(args.reference) and
                os.path.isfile(args.reference)):
            LOGGER.error('Provided reference file does not exist or is ' +
                         'not a file.')
            sys.exit(1)
        aligner = mappy.Aligner(
            str(args.reference), preset=str('map-ont'), best_n=1)
    else:
        aligner = None
        if args.reference is not None:
            LOGGER.warning(
                '[--reference] provided, but no [--outputs] requiring ' +
                'alignment was requested. Argument will be ignored.')
    map_info = mapping.MapInfo(
        aligner=aligner, map_fmt=args.mappings_format, ref_fn=args.reference,
        out_dir=args.output_directory,
        do_output_mappings=mh.MAP_NAME in args.outputs,
        samtools_exec=args.samtools_executable,
        do_sort_mappings=args.sort_mappings, cram_ref_fn=args.cram_reference)
    if map_info.do_output_mappings:
        try:
            map_info.test_open_alignment_out_file()
        except mh.MegaError:
            sys.exit(1)
        if map_info.do_sort_mappings:
            map_info.test_samtools()
    return aligner, map_info


def parse_var_args(args, model_info, aligner, ref_out_info):
    if args.ref_include_variants and mh.PR_VAR_NAME not in args.outputs:
        LOGGER.warning('--ref-include-variants set, so adding ' +
                       '"per_read_vars" to --outputs.')
        args.outputs.append(mh.PR_VAR_NAME)
    if mh.VAR_MAP_NAME in args.outputs and \
       mh.PR_VAR_NAME not in args.outputs:
        LOGGER.warning((
            'Adding "{}" to --outputs since "{}" was requested. For full ' +
            'phased variant pipeline add "{}" or run aggregation after run ' +
            'is complete.').format(
                mh.PR_VAR_NAME, mh.VAR_MAP_NAME, mh.VAR_NAME))
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
    skip_db_index = args.skip_database_index
    if skip_db_index and mh.PR_VAR_NAME in args.outputs:
        LOGGER.warning(
            'Database index skipping is not currently implemented for ' +
            'variants output. Ignoring --skip-database-index.')
        skip_db_index = False
    if skip_db_index and mh.VAR_NAME in args.outputs:
        LOGGER.warning(
            'Cannot skip database indexing when aggregated output ' +
            '"variants" is requested. Ignoring --skip-database-index.')
        skip_db_index = False

    var_calib_fn = mh.get_var_calibration_fn(
        model_info.params.pyguppy.config, args.variant_calibration_filename,
        args.disable_variant_calibration) \
        if mh.PR_VAR_NAME in args.outputs else None
    do_output = mh.VAR_DO_OUTPUT(
        db=mh.PR_VAR_NAME in args.outputs, text=args.write_variants_text,
        var_map=mh.VAR_MAP_NAME in args.outputs)

    try:
        vars_info = variants.VarInfo(
            args.variant_filename, aligner, max_indel_size=args.max_indel_size,
            all_paths=args.variant_all_paths,
            context_bases=args.variant_context_bases,
            vars_calib_fn=var_calib_fn,
            call_mode=variants.HAPLIOD_MODE if args.haploid else
            variants.DIPLOID_MODE, edge_buffer=args.edge_buffer,
            context_min_alt_prob=args.context_min_alt_prob,
            loc_index_in_memory=not args.variant_locations_on_disk,
            variants_are_atomized=args.variants_are_atomized,
            db_safety=args.database_safety, do_output=do_output,
            out_dir=args.output_directory, skip_db_index=skip_db_index)
    except mh.MegaError as e:
        # catch potential errors reading in variant file
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
    return args, vars_info


def parse_mod_args(args, model_info, ref_out_info, map_info):
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
    if model_info.is_cat_mod and len(mh.MOD_OUTPUTS.intersection(
            args.outputs)) == 0:
        LOGGER.warning(
            ('Model supporting modified base calling specified, but no ' +
             'modified base outputs requested. This is fine; just wanted to ' +
             'let you know.'))
    if mh.BC_NAME not in args.outputs and mh.BC_MODS_NAME in args.outputs:
        args.outputs.append(mh.BC_NAME)
    if not args.mod_map_emulate_bisulfite and \
       args.mod_map_base_conv is not None:
        LOGGER.warning(
            '--mod-map-base-conv provided, but --mod-map-emulate-bisulfite ' +
            'not set. --mod-map-base-conv will be ignored.')
    skip_db_index = args.skip_database_index
    if skip_db_index and mh.MOD_NAME in args.outputs:
        LOGGER.warning(
            'Cannot skip database indexing when aggregated output "mods" is ' +
            'requested. Ignoring --skip-database-index.')
        skip_db_index = False

    mod_calib_fn = (mh.get_mod_calibration_fn(
        model_info.params.pyguppy.config, args.mod_calibration_filename,
        args.disable_mod_calibration)
        if mh.PR_MOD_NAME in args.outputs else None)
    if args.mod_aggregate_method == mh.MOD_EM_NAME:
        agg_info = mods.AGG_INFO(mh.MOD_EM_NAME, None)
    elif args.mod_aggregate_method == mh.MOD_EXPIT:
        agg_info = mods.AGG_INFO(mh.MOD_EXPIT, None)
    elif args.mod_aggregate_method == mh.MOD_BIN_THRESH_NAME:
        agg_info = mods.AGG_INFO(
            mh.MOD_BIN_THRESH_NAME, args.mod_binary_threshold)

    do_output = mh.MOD_DO_OUTPUT(
        db=mh.PR_MOD_NAME in args.outputs, text=args.write_mods_text,
        mod_map=mh.MOD_MAP_NAME in args.outputs,
        any=any((mh.PR_MOD_NAME in args.outputs, args.write_mods_text,
                 mh.MOD_MAP_NAME in args.outputs, args.ref_include_mods)))
    mods_info = mods.ModInfo(
        model_info=model_info, all_mod_motifs_raw=args.mod_motif,
        mod_all_paths=args.mod_all_paths,
        mod_context_bases=args.mod_context_bases,
        mods_calib_fn=mod_calib_fn, mod_output_fmts=args.mod_output_formats,
        edge_buffer=args.edge_buffer, agg_info=agg_info,
        mod_thresh=args.ref_mod_threshold,
        do_ann_all_mods=args.ref_include_mods,
        map_emulate_bisulfite=args.mod_map_emulate_bisulfite,
        map_base_conv=args.mod_map_base_conv,
        map_min_prob=args.mod_min_prob,
        mod_db_timeout=args.mod_database_timeout,
        db_safety=args.database_safety, out_dir=args.output_directory,
        skip_db_index=skip_db_index, do_output=do_output)
    if do_output.db:
        # initialize the database tables
        mods.init_mods_db(mods_info, map_info.ref_names_and_lens)
        # load indices and close connection
        mods_db = mods.ModsDb(mods_info.mods_db_fn, read_only=True)
        mods_info.add_mods_db_arrays(mods_db)
        mods_db.close()
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


def parse_ref_out_args(args, model_info, map_info):
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

    # parse reference output filter parameters
    min_len, max_len = (
        (None, None) if args.ref_length_range is None else
        args.ref_length_range)
    ref_filt_params = mh.REF_OUT_FILTER_PARAMS(
        pct_idnt=args.ref_percent_identity_threshold,
        pct_cov=args.ref_percent_coverage_threshold,
        min_len=min_len, max_len=max_len)

    # Determine requested outputs and parse alphabet info
    sig_map_getter = sm_alphabet_info = None
    do_out_csm = do_out_msm = do_out_vsm = False
    if mh.SIG_MAP_NAME in args.outputs:
        LOGGER.info('Loading signal mapping settings.')
        if args.ref_include_mods:
            do_out_msm = True
        elif args.ref_include_variants:
            raise NotImplementedError(
                'Signal mapping with annotated variants not implemented')
            do_out_vsm = True
        else:
            do_out_csm = True
        # import here so that taiyaki is not required unless outputting
        # signal mappings
        from megalodon import signal_mapping
        global signal_mapping
        sig_map_getter = signal_mapping.write_signal_mappings
        sm_alphabet_info = signal_mapping.get_alphabet_info_from_model(
            model_info)
        if args.ref_include_mods and mh.PR_MOD_NAME not in args.outputs:
            LOGGER.warning(('--ref-include-mods set, so adding ' +
                            '"{}" to --outputs.').format(mh.PR_MOD_NAME))
            args.outputs.append(mh.PR_MOD_NAME)
    do_out_cpr = do_out_mpr = do_out_vpr = False
    if mh.PR_REF_NAME in args.outputs:
        LOGGER.info('Loading per-read reference output settings.')
        if args.ref_include_mods:
            do_out_mpr = True
        elif args.ref_include_variants:
            do_out_vpr = True
        else:
            do_out_cpr = True
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
    do_output_pr_refs = any((do_out_cpr, do_out_mpr, do_out_vpr))
    do_output_sig_maps = any((do_out_csm, do_out_msm, do_out_vsm))
    ref_outputs = mh.REF_DO_OUTPUT(
        pr_refs=do_output_pr_refs, can_pr_refs=do_out_cpr,
        mod_pr_refs=do_out_mpr, var_pr_refs=do_out_vpr,
        sig_maps=do_output_sig_maps, can_sig_maps=do_out_csm,
        mod_sig_maps=do_out_msm, var_sig_maps=do_out_vsm)

    # warn if irrelevent arguments are provided
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

    # if motif based mod markup is requested parse here
    per_site_threshs = None
    if args.mod_per_site_threshold is not None:
        LOGGER.info('Loading per-site thresholds.')
        per_site_threshs = mh.parse_bed_scores_np(
            args.mod_per_site_threshold, map_info.ref_names_and_lens)

    ref_mods_all_motifs = None
    if args.ref_mods_all_motifs is not None:
        ref_mods_all_motifs, sm_alphabet_info = parse_ref_mods_all_motifs(
            args.ref_mods_all_motifs, sm_alphabet_info)

    ref_out_info = mh.REF_OUT_INFO(
        do_output=ref_outputs, filt_params=ref_filt_params,
        ref_mods_all_motifs=ref_mods_all_motifs,
        alphabet_info=sm_alphabet_info, out_dir=args.output_directory,
        get_sig_map_func=sig_map_getter, per_site_threshs=per_site_threshs)

    return args, ref_out_info


def parse_basecall_args(args, mods_info):
    if mh.BC_MODS_NAME in args.outputs and mods_info.nmod_base == 0:
        LOGGER.warning('mod_basecalls requested, but specified model does ' +
                       'not support calling modified bases. Removing ' +
                       'mod_basecalls from outputs.')
        args.outputs.remove(mh.BC_MODS_NAME)
    bc_do_output = mh.BASECALL_DO_OUTPUT(
        any=mh.BC_NAME in args.outputs or mh.BC_MODS_NAME in args.outputs,
        basecalls=mh.BC_NAME in args.outputs,
        mod_basecalls=mh.BC_MODS_NAME in args.outputs)
    return mh.BASECALL_INFO(
        do_output=bc_do_output,
        out_dir=args.output_directory,
        bc_fmt=args.basecalls_format, mod_bc_fmt=args.mappings_format,
        mod_bc_min_prob=args.mod_min_prob,
        mod_long_names=mods_info.mod_long_names, rev_sig=args.rna,
        reads_per_batch=args.reads_per_guppy_batch)


def parse_input_args(args):
    return mh.INPUT_INFO(
        fast5s_dir=args.fast5s_dir, recursive=not args.not_recursive,
        num_reads=args.num_reads, read_ids_fn=args.read_ids_filename,
        num_ps=args.processes, do_it_live=args.live_processing,
        num_read_enum_ts=args.num_read_enumeration_threads,
        num_extract_sig_proc=args.num_extract_signal_processes)


def parse_status_args(args):
    return mh.STATUS_INFO(
        suppress_prog_bar=args.suppress_progress_bars,
        suppress_queues=args.suppress_queues_status,
        num_prog_errs=args.verbose_read_progress)


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
    LOGGER.info('Running Megalodon version {}'.format(MEGALODON_VERSION))
    LOGGER.debug('Command: """' + ' '.join(sys.argv) + '"""')
    if _DO_PROFILE:
        LOGGER.warning('Running profiling. This may slow processing.')

    model_info = None
    try:
        status_info = parse_status_args(args)
        input_info = parse_input_args(args)
        backend_params = backends.parse_backend_params(args)
        model_info = backends.ModelInfo(backend_params, args.processes)
        # aligner can take a while to load, so load as late as possible
        aligner, map_info = parse_aligner_args(args)
        # process ref out here as it might add mods or variants to outputs
        args, ref_out_info = parse_ref_out_args(args, model_info, map_info)
        args, mods_info = parse_mod_args(
            args, model_info, ref_out_info, map_info)
        bc_info = parse_basecall_args(args, mods_info)
        args, vars_info = parse_var_args(
            args, model_info, aligner, ref_out_info)
        process_all_reads(
            status_info, input_info, model_info, bc_info, aligner, map_info,
            mods_info, vars_info, ref_out_info)
    except mh.MegaError as e:
        LOGGER.error(str(e))
        if model_info is not None:
            model_info.close()
        sys.exit(1)
    finally:
        if model_info is not None:
            model_info.close()

    if aligner is None:
        # all other tasks require reference mapping
        return
    # delete aligner to free memory
    del aligner

    # start mapping processes before other post-per-read tasks
    map_p, mod_map_ps, var_map_p, var_sort_fn = start_sort_mapping_procs(
        map_info, mods_info, vars_info)

    if mh.VAR_NAME in args.outputs or mh.MOD_NAME in args.outputs:
        aggregate.aggregate_stats(
            args.outputs, args.output_directory, args.processes,
            args.write_vcf_log_probs, args.heterozygous_factors,
            vars_info.call_mode, mods_info.agg_info, args.write_mod_log_probs,
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
    LOGGER.info('Mega Done')


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
