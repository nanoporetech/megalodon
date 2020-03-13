import os
import sys
import h5py
import queue
import argparse
import threading
import traceback
from time import sleep
import multiprocessing as mp
from collections import defaultdict, OrderedDict

import numpy as np
from tqdm import tqdm
from tqdm._utils import _term_move_up

from megalodon import (
    aggregate, backends, fast5_io, logging, mapping, mods,
    variants, megalodon_helper as mh)
from megalodon._version import MEGALODON_VERSION


# set blas library environment variables (without these the cblas calls
# can completely halt processing)
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'

_DO_PROFILE = False
_UNEXPECTED_ERROR_CODE = 'Unexpected error'
_UNEXPECTED_ERROR_FN = 'unexpected_megalodon_errors.{}.err'
_MAX_NUM_UNEXP_ERRORS = 50
DO_INTERPOLATE_SIG_POS = False

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
    return


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
        vars_data, vars_q, mods_q, mods_info, failed_reads_q):
    """ Workhorse per-read megalodon function (connects all the parts)
    """
    # perform basecalling using loaded backend
    (r_seq, r_qual, rl_cumsum, can_post, sig_info, post_w_mods,
     mods_scores) = model_info.basecall_read(
         sig_info, return_post_w_mods=mods_q is not None,
         return_mod_scores=mods_info.do_output_mods,
         update_sig_info=sig_map_q is not None)
    if bc_q is not None:
        bc_q.put((sig_info.read_id, r_seq, r_qual, mods_scores))
    # if no mapping connection return after basecalls are passed out
    if caller_conn is None:
        return

    # map read and record mapping from reference to query positions
    r_ref_seq, r_to_q_poss, r_ref_pos, r_cigar = mapping.map_read(
        r_seq, sig_info.read_id, caller_conn)
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
                  mods_info, mod_sig_map_q, sig_map_res),
            r_vals=(sig_info.read_id, r_ref_pos.chrm, r_ref_pos.strand,
                    r_ref_pos.start, r_ref_seq, len(r_seq),
                    r_ref_pos.q_trim_start, r_ref_pos.q_trim_end, r_cigar),
            out_q=mods_q,
            fast5_fn=sig_info.fast5_fn + ':::' + sig_info.read_id,
            failed_reads_q=failed_reads_q)

    return


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
        bc_fp.flush()

        if do_output_mods:
            try:
                mods_fp.create_dataset(
                    'Reads/' + read_id, data=mods_scores,
                    compression="gzip")
            except RuntimeError:
                # same read_id encountered previously
                pass
        return

    bc_fp = open(mh.get_megalodon_fn(out_dir, mh.BC_NAME) + '.' + bc_fmt, 'w')
    write_fastq = bc_fmt == 'fastq'
    # TODO convert this to writing un-aligned sam with htseq recommended format
    if do_output_mods:
        mods_fp = h5py.File(mh.get_megalodon_fn(out_dir, mh.BC_MODS_NAME))
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

    return


def _process_reads_worker(
        read_file_q, bc_q, vars_q, failed_reads_q, mods_q, caller_conn,
        sig_map_q, ref_out_info, model_info, vars_data, mods_info, device):
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
                LOGGER.debug('Gracefully exiting read worker {}'.format(
                    mp.current_process()))
                break
            LOGGER.debug('Analyzing read {}'.format(read_id))
            sig_info = model_info.extract_signal_info(
                fast5_fn, read_id, sig_map_q is not None)
            process_read(
                sig_info, model_info, bc_q, caller_conn, sig_map_q,
                ref_out_info, vars_data, vars_q, mods_q, mods_info,
                failed_reads_q)
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

    return


if _DO_PROFILE:
    _process_reads_wrapper = _process_reads_worker

    def _process_reads_worker(*args):
        import cProfile
        cProfile.runctx('_process_reads_wrapper(*args)', globals(), locals(),
                        filename='read_processing.prof')
        return


############################
# Post Per-read Processing #
############################

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
    aggregate.aggregate_stats(
        outputs, out_dir, num_ps, write_vcf_lp, het_factors,
        vars_data.call_mode, mods_info.agg_info,
        write_mod_lp, mods_info.mod_output_fmts, supp_prog)
    return


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

    return


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
    def update_prog(reads_called, sig_called, unexp_err_fp, read_called=True):
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
                        bar.format_dict['elapsed']})
                except AttributeError:
                    # sometimes get no format_dict error
                    # so don't include ksample/s if so
                    pass
                if q_bars is not None:
                    for q_name, q_bar in q_bars.items():
                        q_bar.n = getter_qs[q_name].queue.qsize()
                        q_bar.refresh()
                if read_called:
                    bar.update(1)
                if num_update_errors > 0:
                    bar.write(prog_prefix + format_fail_summ(
                        bar_header,
                        [(len(fns), err) for err, fns in failed_reads.items()],
                        reads_called, num_update_errors), file=sys.stderr)
            reads_called += 1

        return reads_called, unexp_err_fp

    LOGGER.info('Processing reads.')
    reads_called, sig_called = 0, 0
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
    if not suppress_progress:
        if q_bars is not None:
            while any(getter_qs[q_name].queue.qsize() > 0
                      for q_name in q_bars.keys()):
                reads_called, unexp_err_fp = update_prog(
                    reads_called, 0, unexp_err_fp, False)
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

    return


#######################
# All read processing #
#######################

def process_all_reads(
        fast5s_dir, recursive, num_reads, read_ids_fn, model_info, outputs,
        out_dir, bc_fmt, aligner, vars_data, num_ps, num_update_errors,
        suppress_progress, mods_info, db_safety, ref_out_info, do_show_qs):
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
        whatshap_map_fn = (
            mh.get_megalodon_fn(out_dir, mh.WHATSHAP_MAP_NAME) + '.' +
            aligner.out_fmt) if mh.WHATSHAP_MAP_NAME in outputs else None
        vars_txt_fn = (mh.get_megalodon_fn(out_dir, mh.PR_VAR_TXT_NAME)
                       if vars_data.write_vars_txt else None)
        getter_qs[mh.PR_VAR_NAME] = mh.create_getter_q(
            variants._get_variants_queue, (
                mh.get_megalodon_fn(out_dir, mh.PR_VAR_NAME),
                vars_txt_fn, db_safety, pr_refs_fn, ref_out_info,
                whatshap_map_fn, aligner.ref_names_and_lens, aligner.ref_fn,
                vars_data.loc_index_in_memory))
    if mh.PR_MOD_NAME in outputs:
        pr_refs_fn = mh.get_megalodon_fn(out_dir, mh.PR_REF_NAME) if (
            mh.PR_REF_NAME in outputs and ref_out_info.annotate_mods) else None
        mods_txt_fn = (mh.get_megalodon_fn(out_dir, mh.PR_MOD_TXT_NAME)
                       if mods_info.write_mods_txt else None)
        getter_qs[mh.PR_MOD_NAME] = mh.create_getter_q(
            mods._get_mods_queue, (
                mh.get_megalodon_fn(out_dir, mh.PR_MOD_NAME), db_safety,
                aligner.ref_names_and_lens, mods_txt_fn,
                pr_refs_fn, ref_out_info, mods_info.pos_index_in_memory,
                mods_info.mod_long_names))
    if mh.SIG_MAP_NAME in outputs:
        alphabet_info = signal_mapping.get_alphabet_info(model_info)
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

    proc_reads_ps, map_conns = [], []
    for device in model_info.process_devices:
        if aligner is None:
            map_conn, caller_conn = None, None
        else:
            map_conn, caller_conn = mp.Pipe()
        map_conns.append(map_conn)
        p = mp.Process(
            target=_process_reads_worker, args=(
                read_file_q, getter_qs[mh.BC_NAME].queue,
                getter_qs[mh.PR_VAR_NAME].queue, fr_prog_getter.queue,
                getter_qs[mh.PR_MOD_NAME].queue, caller_conn,
                getter_qs[mh.SIG_MAP_NAME].queue, ref_out_info, model_info,
                vars_data, mods_info, device))
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

    return


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
        args.outputs.append(mh.PR_VAR_NAME)
    if mh.WHATSHAP_MAP_NAME in args.outputs and \
       mh.VAR_NAME not in args.outputs:
        args.outputs.append(mh.VAR_NAME)
    if mh.VAR_NAME in args.outputs and mh.PR_VAR_NAME not in args.outputs:
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
    return args, vars_data


def parse_mod_args(args, model_info):
    if args.ref_include_mods and mh.PR_MOD_NAME not in args.outputs:
        args.outputs.append(mh.PR_MOD_NAME)
    if mh.PR_MOD_NAME not in args.outputs and mh.MOD_NAME in args.outputs:
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
    if args.mod_aggregate_method == mods.EM_NAME:
        agg_info = mods.AGG_INFO(mods.EM_NAME, None)
    elif args.mod_aggregate_method == mods.BIN_THRESH_NAME:
        agg_info = mods.AGG_INFO(
            mods.BIN_THRESH_NAME, args.mod_binary_threshold)
    mods_info = mods.ModInfo(
        model_info, args.mod_motif, args.mod_all_paths,
        args.write_mods_text, args.mod_context_bases,
        mh.BC_MODS_NAME in args.outputs, mod_calib_fn,
        args.mod_output_formats, args.edge_buffer,
        not args.mod_positions_on_disk, agg_info)
    return args, mods_info


def parse_ref_out_args(args, model_info):
    output_pr_refs = output_sig_maps = False
    if mh.SIG_MAP_NAME in args.outputs or mh.PR_REF_NAME in args.outputs:
        if args.ref_include_variants and args.ref_include_mods:
            LOGGER.error('Cannot output both modified base and variants in ' +
                         'per-read references (remove one of ' +
                         '--refs-include-variants or --refs-include-mods).')
            sys.exit(1)
        if args.ref_include_mods and mh.PR_MOD_NAME not in args.outputs:
            args.outputs.append(mh.PR_MOD_NAME)
            LOGGER.warning('--ref-include-mods set, so adding ' +
                           '"per_read_mods" to --outputs.')
        if args.ref_mods_all_motifs and not args.ref_include_mods:
            LOGGER.warning(
                '--ref-mods-all-motifs but not --ref-include-mods set. ' +
                'Ignoring --ref-mods-all-motifs.')
        if args.ref_mod_threshold and not args.ref_include_mods:
            LOGGER.warning(
                '--ref-mod-threshold but not --ref-include-mods set. ' +
                'Ignoring --ref-mod-threshold.')

    if mh.SIG_MAP_NAME in args.outputs:
        output_sig_maps = True
        from megalodon import signal_mapping
        global signal_mapping
        sig_map_alphabet_info = signal_mapping.get_alphabet_info(model_info)
        sig_map_alphabet = sig_map_alphabet_info.alphabet

    if mh.PR_REF_NAME in args.outputs:
        output_pr_refs = True
        if args.ref_include_variants and mh.PR_VAR_NAME not in args.outputs:
            args.outputs.append(mh.PR_VAR_NAME)
            LOGGER.warning('--refs-include-variants set, so adding ' +
                           'per_read_variants to --outputs.')
    else:
        if args.ref_include_variants:
            LOGGER.warning(
                '{0} set but {1} not requested. Ignoring {0}.'.format(
                    '--ref-include-variants', mh.PR_REF_NAME))

    if mh.SIG_MAP_NAME not in args.outputs and \
       mh.PR_REF_NAME not in args.outputs:
        sig_map_alphabet = None
        for arg_flag, arg_str in (
                (args.ref_mods_all_motifs, '--ref-mods-all-motifs'),
                (args.ref_include_mods, '--ref-include-mods'),
                (args.ref_length_range, '--ref-length-range'),
                (args.ref_mod_threshold, '--ref-mod-thresh'),
                (args.ref_percent_identity_threshold,
                 '--ref-percent-identity-threshold'),
                (args.ref_percent_coverage_threshold,
                 '--ref-percent-coverage-threshold')):
            if arg_flag:
                LOGGER.warning((
                    '{0} set but neither {1} nor {2} requested. Ignoring ' +
                    '{0}.').format(arg_str, mh.SIG_MAP_NAME, mh.PR_REF_NAME))

    min_len, max_len = (args.ref_length_range
                        if args.ref_length_range is not None else
                        (None, None))
    # set mod_thresh to infinity if all sites are to be labeled as
    # modified (--ref-mods-all-motifs)
    mod_thresh = 0.0
    if args.ref_mods_all_motifs:
        mod_thresh = np.inf
        if args.ref_mod_threshold is not None:
            LOGGER.warning(
                '--ref-mods-all-motifs and --ref-mod-threshold ' +
                'both set. Ignoring --ref-mod-threshold.')
    elif args.ref_mod_threshold is not None:
        mod_thresh = args.ref_mod_threshold

    ref_out_info = mh.REF_OUT_INFO(
        pct_idnt=args.ref_percent_identity_threshold,
        pct_cov=args.ref_percent_coverage_threshold,
        min_len=min_len, max_len=max_len, alphabet=sig_map_alphabet,
        annotate_mods=args.ref_include_mods,
        annotate_vars=args.ref_include_variants, mod_thresh=mod_thresh,
        output_pr_refs=output_pr_refs, output_sig_maps=output_sig_maps)

    return args, ref_out_info


########
# Main #
########

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

    pyg_grp = parser.add_argument_group('Guppy Backend Arguments')
    pyg_grp.add_argument(
        '--guppy-config', default=backends.DEFAULT_GUPPY_CFG,
        help='Guppy config. Default: %(default)s')
    pyg_grp.add_argument(
        '--guppy-server-path', default=backends.DEFAULT_GUPPY_SERVER_PATH,
        help='Path to guppy server executable. Default: %(default)s')
    pyg_grp.add_argument(
        '--guppy-server-port', type=int, default=backends.DEFAULT_GUPPY_PORT,
        help='Guppy server port. Default: %(default)d')

    pyg_grp.add_argument(
        '--do-not-use-guppy-server', action='store_true',
        help=hidden_help('Use alternative basecalling backend (either ' +
                         'FAST5 --post_out or taiyaki.'))
    pyg_grp.add_argument(
        '--guppy-params',
        help=hidden_help('Extra guppy server parameters. Main purpose for ' +
                         'optimal performance based on compute environment. ' +
                         'Quote parameters to avoid them being parsed by ' +
                         'megalodon.'))
    pyg_grp.add_argument(
        '--guppy-timeout', type=float, default=backends.DEFAULT_GUPPY_TIMEOUT,
        help=hidden_help('Timeout to wait for guppy server to call a single ' +
                         'read in seconds. Default: %(default)f'))
    pyg_grp.add_argument(
        '--list-supported-guppy-configs', action='store_true',
        help=hidden_help('List guppy configs with sequence variant and ' +
                         '(if applicable) modified base support.'))

    out_grp = parser.add_argument_group('Output Arguments')
    out_grp.add_argument(
        '--outputs', nargs='+',
        default=['basecalls', ], choices=tuple(mh.OUTPUT_DESCS.keys()),
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
        '--basecalls-format', choices=mh.BC_OUT_FMTS,
        default=mh.BC_OUT_FMTS[0],
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
        '--variant-calibration-filename',
        help='File containing emperical calibration for variant scores. ' +
        'See megalodon/scripts/calibrate_variant_llr_scores.py. Default: ' +
        'Load default calibration for specified guppy config.')

    var_grp.add_argument(
        '--context-min-alt-prob', type=float,
        default=mh.DEFAULT_CONTEXT_MIN_ALT_PROB,
        help=hidden_help('Minimum alternative alleles probability to ' +
                         'include variant in computation of nearby variants.' +
                         ' Default: %(default)f'))
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
        '--max-indel-size', type=int, default=mh.DEFAULT_MAX_INDEL_SIZE,
        help=hidden_help('Maximum difference in number of reference and ' +
                         'alternate bases. Default: %(default)d'))
    var_grp.add_argument(
        '--variant-all-paths', action='store_true',
        help=hidden_help('Compute forwards algorithm all paths score. ' +
                         '(Default: Viterbi best-path score)'))
    var_grp.add_argument(
        '--variants-are-atomized', action='store_true',
        help=hidden_help('Input variants have been atomized (with ' +
                         'scripts/atomize_variants.py). This saves compute ' +
                         'time, but has unpredictable behavior if variants ' +
                         'are not atomized.'))
    var_grp.add_argument(
        '--variant-context-bases', type=int, nargs=2,
        default=mh.DEFAULT_VAR_CONTEXT_BASES,
        help=hidden_help('Context bases for single base variant and indel ' +
                         'calling. Default: %(default)s'))
    var_grp.add_argument(
        '--variant-locations-on-disk', action='store_true',
        help=hidden_help('Force sequence variant locations to be stored ' +
                         'only within on disk database table. This option ' +
                         'will reduce the RAM memory requirement, but may ' +
                         'drastically slow processing. Default: Store ' +
                         'locations in memory and on disk.'))
    var_grp.add_argument(
        '--write-variants-text', action='store_true',
        help=hidden_help('Write per-read sequence variant calls out to a ' +
                         'text file. Default: Only ouput to database.'))
    var_grp.add_argument(
        '--write-vcf-log-probs', action='store_true',
        help=hidden_help('Write per-read alt log probabilities out in ' +
                         'non-standard VCF field.'))

    mod_grp = parser.add_argument_group('Modified Base Arguments')
    mod_grp.add_argument(
        '--mod-motif', action="append", nargs=3,
        metavar=('BASE', 'MOTIF', 'REL_POSITION'),
        help='Restrict modified base calls to specified motif(s). For ' +
        'example to restrict to CpG sites use "--mod-motif Z CG 0".')
    mod_grp.add_argument(
        '--mod-calibration-filename',
        help='File containing emperical calibration for modified base ' +
        'scores. See megalodon/scripts/calibrate_mod_llr_scores.py. ' +
        'Default: Load default calibration for specified guppy config.')

    mod_grp.add_argument(
        '--disable-mod-calibration', action='store_true',
        help=hidden_help('Use raw modified base scores from the network. ' +
                         'Default: Calibrate scores as described in ' +
                         '--mod-calibration-filename'))
    mod_grp.add_argument(
        '--mod-aggregate-method', choices=list(mods.AGG_METHOD_NAMES),
        default=mods.BIN_THRESH_NAME,
        help=hidden_help('Modified base aggregation method. ' +
                         'Default: %(default)s'))
    mod_grp.add_argument(
        '--mod-all-paths', action='store_true',
        help=hidden_help('Compute forwards algorithm all paths score for ' +
                         'modified base calls. (Default: Viterbi ' +
                         'best-path score)'))
    mod_grp.add_argument(
        '--mod-binary-threshold', type=float,
        default=mods.DEFAULT_BINARY_THRESH,
        help=hidden_help('Threshold for modified base aggregation ' +
                         '(probability of modified/canonical base). ' +
                         'Only applicable for "--mod-aggregate-method ' +
                         'binary_threshold". Default: %(default)s'))
    mod_grp.add_argument(
        '--mod-context-bases', type=int, default=mh.DEFAULT_MOD_CONTEXT,
        help=hidden_help('Context bases for modified base calling. ' +
                         'Default: %(default)d'))
    mod_grp.add_argument(
        '--mod-output-formats', nargs='+',
        default=[mh.MOD_BEDMETHYL_NAME, ],
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
    mod_grp.add_argument(
        '--write-mods-text', action='store_true',
        help=hidden_help('Write per-read modified bases out to a text ' +
                         'file. Default: Only ouput to database.'))

    tai_grp = parser.add_argument_group('Taiyaki Backend Arguments')
    tai_grp.add_argument(
        '--chunk-size', type=int, default=1000,
        help=hidden_help('Chunk length for base calling. ' +
                         'Default: %(default)d'))
    tai_grp.add_argument(
        '--chunk-overlap', type=int, default=100,
        help=hidden_help('Overlap between chunks to be stitched together. ' +
                         'Default: %(default)d'))
    tai_grp.add_argument(
        '--max-concurrent-chunks', type=int, default=200,
        help=hidden_help('Only process N chunks concurrently per-read (to ' +
                         'avoid GPU memory errors). Default: %(default)d'))
    tai_grp.add_argument(
        '--taiyaki-model-filename',
        help=hidden_help('Taiyaki basecalling model checkpoint file.'))

    sigmap_grp = parser.add_argument_group(
        'Reference/Signal Mapping Output Arguments')
    sigmap_grp.add_argument(
        '--ref-include-mods', action='store_true',
        help=hidden_help('Include modified base calls in signal_mappings/' +
                         'per_read_refs output.'))
    sigmap_grp.add_argument(
        '--ref-include-variants', action='store_true',
        help=hidden_help('Include variant calls in per_read_refs output ' +
                         '(does not apply to signal_mappings output).'))
    sigmap_grp.add_argument(
        '--ref-length-range', type=int, nargs=2,
        metavar=('MIN_LENGTH', 'MAX_LENGTH'),
        help=hidden_help('Only include reads with specified read length ' +
                         'in signal_mappings/per_read_refs output.'))
    sigmap_grp.add_argument(
        '--ref-percent-identity-threshold', type=float,
        help=hidden_help('Only include reads with higher percent identity ' +
                         'in signal_mappings/per_read_refs output.'))
    sigmap_grp.add_argument(
        '--ref-percent-coverage-threshold', type=float,
        help=hidden_help('Only include reads with higher read alignment ' +
                         'coverage in signal_mappings/per_read_refs output.'))
    sigmap_grp.add_argument(
        '--ref-mods-all-motifs', action='store_true',
        help=hidden_help('Annotate all --mod-motif occurences as modified.'))
    sigmap_grp.add_argument(
        '--ref-mod-threshold', type=float,
        help=hidden_help('Threshold (log(can_prob/mod_prob)) used to ' +
                         'annotate a modified bases in signal_mappings/' +
                         'per_read_refs output. See ' +
                         'scripts/compute_mod_thresh_score.py'))

    misc_grp = parser.add_argument_group('Miscellaneous Arguments')
    misc_grp.add_argument(
        '--help-long', help='Show all options.', action='help')
    misc_grp.add_argument(
        '--processes', type=int, default=1,
        help='Number of parallel processes. Default: %(default)d')
    misc_grp.add_argument(
        '--devices', nargs='+',
        help='GPU devices for guppy or taiyaki basecalling backends.')
    misc_grp.add_argument(
        '--verbose-read-progress', type=int, default=3,
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
    misc_grp.add_argument(
        '--suppress-queues-status', action='store_true',
        help=hidden_help('Suppress dynamic status of output queues. Helpful ' +
                         'for diagnosing I/O issues.'))

    return parser


def _main():
    args = get_parser().parse_args()
    if args.list_supported_guppy_configs:
        print('\n' + mh.get_supported_configs_message())
        sys.exit()

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
            args.processes, args.verbose_read_progress, args.suppress_progress,
            mods_info, args.database_safety, ref_out_info,
            not args.suppress_queues_status)

    if aligner is not None:
        ref_fn = aligner.ref_fn
        map_out_fmt = aligner.out_fmt
        del aligner

    if mh.MAP_NAME in args.outputs:
        LOGGER.info('Spawning process to sort mappings')
        map_p = post_process_mapping(
            args.output_directory, map_out_fmt, ref_fn)

    if mh.WHATSHAP_MAP_NAME in args.outputs:
        LOGGER.info('Spawning process to sort whatshap mappings')
        whatshap_sort_fn, whatshap_p = post_process_whatshap(
            args.output_directory, map_out_fmt, ref_fn)

    if mh.VAR_NAME in args.outputs or mh.MOD_NAME in args.outputs:
        post_process_aggregate(
            mods_info, args.outputs, args.output_directory, args.processes,
            args.write_vcf_log_probs, args.heterozygous_factors, vars_data,
            args.write_mod_log_probs, args.suppress_progress)

    if mh.VAR_NAME in args.outputs:
        LOGGER.info('Sorting output variant file')
        variant_fn = mh.get_megalodon_fn(args.output_directory, mh.VAR_NAME)
        sort_variant_fn = mh.add_fn_suffix(variant_fn, 'sorted')
        variants.sort_variants(variant_fn, sort_variant_fn)
        LOGGER.info('Indexing output variant file')
        index_variant_fn = variants.index_variants(sort_variant_fn)

    if mh.WHATSHAP_MAP_NAME in args.outputs:
        if whatshap_p.is_alive():
            LOGGER.info('Waiting for whatshap mappings sort')
            while whatshap_p.is_alive():
                sleep(0.001)
        LOGGER.info(variants.get_whatshap_command(
            index_variant_fn, whatshap_sort_fn,
            mh.add_fn_suffix(variant_fn, 'phased')))

    if mh.MAP_NAME in args.outputs:
        if map_p.is_alive():
            LOGGER.info('Waiting for mappings sort')
            while map_p.is_alive():
                sleep(0.001)

    return


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
