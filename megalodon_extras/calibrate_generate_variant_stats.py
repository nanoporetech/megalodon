import queue
import threading
from time import sleep
from random import choice
import multiprocessing as mp
from collections import defaultdict

import mappy
import numpy as np
from tqdm import tqdm

from megalodon import (
    backends, fast5_io, logging, mapping, megalodon_helper as mh,
    megalodon_multiprocessing as mega_mp, variants)
from ._extras_parsers import get_parser_calibrate_generate_variants_stats


MAX_INDEL_LEN = 5
ALL_PATHS = False
TEST_EVERY_N_LOCS = 5
MAX_POS_PER_READ = 400

CAN_BASES = "ACGT"
CAN_BASES_SET = set(CAN_BASES)

LOGGER = logging.get_logger()

_DO_PROFILE = False


def call_variant(
        can_post, post_mapped_start, r_var_pos, rl_cumsum, r_to_q_poss,
        var_ref_seq, var_alt_seq, context_bases, all_paths,
        np_ref_seq=None, ref_seq=None):
    var_context_bases = (context_bases[0]
                         if len(var_ref_seq) == len(var_alt_seq) else
                         context_bases[1])
    pos_bb = min(var_context_bases, r_var_pos)
    if ref_seq is None:
        pos_ab = min(var_context_bases,
                     np_ref_seq.shape[0] - r_var_pos - len(var_ref_seq))
        pos_ref_seq = np_ref_seq[r_var_pos - pos_bb:
                                 r_var_pos + pos_ab + len(var_ref_seq)]
    else:
        pos_ab = min(var_context_bases,
                     len(ref_seq) - r_var_pos - len(var_ref_seq))
        pos_ref_seq = mh.seq_to_int(ref_seq[
            r_var_pos - pos_bb:r_var_pos + pos_ab + len(var_ref_seq)])

    pos_alt_seq = np.concatenate([
        pos_ref_seq[:pos_bb], mh.seq_to_int(var_alt_seq),
        pos_ref_seq[pos_bb + len(var_ref_seq):]])
    blk_start = rl_cumsum[r_to_q_poss[r_var_pos - pos_bb]]
    blk_end = rl_cumsum[r_to_q_poss[r_var_pos + pos_ab] + 1]

    if blk_end - blk_start < max(len(pos_ref_seq), len(pos_alt_seq)):
        return np.NAN
    loc_ref_score = variants.score_seq(
        can_post, pos_ref_seq, post_mapped_start + blk_start,
        post_mapped_start + blk_end, all_paths)
    loc_alt_score = variants.score_seq(
        can_post, pos_alt_seq, post_mapped_start + blk_start,
        post_mapped_start + blk_end, all_paths)

    return loc_ref_score - loc_alt_score


def call_alt_true_indel(
        indel_size, r_var_pos, true_ref_seq, r_seq, map_thr_buf, context_bases,
        can_post, rl_cumsum, all_paths):
    def run_aligner():
        return next(mappy.Aligner(
            seq=false_ref_seq, preset=str('map-ont'), best_n=1).map(
                str(r_seq), buf=map_thr_buf))

    if indel_size == 0:
        false_base = choice(
            list(set(CAN_BASES).difference(true_ref_seq[r_var_pos])))
        false_ref_seq = (
            true_ref_seq[:r_var_pos] + false_base +
            true_ref_seq[r_var_pos + 1:])
        var_ref_seq = false_base
        var_alt_seq = true_ref_seq[r_var_pos]
    elif indel_size > 0:
        # test alt truth reference insertion
        false_ref_seq = (
            true_ref_seq[:r_var_pos + 1] +
            true_ref_seq[r_var_pos + indel_size + 1:])
        var_ref_seq = true_ref_seq[r_var_pos]
        var_alt_seq = true_ref_seq[r_var_pos:r_var_pos + indel_size + 1]
    else:
        # test alt truth reference deletion
        deleted_seq = ''.join(choice(CAN_BASES) for _ in range(-indel_size))
        false_ref_seq = (
            true_ref_seq[:r_var_pos + 1] + deleted_seq +
            true_ref_seq[r_var_pos + 1:])
        var_ref_seq = true_ref_seq[r_var_pos] + deleted_seq
        var_alt_seq = true_ref_seq[r_var_pos]

    try:
        r_algn = run_aligner()
    except StopIteration:
        raise mh.MegaError('No alignment')

    r_ref_seq = false_ref_seq[r_algn.r_st:r_algn.r_en]
    if r_algn.strand == -1:
        raise mh.MegaError('Indel mapped read mapped to reverse strand.')

    r_to_q_poss = mapping.parse_cigar(r_algn.cigar, r_algn.strand)
    if r_algn.r_st > r_var_pos - context_bases[1] or \
       r_algn.r_en < r_var_pos + context_bases[1]:
        raise mh.MegaError('Indel mapped read clipped variant position.')

    post_mapped_start = rl_cumsum[r_algn.q_st]
    mapped_rl_cumsum = rl_cumsum[
        r_algn.q_st:r_algn.q_en + 1] - post_mapped_start

    score = call_variant(
        can_post, post_mapped_start, r_var_pos, mapped_rl_cumsum, r_to_q_poss,
        var_ref_seq, var_alt_seq, context_bases, all_paths, ref_seq=r_ref_seq)

    return score, var_ref_seq, var_alt_seq


def process_read(
        bc_res, caller_conn, map_thr_buf, do_false_ref,
        context_bases=mh.DEFAULT_VAR_CONTEXT_BASES,
        edge_buffer=mh.DEFAULT_EDGE_BUFFER, max_indel_len=MAX_INDEL_LEN,
        all_paths=ALL_PATHS, every_n=TEST_EVERY_N_LOCS,
        max_pos_per_read=MAX_POS_PER_READ):
    sig_info, called_read, rl_cumsum, can_post, = bc_res
    r_ref_seq, r_to_q_poss, r_ref_pos, _ = mapping.map_read(
        caller_conn, called_read, backends.SIGNAL_DATA(*sig_info))
    np_ref_seq = mh.seq_to_int(r_ref_seq)
    if np_ref_seq.shape[0] < edge_buffer * 2:
        raise NotImplementedError(
            'Mapping too short for calibration statistic computation.')
    # get mapped start in post and run len to mapped bit of output
    post_mapped_start = rl_cumsum[r_ref_pos.q_trim_start]
    mapped_rl_cumsum = rl_cumsum[
        r_ref_pos.q_trim_start:r_ref_pos.q_trim_end + 1] - post_mapped_start

    # candidate variant locations within a read
    var_poss = list(range(
        edge_buffer, np_ref_seq.shape[0] - edge_buffer,
        every_n))[:max_pos_per_read]
    read_var_calls = []

    if do_false_ref:
        # first process reference false calls (need to spoof an incorrect
        # reference for mapping and signal remapping)
        for r_var_pos in var_poss:
            # first test single base swap SNPs
            try:
                score, var_ref_seq, var_alt_seq = call_alt_true_indel(
                    0, r_var_pos, r_ref_seq, called_read.seq, map_thr_buf,
                    context_bases, can_post, rl_cumsum, all_paths)
                read_var_calls.append((False, score, var_ref_seq, var_alt_seq))
            except mh.MegaError:
                # introduced error either causes read not to map,
                # mapping trims the location of interest or invalid ref or alt
                # sequence
                pass
            # then test small indels
            for indel_size in range(1, max_indel_len + 1):
                try:
                    score, var_ref_seq, var_alt_seq = call_alt_true_indel(
                        indel_size, r_var_pos, r_ref_seq, called_read.seq,
                        map_thr_buf, context_bases, can_post, rl_cumsum,
                        all_paths)
                    read_var_calls.append((
                        False, score, var_ref_seq, var_alt_seq))
                except mh.MegaError:
                    pass
                try:
                    score, var_ref_seq, var_alt_seq = call_alt_true_indel(
                        -indel_size, r_var_pos, r_ref_seq, called_read.seq,
                        map_thr_buf, context_bases, can_post, rl_cumsum,
                        all_paths)
                    read_var_calls.append((
                        False, score, var_ref_seq, var_alt_seq))
                except mh.MegaError:
                    pass

    # now test reference correct variants
    for r_var_pos in var_poss:
        if len(set(r_ref_seq[
                r_var_pos - max(context_bases):
                r_var_pos + max_indel_len + 1 +
                max(context_bases)]).difference(CAN_BASES_SET)) > 0:
            # skip reference positions with N's in any context
            continue
        # test simple SNP first
        var_ref_seq = r_ref_seq[r_var_pos]
        for var_alt_seq in CAN_BASES_SET.difference(var_ref_seq):
            try:
                score = call_variant(
                    can_post, post_mapped_start, r_var_pos, mapped_rl_cumsum,
                    r_to_q_poss, var_ref_seq, var_alt_seq, context_bases,
                    all_paths, np_ref_seq=np_ref_seq)
            except mh.MegaError:
                # invalid reference or alternative sequence
                continue
            read_var_calls.append((True, score, var_ref_seq, var_alt_seq))

        # then test indels
        for indel_size in range(1, max_indel_len + 1):
            # test deletion
            var_ref_seq = r_ref_seq[r_var_pos:r_var_pos + indel_size + 1]
            var_alt_seq = r_ref_seq[r_var_pos]
            try:
                score = call_variant(
                    can_post, post_mapped_start, r_var_pos, mapped_rl_cumsum,
                    r_to_q_poss, var_ref_seq, var_alt_seq, context_bases,
                    all_paths, np_ref_seq=np_ref_seq)
            except mh.MegaError:
                # invalid reference or alternative sequence
                continue
            read_var_calls.append((True, score, var_ref_seq, var_alt_seq))

            # test random insertion
            var_ref_seq = r_ref_seq[r_var_pos]
            var_alt_seq = var_ref_seq + ''.join(
                choice(CAN_BASES) for _ in range(indel_size))
            try:
                score = call_variant(
                    can_post, post_mapped_start, r_var_pos, mapped_rl_cumsum,
                    r_to_q_poss, var_ref_seq, var_alt_seq, context_bases,
                    all_paths, np_ref_seq=np_ref_seq)
            except mh.MegaError:
                # invalid reference or alternative sequence
                continue
            read_var_calls.append((True, score, var_ref_seq, var_alt_seq))

    return read_var_calls


def _process_reads_worker(bc_q, var_calls_q, caller_conn, do_false_ref):
    LOGGER.debug('InitWorker')
    map_thr_buf = mappy.ThreadBuffer()

    while True:
        try:
            bc_res = bc_q.get(block=True, timeout=0.1)
        except queue.Empty:
            continue
        if bc_res is None:
            LOGGER.debug('Closing')
            if caller_conn is not None:
                caller_conn.close()
            break

        try:
            read_var_calls = process_read(
                bc_res, caller_conn, map_thr_buf, do_false_ref)
            var_calls_q.put((True, read_var_calls))
        except Exception as e:
            var_calls_q.put((False, str(e)))


if _DO_PROFILE:
    _process_reads_wrapper = _process_reads_worker

    def _process_reads_worker(*args):
        import cProfile
        cProfile.runctx('_process_reads_wrapper(*args)', globals(), locals(),
                        filename='variant_calibration.prof')


def basecall_worker(
        sig_q, bc_q, model_info, device,
        reads_per_batch=mh.DEFAULT_GUPPY_BATCH_SIZE):
    def create_batch_gen():
        for _ in range(reads_per_batch):
            try:
                read_sig_data = sig_q.get(timeout=0.01)
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
            for (sig_info, _, called_read, rl_cumsum,
                 can_post, _, _) in model_info.iter_basecalled_reads(
                     reads_batch_gen, return_post_w_mods=False):
                yield tuple(sig_info), called_read, rl_cumsum, can_post

    LOGGER.debug('InitBasecaller')
    model_info.prep_model_worker(device)
    full_iter_conn, gen_conn = mp.Pipe(duplex=False)
    for bc_res in iter_bc_res():
        bc_q.put(bc_res)


def extract_signal_worker(fn_read_ids_q, sig_q, model_info):
    while True:
        try:
            fn_rids = fn_read_ids_q.get(block=True, timeout=0.1)
        except queue.Empty:
            continue
        if fn_rids is None:
            LOGGER.debug('Closing')
            break

        fast5_fn, read_ids = fn_rids
        with fast5_io.get_fast5_file(fast5_fn) as fast5_fp:
            for read_id in read_ids:
                sig_info, seq_summ_info = model_info.extract_signal_info(
                    fast5_fp, read_id)
                sig_q.put((tuple(sig_info), tuple(seq_summ_info)))


def _get_variant_calls(
        var_calls_q, var_calls_conn, out_fn, getter_num_reads_conn,
        suppress_progress):
    out_fp = open(out_fn, 'w')
    bar = None
    if not suppress_progress:
        bar = tqdm(smoothing=0, dynamic_ncols=True)

    err_types = defaultdict(int)
    while True:
        try:
            valid_res, read_var_calls = var_calls_q.get(block=False)
            if valid_res:
                for var_call in read_var_calls:
                    out_fp.write('{}\t{}\t{}\t{}\n'.format(*var_call))
            else:
                err_types[read_var_calls] += 1
            if not suppress_progress:
                bar.update(1)
        except queue.Empty:
            if bar is not None and bar.total is None:
                if getter_num_reads_conn.poll():
                    bar.total = getter_num_reads_conn.recv()
            else:
                if var_calls_conn.poll():
                    break
            sleep(0.01)
            continue

    while not var_calls_q.empty():
        valid_res, read_var_calls = var_calls_q.get(block=False)
        if valid_res:
            for var_call in read_var_calls:
                out_fp.write('{}\t{}\t{}\t{}\n'.format(*var_call))
        else:
            err_types[read_var_calls] += 1
        if not suppress_progress:
            bar.update(1)
    out_fp.close()
    if not suppress_progress:
        bar.close()

    if len(err_types) > 0:
        LOGGER.info('Failed reads summary:')
        for n_errs, err_str in sorted(
                (v, k) for k, v in err_types.items())[::-1]:
            LOGGER.info('\t{} : {} reads'.format(err_str, n_errs))


def process_all_reads(
        fast5s_dir, recursive, num_reads, read_ids_fn, model_info, aligner,
        num_ps, out_fn, suppress_progress, do_false_ref):
    LOGGER.info('Preparing workers and calling reads.')
    num_extract_sig_ps = num_ps
    num_bc_ps = num_ps
    # read filename queue filler
    fn_read_ids_q = mega_mp.CountingMPQueue()
    sig_q = mega_mp.CountingMPQueue(maxsize=mh._MAX_QUEUE_SIZE)
    bc_q = mega_mp.CountingMPQueue(maxsize=mh._MAX_QUEUE_SIZE)
    num_reads_conn, getter_num_reads_conn = mp.Pipe()
    input_info = mh.INPUT_INFO(
        fast5s_dir=fast5s_dir, recursive=recursive, num_reads=num_reads,
        read_ids_fn=read_ids_fn, num_ps=num_ps, do_it_live=False)
    aux_failed_q = mp.Queue()
    files_p = mp.Process(
        target=fast5_io._fill_files_queue, args=(
            input_info, fn_read_ids_q, num_reads_conn, aux_failed_q),
        daemon=True, name='FileFiller')
    files_p.start()
    extract_sig_ps = []
    for espi in range(num_extract_sig_ps):
        extract_sig_ps.append(mp.Process(
            target=extract_signal_worker,
            args=(fn_read_ids_q, sig_q, model_info),
            daemon=True, name='ExtractSig{:03d}'.format(espi)))
        extract_sig_ps[-1].start()
    bc_ps = []
    for bcpi, device in enumerate(model_info.process_devices):
        bc_ps.append(mp.Process(
            target=basecall_worker, args=(sig_q, bc_q, model_info, device),
            daemon=True, name='Basecaller{:03d}'.format(bcpi)))
        bc_ps[-1].start()

    var_calls_q, var_calls_p, main_sc_conn = mega_mp.create_getter_qpc(
        _get_variant_calls,
        (out_fn, getter_num_reads_conn, suppress_progress))

    proc_reads_ps, map_conns = [], []
    for pnum, device in enumerate(model_info.process_devices):
        if aligner is None:
            map_conn, caller_conn = None, None
        else:
            map_conn, caller_conn = mp.Pipe()
        map_conns.append(map_conn)
        p = mp.Process(
            target=_process_reads_worker,
            args=(bc_q, var_calls_q, caller_conn, do_false_ref),
            name='ReadWorker{:03d}'.format(pnum))
        p.daemon = True
        p.start()
        proc_reads_ps.append(p)
        if caller_conn is not None:
            caller_conn.close()
        del caller_conn
    sleep(0.1)
    map_read_ts = []
    for map_conn in map_conns:
        t = threading.Thread(target=mapping._map_read_worker,
                             args=(aligner, map_conn))
        t.daemon = True
        t.start()
        map_read_ts.append(t)

    LOGGER.debug('Waiting for read filler process to join.')
    files_p.join()
    LOGGER.debug('Read filler process complete.')
    # put extra None values in fn_read_ids_q to indicate completed processing
    for _ in range(num_extract_sig_ps):
        fn_read_ids_q.put(None)
    LOGGER.debug('Waiting for signal extraction processes to join.')
    for extract_sig_p in extract_sig_ps:
        extract_sig_p.join()
    LOGGER.debug('Signal extraction complete.')
    # put None values in sig_q to indicate completed processing
    for _ in range(num_bc_ps):
        sig_q.put(None)
    LOGGER.debug('Waiting for basecall processes to join.')
    for bc_p in bc_ps:
        bc_p.join()
    LOGGER.debug('Basecalling complete.')
    # put None values in bc_q to indicate completed processing
    for _ in range(num_ps):
        bc_q.put(None)
    LOGGER.debug('Waiting for process reads processes to join.')
    for proc_reads_p in proc_reads_ps:
        proc_reads_p.join()
    LOGGER.debug('Process reads processes complete.')
    LOGGER.debug('Waiting for mapping threads to join.')
    if map_read_ts is not None:
        for map_t in map_read_ts:
            map_t.join()
    LOGGER.debug('Mapping threads complete.')
    LOGGER.debug('Waiting for output process to join.')
    if var_calls_p.is_alive():
        main_sc_conn.send(True)
        var_calls_p.join()
    LOGGER.debug('Output process complete.')


def _main(args):
    try:
        mh.mkdir(args.guppy_logs_output_directory, False)
    except mh.MegaError:
        LOGGER.warning(
            'Guppy logs output directory exists. Potentially overwriting ' +
            'guppy logs.')
    logging.init_logger(args.guppy_logs_output_directory)
    # add required attributes for loading guppy, but not valid options for
    # this script.
    args.do_not_use_guppy_server = False
    args.output_directory = args.guppy_logs_output_directory

    LOGGER.info('Loading model.')
    backend_params = backends.parse_backend_params(args)
    with backends.ModelInfo(backend_params, args.processes) as model_info:
        LOGGER.info('Loading reference.')
        aligner = mappy.Aligner(
            str(args.reference), preset=str('map-ont'), best_n=1)

        process_all_reads(
            args.fast5s_dir, not args.not_recursive, args.num_reads,
            args.read_ids_filename, model_info, aligner, args.processes,
            args.output, args.suppress_progress,
            args.compute_false_reference_scores)


if __name__ == '__main__':
    _main(get_parser_calibrate_generate_variants_stats().parse_args())
