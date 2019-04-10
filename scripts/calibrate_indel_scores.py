import sys
import queue
import argparse
import threading
from time import sleep
from random import choice
import multiprocessing as mp

import mappy
import numpy as np
from tqdm import tqdm

from megalodon import (
    decode, megalodon_helper as mh, megalodon, backends, mapping, snps)


CAN_BASES = "ACGT"
CAN_BASES_SET = set(CAN_BASES)

_DO_PROFILE = False


def call_snp(
        r_post, post_mapped_start, r_snp_pos, rl_cumsum, r_to_q_poss,
        snp_ref_seq, snp_alt_seq, context_bases, all_paths,
        np_ref_seq=None, ref_seq=None):
    snp_context_bases = (context_bases[0]
                         if len(snp_ref_seq) == len(snp_alt_seq) else
                         context_bases[1])
    pos_bb = min(snp_context_bases, r_snp_pos)
    if ref_seq is None:
        pos_ab = min(snp_context_bases,
                     np_ref_seq.shape[0] - r_snp_pos - len(snp_ref_seq))
        pos_ref_seq = np_ref_seq[r_snp_pos - pos_bb:
                                 r_snp_pos + pos_ab + len(snp_ref_seq)]
    else:
        pos_ab = min(snp_context_bases,
                     len(ref_seq) - r_snp_pos - len(snp_ref_seq))
        pos_ref_seq = np.array([
            mh.ALPHABET.find(b) for b in ref_seq[
                r_snp_pos - pos_bb:
                r_snp_pos + pos_ab + len(snp_ref_seq)]], dtype=np.uintp)

    pos_alt_seq = np.concatenate([
        pos_ref_seq[:pos_bb],
        np.array([mh.ALPHABET.find(b) for b in snp_alt_seq],
                 dtype=np.uintp),
        pos_ref_seq[pos_bb + len(snp_ref_seq):]])
    blk_start  = rl_cumsum[r_to_q_poss[r_snp_pos - pos_bb]]
    blk_end = rl_cumsum[r_to_q_poss[r_snp_pos + pos_ab] + 1]

    if blk_end - blk_start < max(len(pos_ref_seq), len(pos_alt_seq)):
        return np.NAN
    loc_ref_score = snps.score_seq(
        r_post, pos_ref_seq, post_mapped_start + blk_start,
        post_mapped_start + blk_end, all_paths)
    loc_alt_score = snps.score_seq(
        r_post, pos_alt_seq, post_mapped_start + blk_start,
        post_mapped_start + blk_end, all_paths)

    return loc_ref_score - loc_alt_score

def call_alt_true_indel(
        indel_size, r_snp_pos, true_ref_seq, r_seq, map_thr_buf, context_bases,
        r_post, rl_cumsum, all_paths):
    def run_aligner():
        return next(mappy.Aligner(
            seq=false_ref_seq, preset=str('map-ont'), best_n=1).map(
                str(r_seq), buf=map_thr_buf))


    if indel_size > 0:
        # test alt truth reference insertion
        false_ref_seq = (
            true_ref_seq[:r_snp_pos + 1] +
            true_ref_seq[r_snp_pos + indel_size + 1:])
        snp_ref_seq = true_ref_seq[r_snp_pos]
        snp_alt_seq = true_ref_seq[r_snp_pos:r_snp_pos + indel_size + 1]
    else:
        # test alt truth reference deletion
        deleted_seq = ''.join(choice(CAN_BASES) for _ in range(-indel_size))
        false_ref_seq = (
            true_ref_seq[:r_snp_pos + 1] + deleted_seq +
            true_ref_seq[r_snp_pos + 1:])
        snp_ref_seq = true_ref_seq[r_snp_pos] + deleted_seq
        snp_alt_seq = true_ref_seq[r_snp_pos]

    try:
        r_algn = run_aligner()
    except StopIteration:
        raise mh.MegaError('No alignment')

    r_ref_seq = false_ref_seq[r_algn.r_st:r_algn.r_en]
    if r_algn.strand == -1:
        raise mh.MegaError('Indel mapped read mapped to reverse strand.')

    r_to_q_poss = mapping.parse_cigar(r_algn.cigar, r_algn.strand)
    if (r_algn.r_st > r_snp_pos - context_bases[1] or
        r_algn.r_en < r_snp_pos + context_bases[1]):
        raise mh.MegaError('Indel mapped read clipped snp position.')

    post_mapped_start = rl_cumsum[r_algn.q_st]
    mapped_rl_cumsum = rl_cumsum[
        r_algn.q_st:r_algn.q_en + 1] - post_mapped_start

    score = call_snp(
        r_post, post_mapped_start, r_snp_pos, rl_cumsum, r_to_q_poss,
        snp_ref_seq, snp_alt_seq, context_bases, all_paths, ref_seq=r_ref_seq)

    return score, snp_ref_seq, snp_alt_seq

def process_read(
        raw_sig, read_id, model_info, alphabet_info, caller_conn,
        map_thr_buf, context_bases=[10,30], edge_buffer=100, max_indel_len=3,
        all_paths=False, every_n=5):
    if model_info.is_cat_mod:
        bc_weights, mod_weights = model_info.run_model(
            raw_sig, n_can_state=alphabet_info.n_can_state)
    else:
        bc_weights = model_info.run_model(raw_sig)

    r_post = decode.crf_flipflop_trans_post(bc_weights, log=True)
    r_seq, score, runlen = decode.decode_post(
        r_post, alphabet_info.alphabet)
    rl_cumsum = np.cumsum(np.concatenate([[0], runlen]))

    r_ref_seq, r_to_q_poss, r_ref_pos = mapping.map_read(
        r_seq, read_id, caller_conn)
    np_ref_seq = np.array([
        mh.ALPHABET.find(b) for b in r_ref_seq], dtype=np.uintp)
    if np_ref_seq.shape[0] < edge_buffer * 2:
        raise NotImplementedError(
            'Mapping too short for calibration statistic computation.')

    read_snp_calls = []
    for r_snp_pos in range(
            edge_buffer, np_ref_seq.shape[0] - edge_buffer, every_n):
        for indel_size in range(1, max_indel_len + 1):
            try:
                score, snp_ref_seq, snp_alt_seq = call_alt_true_indel(
                    indel_size, r_snp_pos, r_ref_seq, r_seq, map_thr_buf,
                    context_bases, r_post, rl_cumsum, all_paths)
                read_snp_calls.append((False, score, snp_ref_seq, snp_alt_seq))
            except mh.MegaError:
                pass
            try:
                score, snp_ref_seq, snp_alt_seq = call_alt_true_indel(
                    -indel_size, r_snp_pos, r_ref_seq, r_seq, map_thr_buf,
                    context_bases, r_post, rl_cumsum, all_paths)
                read_snp_calls.append((False, score, snp_ref_seq, snp_alt_seq))
            except mh.MegaError:
                pass

    # get mapped start in post and run len to mapped bit of output
    post_mapped_start = sum(runlen[:r_ref_pos.q_trim_start])
    mapped_rl_cumsum = rl_cumsum[
        r_ref_pos.q_trim_start:
        r_ref_pos.q_trim_end + 1] - rl_cumsum[r_ref_pos.q_trim_start]

    for r_snp_pos in range(
            edge_buffer, np_ref_seq.shape[0] - edge_buffer, every_n):
        # test simple SNP first
        snp_ref_seq = r_ref_seq[r_snp_pos]
        snp_alt_seq = choice(list(CAN_BASES_SET.difference(snp_ref_seq)))
        score = call_snp(
            r_post, post_mapped_start, r_snp_pos, mapped_rl_cumsum,
            r_to_q_poss, snp_ref_seq, snp_alt_seq, context_bases, all_paths,
            np_ref_seq=np_ref_seq)
        read_snp_calls.append((True, score, snp_ref_seq, snp_alt_seq))

        # then test indels
        for indel_size in range(1, max_indel_len + 1):
            # test deletion
            snp_ref_seq = r_ref_seq[r_snp_pos:r_snp_pos + indel_size + 1]
            snp_alt_seq = r_ref_seq[r_snp_pos]
            score = call_snp(
                r_post, post_mapped_start, r_snp_pos, mapped_rl_cumsum,
                r_to_q_poss, snp_ref_seq, snp_alt_seq, context_bases,
                all_paths, np_ref_seq=np_ref_seq)
            read_snp_calls.append((True, score, snp_ref_seq, snp_alt_seq))

            # test random insertion
            snp_ref_seq = r_ref_seq[r_snp_pos]
            snp_alt_seq = snp_ref_seq + ''.join(
                choice(CAN_BASES) for _ in range(indel_size))
            score = call_snp(
                r_post, post_mapped_start, r_snp_pos, mapped_rl_cumsum,
                r_to_q_poss, snp_ref_seq, snp_alt_seq, context_bases,
                all_paths, np_ref_seq=np_ref_seq)
            read_snp_calls.append((True, score, snp_ref_seq, snp_alt_seq))

    return read_snp_calls

def _process_reads_worker(
        fast5_q, snp_calls_q, caller_conn, model_info, alphabet_info, device):
    model_info.prep_model_worker(device)
    map_thr_buf = mappy.ThreadBuffer()

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
            read_snp_calls = process_read(
                raw_sig, read_id, model_info, alphabet_info, caller_conn,
                map_thr_buf)
            snp_calls_q.put(read_snp_calls)
        except:
            #raise
            pass

    return

if _DO_PROFILE:
    _process_reads_wrapper = _process_reads_worker
    def _process_reads_worker(*args):
        import cProfile
        cProfile.runctx('_process_reads_wrapper(*args)', globals(), locals(),
                        filename='snp_calibration.prof')
        return


def _get_snp_calls(
        snp_calls_q, snp_calls_conn, out_fn, total_reads, suppress_progress):
    out_fp = open(out_fn, 'w')
    if not suppress_progress:
        bar = tqdm(total=total_reads, smoothing=0)

    while True:
        try:
            read_snp_calls = snp_calls_q.get(block=False)
            for snp_call in read_snp_calls:
                out_fp.write('{}\t{}\t{}\t{}\n'.format(*snp_call))
            out_fp.flush()
            if not suppress_progress:
                bar.update(1)
        except queue.Empty:
            if snp_calls_conn.poll():
                break
            sleep(0.1)
            continue

    while not snp_calls_q.empty():
        read_snp_calls = snp_calls_q.get(block=False)
        for snp_call in read_snp_calls:
            out_fp.write('{}\t{}\t{}\t{}\n'.format(*snp_call))
        out_fp.flush()
        if not suppress_progress:
            bar.update(1)
    out_fp.close()
    if not suppress_progress:
        bar.close()

    return


def process_all_reads(
        fast5s_dir, num_reads, model_info, aligner, num_ps, alphabet_info,
        out_fn, suppress_progress):
    sys.stderr.write('Searching for reads.\n')
    fast5_fns = megalodon.get_read_files(fast5s_dir)
    if num_reads is not None:
        fast5_fns = fast5_fns[:num_reads]


    sys.stderr.write('Preparing workers and calling reads.\n')
    # read filename queue filler
    fast5_q = mp.Queue(maxsize=mh._MAX_QUEUE_SIZE)
    files_p = mp.Process(
        target=megalodon._fill_files_queue, args=(fast5_q, fast5_fns, num_ps),
        daemon=True)
    files_p.start()

    snp_calls_q, snp_calls_p, main_sc_conn = mh.create_getter_q(
            _get_snp_calls, (out_fn, len(fast5_fns), suppress_progress))

    proc_reads_ps, map_conns = [], []
    for device in model_info.process_devices:
        if aligner is None:
            map_conn, caller_conn = None, None
        else:
            map_conn, caller_conn = mp.Pipe()
        map_conns.append(map_conn)
        p = mp.Process(
            target=_process_reads_worker, args=(
                fast5_q, snp_calls_q, caller_conn, model_info, alphabet_info,
                device))
        p.daemon = True
        p.start()
        proc_reads_ps.append(p)
    sleep(0.1)
    map_read_ts = []
    for map_conn in map_conns:
        t = threading.Thread(
            target=mapping._map_read_worker,
            args=(aligner, map_conn, None, False))
        t.daemon = True
        t.start()
        map_read_ts.append(t)

    files_p.join()
    for proc_reads_p in proc_reads_ps:
        proc_reads_p.join()
    if map_read_ts is not None:
        for map_t in map_read_ts:
            map_t.join()
    if snp_calls_p.is_alive():
        main_sc_conn.send(True)
        snp_calls_p.join()

    return


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

    map_grp = parser.add_argument_group('Mapping Arguments')
    map_grp.add_argument(
        '--reference',
        help='Reference FASTA file used for mapping called reads.')

    out_grp = parser.add_argument_group('Output Arguments')
    out_grp.add_argument(
        '--output', default='snp_calibration_statistics.txt',
        help='Filename to output statistics. Default: %(default)s')
    out_grp.add_argument(
        '--num-reads', type=int,
        help='Number of reads to process. Default: All reads')

    tai_grp = parser.add_argument_group('Taiyaki Signal Chunking Arguments')
    tai_grp.add_argument(
        '--devices', type=int, nargs='+',
        help='CUDA GPU devices to use (only valid for taiyaki), default: CPU')
    tai_grp.add_argument(
        '--chunk_size', type=int, default=1000,
        help='Chunk length for base calling. Default: %(default)d')
    tai_grp.add_argument(
        '--chunk_overlap', type=int, default=100,
        help='Overlap between chunks to be stitched together. ' +
        'Default: %(default)d')
    tai_grp.add_argument(
        '--max_concurrent_chunks', type=int, default=50,
        help='Only process N chunks concurrently per-read (to avoid GPU ' +
        'memory errors). Default: %(default)d')

    misc_grp = parser.add_argument_group('Miscellaneous Arguments')
    misc_grp.add_argument(
        '--processes', type=int, default=1,
        help='Number of parallel processes. Default: %(default)d')
    misc_grp.add_argument(
        '--suppress-progress', action='store_true',
        help='Suppress progress bar.')

    return parser

def main():
    args = get_parser().parse_args()

    model_info = backends.ModelInfo(
        args.flappie_model_name, args.taiyaki_model_filename, args.devices,
        args.processes, args.chunk_size, args.chunk_overlap,
        args.max_concurrent_chunks)
    alphabet_info = megalodon.AlphabetInfo(
        model_info, None, False, None, None)
    aligner = mapping.alignerPlus(
        str(args.reference), preset=str('map-ont'), best_n=1)

    process_all_reads(
        args.fast5s_dir, args.num_reads, model_info, aligner, args.processes,
        alphabet_info, args.output, args.suppress_progress)

    return

if __name__ == '__main__':
    main()
