from random import choice

import mappy
import numpy as np

from megalodon import decode, megalodon_helper as mh, megalodon, backends


CAN_BASES = "ACGT"
CAN_BASES_SET = set(CAN_BASES)


def call_snp(
        r_post, post_mapped_start, r_snp_pos, r_ref_seq, rl_cumsum, r_to_q_poss,
        snp_ref_seq, snp_alt_seq, context_bases, all_paths):
    snp_context_bases = (context_bases[0]
                         if len(snp_ref_seq) == len(snp_alt_seq) else
                         context_bases[1])
    pos_bb = min(snp_context_bases, r_snp_pos)
    pos_ab = min(snp_context_bases,
                 r_ref_seq.shape[0] - r_snp_pos - len(snp_ref_seq))
    pos_ref_seq = r_ref_seq[r_snp_pos - pos_bb:
                            r_snp_pos + pos_ab + len(snp_ref_seq)]
    pos_alt_seq = np.concatenate([
        pos_ref_seq[:pos_bb],
        np.array([mh.ALPHABET.find(b) for b in snp_alt_seq],
                 dtype=np.uintp),
        pos_ref_seq[pos_bb + len(snp_ref_seq):]])
    blk_start  = rl_cumsum[r_to_q_poss[r_snp_pos - pos_bb]]
    blk_end = rl_cumsum[r_to_q_poss[r_snp_pos + pos_ab]]

    if blk_end - blk_start < max(len(pos_ref_seq), len(pos_alt_seq)):
        return np.NAN
    loc_ref_score = score_seq(
        r_post, pos_ref_seq, post_mapped_start + blk_start,
        post_mapped_start + blk_end, all_paths)
    loc_alt_score = score_seq(
        r_post, pos_alt_seq, post_mapped_start + blk_start,
        post_mapped_start + blk_end, all_paths)

    return loc_ref_score - loc_alt_score

def process_read(
        raw_sig, read_id, model_info, caller_conn, alphabet_info,
        context_bases=[10,30], edge_buffer=100, max_indel_len=3,
        all_paths=False):
    if model_info.is_cat_mod:
        bc_weights, mod_weights = model_info.run_model(
            raw_sig, alphabet_info.n_can_state)
    else:
        bc_weights = model_info.run_model(raw_sig)

    r_post = decode.crf_flipflop_trans_post(bc_weights, log=True)
    r_seq, score, runlen = decode.decode_post(
        r_post, alphabet_info.collapse_alphabet)

    r_ref_seq, r_to_q_poss, r_ref_pos = mapping.map_read(
        r_seq, read_id, caller_conn)

    # TODO need to re-map to new reference in order to test ground truth
    # involves creating new mappy aligner

    np_ref_seq = np.array([
        mh.ALPHABET.find(b) for b in r_ref_seq], dtype=np.uintp)
    # get mapped start in post and run len to mapped bit of output
    post_mapped_start = sum(runlen[:r_ref_pos.q_trim_start])
    rl_cumsum = np.cumsum(np.concatenate([
        [0], runlen[r_ref_pos.q_trim_start:r_ref_pos.q_trim_end]]))

    if np_ref_seq.shape[0] < edge_buffer * 2:
        raise NotImplementedError('Too few reads.')
    for r_snp_pos in range(edge_buffer, np_ref_seq.shape[0]):
        # test simple SNP first
        snp_ref_seq = r_ref_seq[r_snp_pos]
        snp_alt_seq = choice(CAN_BASES.difference(snp_ref_seq))
        score = call_snp(
            r_post, post_mapped_start, r_snp_pos, r_ref_seq, rl_cumsum,
            r_to_q_poss, snp_ref_seq,
            snp_alt_seq, context_bases, all_paths)
        # TODO change from printing to passing back through a queue
        print('TRUE\t{}\t{}\t{}\n'.format(score, snp_ref_seq, snp_alt_seq))

        # then test indels
        for indel_size in range(1, max_indel_len + 1):
            # test deletion
            snp_ref_seq = r_ref_seq[r_snp_pos:r_snp_pos + indel_size]
            snp_alt_seq = r_ref_seq[r_snp_pos]
            score = call_snp(
                r_post, post_mapped_start, r_snp_pos, r_ref_seq, rl_cumsum,
                r_to_q_poss, snp_ref_seq, snp_alt_seq, context_bases, all_paths)
            print('TRUE\t{}\t{}\t{}\n'.format(score, snp_ref_seq, snp_alt_seq))

            # test random insertion
            snp_ref_seq = r_ref_seq[r_snp_pos]
            snp_alt_seq = snp_ref_seq + ''.join(
                choice(CAN_BASES) for _ in range(indel_size))
            score = call_snp(
                r_post, post_mapped_start, r_snp_pos, r_ref_seq, rl_cumsum,
                r_to_q_poss, snp_ref_seq, snp_alt_seq, context_bases, all_paths)
            print('TRUE\t{}\t{}\t{}\n'.format(score, snp_ref_seq, snp_alt_seq))

    return

def _process_reads_worker(fast5_q, caller_conn, model_info, alphabet_info):
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
                raw_sig, read_id, model_info, caller_conn, alphabet_info)
        except:
            pass

    return


def process_all_reads(
        fast5s_dir, num_reads, model_info, aligner, num_ps, alphabet_info):
    sys.stderr.write('Searching for reads.\n')
    fast5_fns = megalodon.get_read_files(fast5s_dir)
    if num_reads is not None:
        fast5_fns = fast5_fns[:num_reads]


    sys.stderr.write('Preparing workers and calling reads.\n')
    # read filename queue filler
    fast5_q = mp.Queue(maxsize=mh._MAX_QUEUE_SIZE)
    files_p = mp.Process(
        target=_fill_files_queue, args=(fast5_q, fast5_fns, num_ps),
        daemon=True)
    files_p.start()


    proc_reads_ps, map_conns = [], []
    for _ in range(num_ps):
        if aligner is None:
            map_conn, caller_conn = None, None
        else:
            map_conn, caller_conn = mp.Pipe()
        map_conns.append(map_conn)
        p = mp.Process(
            target=_process_reads_worker, args=(
                fast5_q, caller_conn, model_info, alphabet_info))
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
        '--num-reads', type=int,
        help='Number of reads to process. Default: All reads')

    misc_grp = parser.add_argument_group('Miscellaneous Arguments')
    misc_grp.add_argument(
        '--processes', type=int, default=1,
        help='Number of parallel processes. Default: %(default)d')
    misc_grp.add_argument(
        '--device', default=None, type=int,
        help='CUDA device to use (only valid for taiyaki), or None to use CPU')

    return parser

def main():
    args = get_parser().parse_args()

    model_info = backends.ModelInfo(
        args.flappie_model_name, args.taiyaki_model_filename, args.device)
    alphabet_info = megalodon.AlphabetInfo(
        model_info, None, False, None, None, None)
    aligner = mapping.alignerPlus(
        str(args.reference), preset=str('map-ont'), best_n=1)

    process_all_reads(
        args.fast5s_dir, args.num_reads, model_info, aligner, args.processes,
        alphabet_info)

    return

if __name__ == '__main__':
    main()
