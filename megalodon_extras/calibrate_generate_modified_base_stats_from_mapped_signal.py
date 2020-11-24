import queue
from time import sleep
import multiprocessing as mp

import numpy as np
from tqdm import tqdm

from taiyaki import mapped_signal_files
from megalodon import backends, logging, megalodon_helper as mh, mods
from ._extras_parsers import \
    get_parser_calibrate_generate_mod_diff_ctc_stats


LOGGER = logging.get_logger()


def output_mods_data(
        all_mod_llrs, all_can_llrs, mod_base_set, exclude_mod_bases, out_fn):
    LOGGER.info('Merging modified base data')
    all_mod_bases = list(set(all_mod_llrs.keys()).intersection(
        all_can_llrs.keys()))
    if len(set(all_mod_llrs.keys()).difference(all_mod_bases)) > 0:
        LOGGER.warning(
            'Modified base(s) found in modified dataset which were not ' +
            'found in canonical dataset: {}'.format(','.join(
                set(all_mod_llrs.keys()).difference(all_mod_bases))))
    if len(set(all_can_llrs.keys()).difference(all_mod_bases)) > 0:
        LOGGER.warning(
            'Modified base(s) found in modified dataset which were ' +
            'not found in canonical dataset: {}'.format(','.join(
                set(all_mod_llrs.keys()).difference(all_mod_bases))))
    if mod_base_set is not None:
        all_mod_bases = list(set(all_mod_bases).intersection(mod_base_set))
        if len(all_mod_bases) == 0:
            LOGGER.error((
                'No modified bases to process.\n\tModified bases from ' +
                'results: {}\n\tModified base set: {}').format(
                    ','.join(all_mod_bases), ','.join(mod_base_set)))
    if exclude_mod_bases is not None:
        all_mod_bases = list(set(all_mod_bases).difference(exclude_mod_bases))
        if len(all_mod_bases) == 0:
            LOGGER.error((
                'No modified bases to process.\n\tModified bases from ' +
                'results: {}\n\tExcluded modified bases: {}').format(
                    ','.join(all_mod_bases), ','.join(exclude_mod_bases)))
    mod_base_stats = {mods.GT_ALL_MOD_BASE_STR: all_mod_bases}
    for mod_base in all_mod_bases:
        mod_base_stats[mods.GT_MOD_LLR_STR.format(
            mod_base)] = all_mod_llrs[mod_base]
        mod_base_stats[mods.GT_CAN_LLR_STR.format(
            mod_base)] = all_can_llrs[mod_base]
    np.savez(out_fn, **mod_base_stats)


def process_read(read, model_info, edge_buffer, context_bases):
    sig_info = backends.SIGNAL_DATA(
        dacs=read.Dacs, raw_len=read.Dacs.shape[0], fast5_fn='',
        read_id=read.read_id, stride=model_info.stride)
    post_w_mods = model_info.basecall_read(
        sig_info, return_post_w_mods=True)[5]
    # TODO convert read.Ref_to_signal to blocks coordinates with
    # model_info.stride
    # TOOD loop over motif hits in read.Reference
    # then run mods.score_mod_seq on extracted locations
    for ref_pos in iter_motifs():
        pose_st, post_en = get_block_positions()
        if post_en - post_st < can_seq.shape[0]:
            continue
        mod_score = (
            mods.score_mod_seq(
                mapped_post, can_seq, ref_cats, can_post_indices,
                post_st, post_en) -
            mods.score_mod_seq(
                mapped_post, can_seq, mod_cats, can_post_indices,
                post_st, post_en))
        yield ref_base in can_base_indices, mod_base, mod_score


def _process_reads_worker(
        read_q, mod_calls_q, model_info, device, edge_buffer, context_bases):
    model_info.prep_model_worker(device)

    while True:
        try:
            read = read_q.get(block=False)
        except queue.Empty:
            sleep(0.001)
            continue
        if read is None:
            break

        read_mod_scores = []
        for can_truth, mod_base, mod_score in process_read(
                read, model_info, edge_buffer, context_bases):
            read_mod_scores.append((not can_truth, mod_base, mod_score))
        mod_calls_q.put(read_mod_scores)


def fill_reads_queue(
        read_q, read_filler_conn, ms_fn, num_reads_limit, num_proc):
    msf = mapped_signal_files.HDF5Reader(ms_fn)
    num_reads = 0
    for read in msf:
        read_q.put(read)
        num_reads += 1
        if num_reads_limit is not None and num_reads >= num_reads_limit:
            break
    read_filler_conn.send(num_reads)
    for _ in num_proc:
        read_q.put(None)
    msf.close()


def check_map_sig_alphabet(model_info, ms_fn):
    # read filename queue filler
    msf = mapped_signal_files.HDF5Reader(ms_fn)
    tai_alph_info = msf.get_alphabet_information()
    msf.close()
    if model_info.output_alphabet != tai_alph_info.alphabet:
        raise mh.MegaError((
            'Different alphabets specified in model ({}) and mapped ' +
            'signal file ({})').format(model_info.output_alphabet,
                                       tai_alph_info.alphabet))
    if set(model_info.can_alphabet) != set(tai_alph_info.collapse_alphabet):
        raise mh.MegaError((
            'Different canonical alphabets specified in model ({}) and ' +
            'mapped signal file ({})').format(
                model_info.can_alphabet,
                tai_alph_info.collapse_alphabet))
    if model_info.mod_long_names != tai_alph_info.mod_long_names:
        raise mh.MegaError((
            'Different modified base long names specified in model ({}) and ' +
            'mapped signal file ({})').format(
                ', '.join(model_info.mod_long_names),
                ', '.join(tai_alph_info.mod_long_names)))


def compute_diff_scores(
        ms_fn, model_info, context_bases, edge_buffer, num_reads_limit):
    # make edge_buffer >= context_bases to simplify processing
    edge_buffer = max(context_bases, edge_buffer)

    LOGGER.info('Processing reads')
    check_map_sig_alphabet(model_info, ms_fn)
    # TODO parse motifs using model_info into integer encoded motifs

    read_q = mp.Queue(maxsize=mh._MAX_QUEUE_SIZE)
    read_filler_conn, num_reads_conn = mp.Pipe(duplex=False)
    read_fill_p = mp.Process(
        target=fill_reads_queue, args=(
            read_q, read_filler_conn, ms_fn, num_reads_limit,
            model_info.process_devices), daemon=True)
    read_fill_p.start()

    mod_calls_q = mp.Queue()
    LOGGER.info('Processing mapped signal reads.')
    proc_reads_ps = []
    for device in model_info.process_devices:
        p = mp.Process(
            target=_process_reads_worker, args=(
                read_q, mod_calls_q, model_info, device, edge_buffer,
                context_bases))
        p.daemon = True
        p.start()
        proc_reads_ps.append(p)

    bar = tqdm(smoothing=0, dynamic_ncols=True, unit='reads')
    all_mod_llrs = dict((mod_base, [])
                        for mod_base in model_info.mod_base_to_can)
    all_can_llrs = dict((mod_base, [])
                        for mod_base in model_info.mod_base_to_can)
    while any(p.is_alive() for p in proc_reads_ps):
        try:
            read_calls = mod_calls_q.get(block=False, timeout=0.1)
            for is_mod, mod_base, mod_score in read_calls:
                if is_mod:
                    all_mod_llrs[mod_base].append(mod_score)
                else:
                    all_can_llrs[mod_base].append(mod_score)
            bar.update(1)
            if bar.total is None:
                if num_reads_conn.poll():
                    bar.total = num_reads_conn.recv()
        except queue.Empty:
            continue

    sleep(0.1)
    while not mod_calls_q.empty():
        read_calls = mod_calls_q.get(block=False)
        for is_mod, mod_base, mod_score in read_calls:
            if is_mod:
                all_mod_llrs[mod_base].append(mod_score)
            else:
                all_can_llrs[mod_base].append(mod_score)
        bar.update(1)
    bar.close()

    return dict(all_mod_llrs), dict(all_can_llrs)


def _main(args):
    logging.init_logger(quiet=args.quiet)
    # add required attributes for loading guppy, but not valid options for
    # this script.
    args.do_not_use_guppy_server = False
    args.output_directory = args.guppy_logs_output_directory
    try:
        mh.mkdir(args.output_directory, False)
    except mh.MegaError:
        LOGGER.warning(
            'Guppy logs output directory exists. Potentially overwriting ' +
            'guppy logs.')
    # TODO add option for motif

    LOGGER.info('Loading model.')
    backend_params = backends.parse_backend_params(args)
    # TODO add no trim option to guppy to align caller output
    with backends.ModelInfo(backend_params, args.processes) as model_info:
        all_mod_llrs, all_can_llrs = compute_diff_scores(
            args.mapped_signal_file, model_info, args.mod_context_bases,
            args.edge_buffer, args.num_reads)

    mod_summary = [
        (mod, len(all_mod_llrs[mod]) if mod in all_mod_llrs else 0,
         len(all_can_llrs[mod]) if mod in all_can_llrs else 0)
        for mod in set(all_mod_llrs).union(all_can_llrs)]
    LOGGER.info(
        'Data summary:\n\tmod\tmod_N\tcan_N\n' + '\n'.join(
            '\t' + '\t'.join(map(str, x)) for x in mod_summary))
    output_mods_data(
        all_mod_llrs, all_can_llrs, args.modified_bases_set,
        args.exclude_modified_bases, args.out_filename)


if __name__ == '__main__':
    _main(get_parser_calibrate_generate_mod_diff_ctc_stats().parse_args())
