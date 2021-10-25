import queue
from time import sleep
from random import choice
import multiprocessing as mp

import numpy as np
from tqdm import tqdm

from taiyaki import mapped_signal_files
from megalodon import backends, logging, megalodon_helper as mh, mods
from ._extras_parsers import get_parser_calibrate_generate_mod_stats_from_msf


LOGGER = logging.get_logger()


def output_mods_data(all_mod_llrs, all_can_llrs, out_fn):
    LOGGER.info("Merging modified base data")
    all_mod_bases = list(
        set(all_mod_llrs.keys()).intersection(all_can_llrs.keys())
    )
    if len(set(all_mod_llrs.keys()).difference(all_mod_bases)) > 0:
        LOGGER.warning(
            "Modified base(s) found in modified dataset which were not "
            + "found in canonical dataset: {}".format(
                ",".join(set(all_mod_llrs.keys()).difference(all_mod_bases))
            )
        )
    if len(set(all_can_llrs.keys()).difference(all_mod_bases)) > 0:
        LOGGER.warning(
            "Modified base(s) found in modified dataset which were "
            + "not found in canonical dataset: {}".format(
                ",".join(set(all_mod_llrs.keys()).difference(all_mod_bases))
            )
        )
    mod_base_stats = {mods.GT_ALL_MOD_BASE_STR: all_mod_bases}
    for mod_base in all_mod_bases:
        mod_base_stats[mods.GT_MOD_LLR_STR.format(mod_base)] = all_mod_llrs[
            mod_base
        ]
        mod_base_stats[mods.GT_CAN_LLR_STR.format(mod_base)] = all_can_llrs[
            mod_base
        ]
    np.savez(out_fn, **mod_base_stats)


def iter_motifs(motifs, alphabet, can_labs, int_ref_seq, context_bases):
    for motif, rel_pos in motifs:
        for motif_match in motif.finditer(
            "".join(alphabet[bi] for bi in can_labs[int_ref_seq])
        ):
            ref_pos = motif_match.start() + rel_pos
            if context_bases <= ref_pos < int_ref_seq.shape[0] - context_bases:
                yield ref_pos


def process_read(
    read,
    motifs,
    can_labs,
    mod_labs,
    can_post_indices,
    model_info,
    edge_buffer,
    context_bases,
):
    sig_info = backends.SIGNAL_DATA(
        dacs=read.Dacs,
        raw_len=read.Dacs.shape[0],
        fast5_fn="",
        read_id=read.read_id,
        stride=model_info.stride,
        channel_info={
            mh.CHAN_INFO_OFFSET: read.offset,
            mh.CHAN_INFO_RANGE: read.range,
            mh.CHAN_INFO_DIGI: read.digitisation,
        },
    )
    seq_summ_info = mh.SEQ_SUMM_INFO("", read.read_id, "", "", "", "", 0, 0)
    post_w_mods = next(
        model_info.iter_basecalled_reads([(sig_info, seq_summ_info)])
    )[5]
    # convert read.Ref_to_signal to blocks coordinates with
    # model_info.stride
    ref_block_pos = np.around(read.Ref_to_signal / model_info.stride).astype(
        int
    )
    # loop over motif hits in read.Reference
    # then run mods.score_mod_seq on extracted locations
    for ref_pos in iter_motifs(
        motifs, model_info.can_alphabet, can_labs, read.Reference, context_bases
    ):
        seq_st, seq_en = ref_pos - context_bases, ref_pos + context_bases + 1
        pos_seq = read.Reference[seq_st:seq_en]
        post_st, post_en = ref_block_pos[seq_st], ref_block_pos[seq_en]
        if post_en - post_st < pos_seq.shape[0]:
            continue
        can_seq, mod_seq = can_labs[pos_seq], mod_labs[pos_seq]
        alt_mod_seq = mod_seq.copy()
        if mod_labs[read.Reference[ref_pos]] == 0:
            can_base = model_info.output_alphabet[read.Reference[ref_pos]]
            mod_base = choice(model_info.can_base_mods[can_base])
            alt_mod_lab = model_info.str_to_int_mod_labels[mod_base]
            mod_base = model_info.output_alphabet[
                int(read.Reference[ref_pos]) + alt_mod_lab
            ]
            alt_mod_seq[context_bases] = alt_mod_lab
            mod_score = mods.score_mod_seq(
                post_w_mods,
                can_seq,
                mod_seq,
                can_post_indices,
                post_st,
                post_en,
            ) - mods.score_mod_seq(
                post_w_mods,
                can_seq,
                alt_mod_seq,
                can_post_indices,
                post_st,
                post_en,
            )
            yield True, mod_base, mod_score
        else:
            mod_base = model_info.output_alphabet[read.Reference[ref_pos]]
            alt_mod_seq[context_bases] = 0
            mod_score = mods.score_mod_seq(
                post_w_mods,
                can_seq,
                alt_mod_seq,
                can_post_indices,
                post_st,
                post_en,
            ) - mods.score_mod_seq(
                post_w_mods,
                can_seq,
                mod_seq,
                can_post_indices,
                post_st,
                post_en,
            )
            yield False, mod_base, mod_score


def _process_reads_worker(
    read_q,
    mod_calls_q,
    motifs,
    can_labs,
    mod_labs,
    can_post_indices,
    model_info,
    device,
    edge_buffer,
    context_bases,
):
    model_info.prep_model_worker(device)

    while True:
        try:
            read = read_q.get(timeout=0.1)
        except queue.Empty:
            continue
        if read is None:
            break

        read_mod_scores = []
        for can_truth, mod_base, mod_score in process_read(
            read,
            motifs,
            can_labs,
            mod_labs,
            can_post_indices,
            model_info,
            edge_buffer,
            context_bases,
        ):
            read_mod_scores.append((not can_truth, mod_base, mod_score))
        mod_calls_q.put(read_mod_scores)


def fill_reads_queue(
    read_q, read_filler_conn, ms_fn, num_reads_limit, num_proc
):
    msf = mapped_signal_files.HDF5Reader(ms_fn)
    num_reads = 0
    for read in msf:
        read_q.put(read)
        num_reads += 1
        if num_reads_limit is not None and num_reads >= num_reads_limit:
            break
    read_filler_conn.send(num_reads)
    msf.close()
    for _ in num_proc:
        read_q.put(None)


def compute_diff_scores(
    ms_fn,
    model_info,
    context_bases,
    edge_buffer,
    num_reads_limit,
    motifs,
    can_labs,
    mod_labs,
    can_post_indices,
):
    LOGGER.info("Processing reads")
    read_q = mp.Queue(maxsize=mh._MAX_QUEUE_SIZE)
    num_reads_conn, read_filler_conn = mp.Pipe(duplex=False)
    read_fill_p = mp.Process(
        target=fill_reads_queue,
        args=(
            read_q,
            read_filler_conn,
            ms_fn,
            num_reads_limit,
            model_info.process_devices,
        ),
        daemon=True,
    )
    read_fill_p.start()

    mod_calls_q = mp.Queue()
    LOGGER.info("Processing mapped signal reads.")
    proc_reads_ps = []
    for device in model_info.process_devices:
        p = mp.Process(
            target=_process_reads_worker,
            args=(
                read_q,
                mod_calls_q,
                motifs,
                can_labs,
                mod_labs,
                can_post_indices,
                model_info,
                device,
                edge_buffer,
                context_bases,
            ),
        )
        p.daemon = True
        p.start()
        proc_reads_ps.append(p)

    bar = tqdm(smoothing=0, dynamic_ncols=True, unit="reads")
    all_mod_llrs = dict(
        (mod_base, []) for mod_base in model_info.mod_base_to_can
    )
    all_can_llrs = dict(
        (mod_base, []) for mod_base in model_info.mod_base_to_can
    )
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


def extract_label_conversions(model_info):
    can_labels, mod_labels = [0], [0]
    curr_can_lab = 0
    for base in model_info.output_alphabet[1:]:
        if base in model_info.can_alphabet:
            curr_can_lab += 1
            can_labels.append(curr_can_lab)
            mod_labels.append(0)
        else:
            can_labels.append(curr_can_lab)
            mod_labels.append(model_info.str_to_int_mod_labels[base])

    return (
        np.array(can_labels, dtype=np.uintp),
        np.array(mod_labels, dtype=np.uintp),
    )


def parse_motifs(raw_motifs, model_info, mod_set):
    if mod_set is not None:
        mod_set = set(mod_set)
    can_bases = set()
    for mod_base, can_base in model_info.mod_base_to_can.items():
        if mod_set is not None and mod_base not in mod_set:
            continue
        can_bases.add(can_base)
    if raw_motifs is None:
        return [(mh.compile_motif_pat(can_base), 0) for can_base in can_bases]
    motifs = []
    for raw_motif, rel_pos in raw_motifs:
        rel_pos = int(rel_pos)
        if raw_motif[rel_pos] not in can_bases:
            LOGGER.warning(
                (
                    'Motif "{}" not associated with requested modified '
                    + "bases."
                ).format(raw_motif)
            )
            continue
        motifs.append((mh.compile_motif_pat(raw_motif), rel_pos))
    return motifs


def check_map_sig_alphabet(model_info, ms_fn):
    # read filename queue filler
    msf = mapped_signal_files.HDF5Reader(ms_fn)
    tai_alph_info = msf.get_alphabet_information()
    msf.close()
    if model_info.output_alphabet != tai_alph_info.alphabet:
        raise mh.MegaError(
            (
                "Different alphabets specified in model ({}) and mapped "
                + "signal file ({})"
            ).format(model_info.output_alphabet, tai_alph_info.alphabet)
        )
    if set(model_info.can_alphabet) != set(tai_alph_info.collapse_alphabet):
        raise mh.MegaError(
            (
                "Different canonical alphabets specified in model ({}) and "
                + "mapped signal file ({})"
            ).format(model_info.can_alphabet, tai_alph_info.collapse_alphabet)
        )
    if model_info.ordered_mod_long_names != tai_alph_info.mod_long_names:
        raise mh.MegaError(
            (
                "Different modified base long names specified in model ({}) and "
                + "mapped signal file ({})"
            ).format(
                ", ".join(model_info.ordered_mod_long_names),
                ", ".join(tai_alph_info.mod_long_names),
            )
        )


def add_trim_guppy_none(args):
    if args.guppy_params is None:
        args.guppy_params = "--trim_strategy none"
    else:
        args.guppy_params += " --trim_strategy none"
    return args


def _main(args):
    logging.init_logger(log_fn=args.log_filename, quiet=args.quiet)
    # add required attributes for loading guppy, but not valid options for
    # this script.
    args.do_not_use_guppy_server = False
    args.output_directory = args.guppy_logs_output_directory
    try:
        mh.mkdir(args.output_directory, False)
    except mh.MegaError:
        LOGGER.warning(
            "Guppy logs output directory exists. Potentially overwriting "
            + "guppy logs."
        )
    args = add_trim_guppy_none(args)
    args.outputs = [mh.PR_MOD_NAME]
    # make edge_buffer >= context_bases to simplify processing
    if args.edge_buffer < args.mod_context_bases:
        LOGGER.warning(
            "[--edge-buffer] less than [--mod-context-bases]. Setting "
            + "[--edge-buffer] to value from [--mod-context-bases]"
        )
        args.edge_buffer = args.mod_context_bases

    LOGGER.info("Loading model.")
    backend_params = backends.parse_backend_params(args)
    with backends.ModelInfo(backend_params, args.processes) as model_info:
        check_map_sig_alphabet(model_info, args.mapped_signal_file)
        motifs = parse_motifs(args.motif, model_info, args.modified_bases_set)
        can_labs, mod_labs = extract_label_conversions(model_info)
        can_post_indices = model_info.can_indices.astype(np.uintp)
        all_mod_llrs, all_can_llrs = compute_diff_scores(
            args.mapped_signal_file,
            model_info,
            args.mod_context_bases,
            args.edge_buffer,
            args.num_reads,
            motifs,
            can_labs,
            mod_labs,
            can_post_indices,
        )

    mod_summary = [
        (
            mod,
            len(all_mod_llrs[mod]) if mod in all_mod_llrs else 0,
            len(all_can_llrs[mod]) if mod in all_can_llrs else 0,
        )
        for mod in set(all_mod_llrs).union(all_can_llrs)
    ]
    LOGGER.info(
        "Data summary:\n\tmod\tmod_N\tcan_N\n"
        + "\n".join("\t" + "\t".join(map(str, x)) for x in mod_summary)
    )
    output_mods_data(all_mod_llrs, all_can_llrs, args.out_filename)


if __name__ == "__main__":
    _main(get_parser_calibrate_generate_mod_stats_from_msf().parse_args())
