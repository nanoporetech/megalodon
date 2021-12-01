import sys
import queue
import traceback
from collections import namedtuple

import numpy as np

from megalodon import fast5_io, megalodon_helper as mh, logging

try:
    import taiyaki

    # appease flake8
    taiyaki
except ImportError:
    raise mh.MegaError(
        "Taiyaki installation required for signal mapping not found."
    )
try:
    from taiyaki import (
        alphabet,
        fast5utils,
        signal_mapping as tai_mapping,
        prepare_mapping_funcs,
        signal as tai_signal,
    )
except ImportError:
    raise mh.MegaError(
        "Taiyaki modules could not be loaded. Signal mappings require "
        + 'Taiyaki version >= 5.2. Full error:\n"""\n{}\n"""'.format(
            traceback.format_exc()
        )
    )


LOGGER = logging.get_logger()
SIG_MAP_RESULT = namedtuple(
    "SIG_MAP_RESULT",
    (
        "pass_filts",
        "fast5_fn",
        "dacs",
        "scale_params",
        "ref_seq",
        "stride",
        "read_id",
        "r_to_q_poss",
        "rl_cumsum",
        "ref_pos",
        "ref_out_info",
    ),
)


def set_all_motif_mods(int_ref, ref_mods_all_motifs):
    ref_mod_pos, ref_mods = [], []
    for mod_base, int_mod_base, mln, int_motif, rel_pos in ref_mods_all_motifs:
        for pos in np.where(
            np.all(
                mh.rolling_window(int_ref, int_motif.shape[0]) == int_motif,
                axis=1,
            )
        )[0]:
            ref_mod_pos.append(pos + rel_pos)
            ref_mods.append(int_mod_base)
    if len(ref_mod_pos) > 0:
        int_ref[ref_mod_pos] = ref_mods
    return int_ref


def get_remapping(
    sig_fn,
    dacs,
    scale_params,
    ref_seq,
    stride,
    read_id,
    r_to_q_poss,
    rl_cumsum,
    r_ref_pos,
    ref_out_info,
):
    read = fast5_io.get_read(sig_fn, read_id)
    channel_info = dict(fast5utils.get_channel_info(read).items())
    read_params = {
        "trim_start": 0,
        "trim_end": 0,
        "shift": scale_params[0],
        "scale": scale_params[1],
    }
    sig = tai_signal.Signal(
        dacs=dacs,
        channel_info=channel_info,
        read_id=read_id,
        read_params=read_params,
    )

    ref_to_sig = np.empty(len(ref_seq) + 1, dtype=np.int32)
    # skip last value since this is where the two seqs end
    for ref_pos, q_pos in enumerate(r_to_q_poss):
        ref_to_sig[ref_pos] = rl_cumsum[q_pos + r_ref_pos.q_trim_start] * stride
    try:
        int_ref = tai_mapping.SignalMapping.get_integer_reference(
            ref_seq, ref_out_info.alphabet_info.alphabet
        )
    except Exception:
        raise mh.MegaError("Invalid reference sequence encountered")
    sig_mapping = tai_mapping.SignalMapping(ref_to_sig, int_ref, signalObj=sig)

    # annotate mod motifs
    if ref_out_info.ref_mods_all_motifs is not None:
        # annotate all mod base motif positions with alts
        int_ref = set_all_motif_mods(int_ref, ref_out_info.ref_mods_all_motifs)
        # set new Reference with mods annotated
        sig_mapping.Reference = int_ref

    return (
        sig_mapping.get_read_dictionary(),
        prepare_mapping_funcs.RemapResult.SUCCESS,
    )


def get_alphabet_info_from_model(model_info):
    flat_alphabet = model_info.output_alphabet[0]
    can_base = model_info.output_alphabet[0]
    for base in model_info.output_alphabet[1:]:
        if base in model_info.can_alphabet:
            can_base = base
        flat_alphabet += can_base
    mod_long_names = (
        []
        if len(model_info.mod_long_names) == 0
        else list(zip(*model_info.mod_long_names))[1]
    )
    return alphabet.AlphabetInfo(
        model_info.output_alphabet,
        flat_alphabet,
        mod_long_names,
        do_reorder=True,
    )


def get_alphabet_info(output_alphabet, collapse_alphabet, mod_long_names):
    return alphabet.AlphabetInfo(
        output_alphabet, collapse_alphabet, mod_long_names, do_reorder=True
    )


def write_signal_mappings(sig_map_q, sig_map_conn, ref_out_info, aux_failed_q):
    def apply_sig_map_offset(read_mapping):
        """Apply signal mapping shift to center coarse mappings to a registered
        signal based mapping.
        """
        if (
            ref_out_info.sig_map_offset is not None
            and ref_out_info.sig_map_offset != 0
        ):
            if ref_out_info.sig_map_offset > 0:
                # clip beginning of signal mapping and end of reference to
                # shift signal assignments to the left
                read_mapping[0]["Ref_to_signal"] = read_mapping[0][
                    "Ref_to_signal"
                ][ref_out_info.sig_map_offset :]
                read_mapping[0]["Reference"] = read_mapping[0]["Reference"][
                    : -ref_out_info.sig_map_offset
                ]
            else:
                # clip end of signal mapping and beginning of reference to
                # shift signal assignments to the right
                read_mapping[0]["Ref_to_signal"] = read_mapping[0][
                    "Ref_to_signal"
                ][: ref_out_info.sig_map_offset]
                read_mapping[0]["Reference"] = read_mapping[0]["Reference"][
                    -ref_out_info.sig_map_offset :
                ]
        return read_mapping

    def iter_mappings():
        workers_active = True
        LOGGER.debug("GetterInitComplete")
        while workers_active or not sig_map_q.empty():
            try:
                read_mapping = sig_map_q.get(timeout=0.1)
                yield apply_sig_map_offset(read_mapping)
            except queue.Empty:
                if sig_map_conn.poll():
                    workers_active = False

    try:
        LOGGER.debug("GetterStarting")
        prepare_mapping_funcs.generate_output_from_results(
            iter_mappings(),
            mh.get_megalodon_fn(ref_out_info.out_dir, mh.SIG_MAP_NAME),
            ref_out_info.alphabet_info,
            verbose=False,
        )
        LOGGER.debug("GetterClosing")
    except Exception as e:
        aux_failed_q.put(
            ("SigMapProcessingError", str(e), traceback.format_exc())
        )


if __name__ == "__main__":
    sys.stderr.write("This is a module. See commands with `megalodon -h`")
    sys.exit(1)
