import sys
import queue
from time import sleep
from collections import namedtuple

import numpy as np

from ont_fast5_api import fast5_interface
from megalodon import megalodon_helper as mh, logging
from taiyaki import (
    alphabet, fast5utils, mapping as tai_mapping, prepare_mapping_funcs,
    signal as tai_signal)


LOGGER = logging.get_logger()
SIG_MAP_RESULT = namedtuple('SIG_MAP_RESULT', (
    'pass_filts', 'fast5_fn', 'dacs', 'scale_params', 'ref_seq', 'stride',
    'read_id', 'r_to_q_poss', 'rl_cumsum', 'ref_pos', 'ref_out_info'))


def get_remapping(
        sig_fn, dacs, scale_params, ref_seq, stride, read_id, r_to_q_poss,
        rl_cumsum, r_ref_pos, ref_out_info):
    read = fast5_interface.get_fast5_file(sig_fn, 'r').get_read(read_id)
    channel_info = dict(fast5utils.get_channel_info(read).items())
    rd_factor = channel_info['range'] / channel_info['digitisation']
    shift_from_pA = (scale_params[0] + channel_info['offset']) * rd_factor
    scale_from_pA = scale_params[1] * rd_factor
    read_attrs = dict(fast5utils.get_read_attributes(read).items())

    # prepare taiyaki signal object
    sig = tai_signal.Signal(dacs=dacs)
    sig.channel_info = channel_info
    sig.read_attributes = read_attrs
    sig.offset = channel_info['offset']
    sig.range = channel_info['range']
    sig.digitisation = channel_info['digitisation']

    path = np.full((dacs.shape[0] // stride) + 1, -1)
    # skip last value since this is where the two seqs end
    for ref_pos, q_pos in enumerate(r_to_q_poss[:-1]):
        # if the query position maps to the end of the mapping skip it
        if rl_cumsum[q_pos + r_ref_pos.q_trim_start] >= path.shape[0]:
            continue
        path[rl_cumsum[q_pos + r_ref_pos.q_trim_start]] = ref_pos
    remapping = tai_mapping.Mapping.from_remapping_path(
        sig, path, ref_seq, stride)
    try:
        remapping.add_integer_reference(ref_out_info.alphabet)
    except Exception:
        raise mh.MegaError('Invalid reference sequence encountered')

    return (remapping.get_read_dictionary(
        shift_from_pA, scale_from_pA, read_id),
            prepare_mapping_funcs.RemapResult.SUCCESS)


def get_alphabet_info(model_info):
    flat_alphabet = model_info.output_alphabet[0]
    can_base = model_info.output_alphabet[0]
    for base in model_info.output_alphabet[1:]:
        if base in model_info.can_alphabet:
            can_base = base
        flat_alphabet += can_base
    mod_long_names = [] if len(model_info.mod_long_names) == 0 else \
        list(zip(*model_info.mod_long_names))[1]
    return alphabet.AlphabetInfo(
        model_info.output_alphabet, flat_alphabet,
        mod_long_names, do_reorder=True)


def write_signal_mappings(sig_map_q, sig_map_conn, sig_map_fn, alphabet_info):
    def iter_mappings():
        while True:
            try:
                read_mapping = sig_map_q.get(block=False)
                yield read_mapping
            except queue.Empty:
                if sig_map_conn.poll():
                    break
                sleep(0.001)
                continue

        while not sig_map_q.empty():
            read_mapping = sig_map_q.get(block=False)
            yield read_mapping

    prepare_mapping_funcs.generate_output_from_results(
        iter_mappings(), sig_map_fn, alphabet_info, verbose=False)

    return


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
