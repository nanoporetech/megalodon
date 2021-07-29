import numpy as np

from megalodon import banding, decode, logging, megalodon_helper as mh

LOGGER = logging.get_logger()


def construct_allowed_bases(seq):
    """Construct numpy allowed bases array from string sequence containing
    ambiguous bases.
    """
    allowed_bases = np.zeros((len(seq), 4), dtype=np.short)
    for seq_pos, base in enumerate(seq):
        try:
            pos_allowed_bases = mh.SINGLE_LETTER_CODE[base]
        except KeyError:
            raise mh.MegaError(
                f"Invalid IUPAC code ({base}) found wheile performing "
                "constrained basecalling."
            )
        for pos_allowed_base in pos_allowed_bases:
            allowed_bases[seq_pos, mh.ALPHABET.find(pos_allowed_base)] = 1
    return allowed_bases


def constrained_basecall(
    reference,
    trans_logprobs,
    ref_to_block,
    half_bandwidth=mh.DEFAULT_CONSTRAINED_HALF_BW,
):
    """Perform constrained basecalling from initial sequence to ambiguous
    reference.

    Args:
        reference (str): Reference sequence containing ambiguous bases
        trans_logprob (np.array): 2D Float array containing flip-flop transition
            log probabilities. Shape should be num_blocks by num_transitions.
            num_blocks is signal // stride and num_transitions is the number of
            flip-flop transitions (40 for 4 canonical bases).
        ref_to_block (np.array): Containing initial path coordinates from
            reference bases to block coordinates in trans_logprob
        half_bandwidth (int): Half bandwidth over which to restrict path
            between sequence and blocks. Band will be constructed along
            block/signal dimension.
    """
    # if initial mapping starts within trans_logprobs trim and shift mapping
    if ref_to_block[0] != 0:
        trans_logprobs = trans_logprobs[ref_to_block[0] :]
        ref_to_block = ref_to_block - ref_to_block[0]
    # if mapping ends before end of trans_logprobs trim
    if ref_to_block.shape[0] > ref_to_block[-1]:
        trans_logprobs = trans_logprobs[: ref_to_block[-1]]
    allowed_bases = construct_allowed_bases(reference)
    sig_band = banding.compute_sig_band(
        ref_to_block, np.zeros(len(reference)), half_bandwidth
    )
    seq_band = banding.convert_to_seq_band(sig_band)
    int_seq = decode.flipflop_constrain_decode(
        trans_logprobs, allowed_bases, seq_band
    )
    constrained_seq = "".join(mh.ALPHABET[base % 4] for base in int_seq)
    return constrained_seq
