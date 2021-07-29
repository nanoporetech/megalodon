import numpy as np

from megalodon import megalodon_helper as mh, logging

LOGGER = logging.get_logger()


def compute_sig_band(bps, levels, bhw=mh.DEFAULT_CONSTRAINED_HALF_BW):
    """Compute band over which to explore possible paths. Band is represented
    in sequence/level coordinates at each signal position.

    Args:
        bps (np.ndarray): Integer array containing breakpoints
        levels (np.ndarray): float array containing expected signal levels. May
            contain np.NAN values. Band will be constructed to maintain path
            through NAN regions.
        bhw (int): Band half width. If None, full matrix is used.

    Returns:
        int32 np.ndarray with shape (2, sig_len = bps[-1] - bps[0]). The first
        row contains the lower band boundaries in sequence coordinates and the
        second row contains the upper boundaries in sequence coordinates.
    """
    seq_len = levels.shape[0]
    if bps.shape[0] - 1 != seq_len:
        raise mh.MegaError("Breakpoints must be one longer than levels.")
    sig_len = bps[-1] - bps[0]
    seq_indices = np.repeat(np.arange(seq_len), np.diff(bps))

    # Calculate bands
    # The 1st row consists of the start indices (inc) and the 2nd row
    # consists of the end indices (exc) of the valid rows for each col.
    band = np.empty((2, sig_len), dtype=np.int32)
    if bhw is None:
        # specify entire input matrix
        band[0, :] = 0
        band[1, :] = seq_len
    else:
        # use specific band defined by bhw
        band[0, :] = np.maximum(seq_indices - bhw, 0)
        band[1, :] = np.minimum(seq_indices + bhw + 1, seq_len)

    # Modify bands based on invalid levels
    nan_mask = np.isin(seq_indices, np.nonzero(np.isnan(levels)))
    nan_sig_indices = np.where(nan_mask)[0]
    nan_seq_indices = seq_indices[nan_mask]
    band[0, nan_sig_indices] = nan_seq_indices
    band[1, nan_sig_indices] = nan_seq_indices + 1
    # Modify bands close to invalid levels so monotonically increasing
    band[0, :] = np.maximum.accumulate(band[0, :])
    band[1, :] = np.minimum.accumulate(band[1, ::-1])[::-1]

    # expand band around large deletions to ensure valid paths
    invalid_indices = np.where(band[0, 1:] >= band[1, :-1])[0]
    while invalid_indices.shape[0] > 0:
        band[0, invalid_indices + 1] = np.maximum(
            band[0, invalid_indices + 1] - 1, 0
        )
        band[1, invalid_indices] = np.minimum(
            band[1, invalid_indices] + 1, seq_len
        )
        invalid_indices = np.where(band[0, 1:] >= band[1, :-1])[0]

    return band


def convert_to_seq_band(sig_band):
    """Convert band with sig_len entries containing upper and lower band
    boundaries in base coordinates to a seq_len entries contraining upper and
    lower band boundaries in signal space.

    Args:
        sig_band (np.array): int32 array with shape (2, sig_len). The first row
            contains the lower band boundaries in sequence coordinates and the
            second row contains the upper boundaries in sequence coordinates.

    Returns:
        int32 np.ndarray with shape (2, seq_len = sig_band[1, -1]). The first
        row contains the lower band boundaries in signal coordinates and the
        second row contains the upper boundaries in signal coordinates.
    """
    sig_len = sig_band.shape[1]
    seq_len = sig_band[1, -1]
    seq_band = np.zeros((2, seq_len), dtype=np.int32)
    seq_band[1, :] = sig_len

    # upper signal coordinates define lower sequence boundaries
    lower_sig_pos = np.nonzero(np.ediff1d(sig_band[1, :], to_begin=0))[0]
    lower_base_pos = sig_band[1, lower_sig_pos - 1]
    seq_band[0, lower_base_pos] = lower_sig_pos
    seq_band[0, :] = np.maximum.accumulate(seq_band[0, :])

    upper_sig_pos = np.nonzero(np.ediff1d(sig_band[0, :], to_begin=0))[0]
    upper_base_pos = sig_band[0, upper_sig_pos]
    seq_band[1, upper_base_pos - 1] = upper_sig_pos
    seq_band[1, :] = np.minimum.accumulate(seq_band[1, ::-1])[::-1]

    return seq_band


def validate_band(band, sig_len=None, seq_len=None, is_sig_band=True):
    """Validate that band is valid and agrees with input data.

    Args:
        band (np.array): int32 array with shape (2, sig_len or seq_len). The
            first row contains the lower band boundaries and the second row
            contains the upper boundaries.
        sig_len (int): Length of signal associated with band
        seq_len (int): Length of sequence/levels associated with band
        is_sig_band (bool): Does the provided band specify sequence/level
            positions for each signal position? If not it is assumed that the
            band contains signal positions for each sequence/level position.

    Raises:
        MegaError if any portion of the band is determined to be invalid.
    """
    # first coordinate 0, last coordinate signal length
    if band[0, 0] != 0:
        raise mh.MegaError("Band does not start with 0 coordinate.")

    # ends all greater than starts
    if np.diff(band, axis=0)[0].min() <= 0:
        raise mh.MegaError("Band contains 0-length region")
    # monotonic start and end postions
    if np.diff(band[0]).min() < 0:
        raise mh.MegaError(
            "Band start positions are not monotonically increasing"
        )
    if np.diff(band[1]).min() < 0:
        raise mh.MegaError(
            "Band end positions are not monotonically increasing"
        )

    # if provided check that start and end coordinates agree with signal and
    # levels.
    if is_sig_band:
        if sig_len is not None and band.shape[1] != sig_len:
            LOGGER.debug(f"Invalid sig_band length: {band.shape[1]} {sig_len}")
            raise mh.MegaError("Invalid sig_band length")
        if seq_len is not None and band[1, -1] != seq_len:
            LOGGER.debug(
                f"Invalid sig_band end coordinate: {band[1, -1]} {seq_len}"
            )
            raise mh.MegaError("Invalid sig_band end coordinate")
    else:
        if sig_len is not None and band[1, -1] != sig_len:
            LOGGER.debug(
                f"Invalid seq_band end coordinate: {band[1, -1]} {sig_len}"
            )
            raise mh.MegaError("Invalid seq_band end coordinate")
        if seq_len is not None and band.shape[1] != seq_len:
            LOGGER.debug(f"Invalid sig_band length: {band.shape[1]} {seq_len}")
            raise mh.MegaError("Invalid sig_band length")
