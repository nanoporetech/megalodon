import os
import sys
from glob import glob
from time import sleep

import numpy as np
from ont_fast5_api.fast5_interface import get_fast5_file as ont_get_fast5_file

from megalodon import logging, megalodon_helper as mh


LOGGER = logging.get_logger()

LIVE_COMP_FN_START = "final_summary_"
LIVE_SLEEP = 0.1


class LiveDoneError(Exception):
    """ Custom error to indicate that live processing has completed
    """
    pass


def get_read_ids(fast5_fn):
    with ont_get_fast5_file(fast5_fn, 'r') as fast5_fp:
        return fast5_fp.get_read_ids()


def iterate_fast5_filenames(input_path, recursive=True, do_it_live=False):
    """ Iterate over fast5 file from the base directory

    Args:
        input_path (str): Path to root reads directory
        recursive (bool): Search recursively
        do_it_live (bool): Raise error when file starting with
            LIVE_COMP_FN_START is found

    Raises:
        megalodon.fast5_io.LiveDoneError: When do_it_live is set and
            LIVE_COMP_FN_START is found
    """
    if recursive:
        for root, _, fns in os.walk(input_path, followlinks=True):
            for fn in fns:
                if do_it_live and fn.startswith(LIVE_COMP_FN_START):
                    raise LiveDoneError
                if not fn.endswith('.fast5'):
                    continue
                yield os.path.join(root, fn)
    else:
        for fn in glob(os.path.join(input_path, '*.fast5')):
            yield fn


def _iterate_fast5_reads_core(
        fast5s_dir, limit=None, recursive=True, do_it_live=False,
        skip_fns=None):
    """ Yield reads in a directory of fast5 files.

    Each read is specified by a tuple (filepath, read_id)
    Files may be single or multi-read fast5s

    Args:
        fast5s_dir (str): Directory containing fast5 files
        limit (int): Limit number of reads to consider
        recursive (bool): Search path recursively for fast5 files.
        do_it_live (bool): Raise RuntimeError when file starting with
            LIVE_COMP_FN_START is found.
        skip_fns (set): Set of filenames to skip for read iteration
    """
    nreads = 0
    for fast5_fn in iterate_fast5_filenames(fast5s_dir, recursive, do_it_live):
        if skip_fns is not None and fast5_fn in skip_fns:
            continue
        for read_id in get_read_ids(fast5_fn):
            yield fast5_fn, read_id
            nreads += 1
            if limit is not None and nreads >= limit:
                return


def iterate_fast5_reads(
        fast5s_dir, limit=None, recursive=True, do_it_live=False):
    """ Iterate fast5 file path and read id combinations in directory

    Args:
        fast5s_dir (str): Directory containing fast5 files
        limit (int): Limit number of reads to consider
        recursive (bool): Search path recursively for fast5 files.
        do_it_live (bool): Continue searching for reads until file starting
            with LIVE_COMP_FN_START is found.
    """
    if do_it_live:
        LOGGER.debug('LiveProcessingStarting')
        # track files that have been searched for read_ids so they aren't
        # opened again
        used_fns = set()
        # keep track of limit here since inner loop won't keep total count
        nreads = 0
        try:
            while True:
                # freeze set so it doesn't update in inner function
                iter_skip_fns = frozenset(used_fns)
                for fast5_fn, read_id in _iterate_fast5_reads_core(
                        fast5s_dir, recursive=recursive, do_it_live=True,
                        skip_fns=iter_skip_fns):
                    yield fast5_fn, read_id
                    used_fns.add(fast5_fn)
                    nreads += 1
                    if limit is not None and nreads >= limit:
                        return
                # sleep so file searching doesn't take too much processing
                sleep(LIVE_SLEEP)
        except LiveDoneError:
            LOGGER.debug('LiveProcessingComplete')
            # search for any files that were missed before run ended
            for fast5_fn, read_id in _iterate_fast5_reads_core(
                    fast5s_dir, recursive=recursive, skip_fns=used_fns):
                yield fast5_fn, read_id
                nreads += 1
                if limit is not None and nreads >= limit:
                    return
    else:
        for fast5_fn, read_id in _iterate_fast5_reads_core(
                fast5s_dir, limit=limit, recursive=recursive):
            yield fast5_fn, read_id


def get_read(fast5_fn, read_id):
    return ont_get_fast5_file(fast5_fn, mode="r").get_read(read_id)


def get_fast5_file(fast5_fn):
    return ont_get_fast5_file(fast5_fn, mode="r")


def get_signal(read, scale=True):
    """ Get raw signal from read.
    """
    try:
        raw_sig = read.get_raw_data()
    except IOError:
        raise mh.MegaError('Error extracting raw data. Ensure VBZ plugin is ' +
                           'installed (if applicable).')

    if scale:
        med, mad = mh.med_mad(raw_sig)
        raw_sig = ((raw_sig - med) / mad).astype(np.float32)

    return raw_sig


def get_posteriors(read):
    # extract guppy StateData and calls
    latest_basecall = read.get_latest_analysis('Basecall_1D')
    state_data = read.get_analysis_dataset(
        latest_basecall + '/BaseCalled_template', 'StateData')
    state_attrs = read.get_analysis_attributes(
        latest_basecall + '/BaseCalled_template/StateData')
    if state_data is None or state_attrs is None:
        raise mh.MegaError(
            'StateData not found in FAST5 file. Ensure --fast5_out and ' +
            '--post_out were set when running guppy.')
    # convert state data from integers to float values
    posteriors = (state_data.astype(np.float32) + state_attrs['offset']) * \
        state_attrs['scale']

    return posteriors


def get_stride(read):
    latest_basecall = read.get_latest_analysis('Basecall_1D')
    return read.get_summary_data(
        latest_basecall)['basecall_1d_template']['block_stride']


def get_mod_base_info(read):
    latest_basecall = read.get_latest_analysis('Basecall_1D')
    mod_attrs = read.get_analysis_attributes(
        latest_basecall + '/BaseCalled_template/ModBaseProbs')
    if mod_attrs is None:
        return '', mh.ALPHABET
    mod_base_long_names = mod_attrs['modified_base_long_names']
    mod_alphabet = mod_attrs['output_alphabet']
    return mod_base_long_names, mod_alphabet


def get_signal_trim_coordiates(read):
    seg_name = read.get_latest_analysis('Segmentation')
    trim_start = read.get_summary_data(seg_name)['segmentation'][
        'first_sample_template']
    trim_len = read.get_summary_data(seg_name)['segmentation'][
        'duration_template']
    return trim_start, trim_len


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
