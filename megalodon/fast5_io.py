import os
import sys
import queue
import threading
import traceback
from glob import glob
from time import sleep
import multiprocessing as mp

import numpy as np
from ont_fast5_api.fast5_interface import get_fast5_file as ont_get_fast5_file

from megalodon import logging, megalodon_helper as mh


LOGGER = logging.get_logger()

LIVE_COMP_FN_START = "final_summary_"
LIVE_SLEEP = 0.1

_PROFILE_EXTRACT_SIGNAL = False


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


####################################
# Threaded Read/Signal Enumeration #
####################################

def _extract_signal_worker(
        fn_read_ids_q, signal_q, model_info, extract_dacs, aux_failed_q):
    LOGGER.debug('Starting')
    try:
        # while there are reads to process continue extracting signal
        while True:
            try:
                fn_rids = fn_read_ids_q.get(timeout=0.01)
            except queue.Empty:
                continue
            # None indicates end of enumerated reads
            if fn_rids is None:
                LOGGER.debug('Closing')
                break

            fast5_fn, read_ids = fn_rids
            with get_fast5_file(fast5_fn) as fast5_fp:
                for read_id in read_ids:
                    sig_info, seq_summ_info = model_info.extract_signal_info(
                        fast5_fp, read_id, extract_dacs)
                    signal_q.put((tuple(sig_info), tuple(seq_summ_info)))
    except Exception as e:
        aux_failed_q.put((
            'ExtractSignalProcessingError', str(e), traceback.format_exc()))


if _PROFILE_EXTRACT_SIGNAL:
    _extract_signal_wrapper = _extract_signal_worker

    def _extract_signal_worker(*args):
        import cProfile
        cProfile.runctx('_extract_signal_wrapper(*args)', globals(), locals(),
                        filename='megalodon_extract_signal.prof')


def _extract_signal(
        fn_read_ids_q, signal_q, aux_failed_q, input_info, model_info,
        extract_dacs):
    try:
        LOGGER.debug('Starting')
        sig_extract_ts = list()
        for se_i in range(input_info.num_read_enum_ts):
            sig_extract_ts.append(threading.Thread(
                target=_extract_signal_worker,
                args=(fn_read_ids_q, signal_q, model_info, extract_dacs,
                      aux_failed_q),
                daemon=True, name='ExtractSigThread{:03d}'.format(se_i)))
            sig_extract_ts[-1].start()
        LOGGER.debug('InitComplete')
    except KeyboardInterrupt:
        raise
    except Exception as e:
        aux_failed_q.put((
            'ExtractSignalInitError', str(e), traceback.format_exc()))
        return

    try:
        # while there are reads to process continue extracting signal
        while any(t.is_alive() for t in sig_extract_ts):
            sleep(0.01)
        LOGGER.debug('Closing')
    except KeyboardInterrupt:
        raise
    except Exception as e:
        aux_failed_q.put((
            'ExtractSignalProcessingError', str(e), traceback.format_exc()))


def _fill_files_file_enum_worker(
        fast5_fn_q, file_reads_q, is_below_reads_limit):
    LOGGER.debug('Starting')
    while is_below_reads_limit():
        try:
            fast5_fn = fast5_fn_q.get(timeout=0.01)
        except queue.Empty:
            continue
        # None filename sent signifies end of file enumeration
        if fast5_fn is None:
            break
        file_read_ids = get_read_ids(fast5_fn)
        file_reads_q.put((fast5_fn, file_read_ids))
        LOGGER.debug('ReadIDsExtractedFrom: {} {}'.format(
            fast5_fn, len(file_read_ids)))
    LOGGER.debug('Closing')


def _fill_files_queue_worker(fast5_fn_q, input_info, is_below_reads_limit):
    LOGGER.debug('Starting')
    input_info = mh.INPUT_INFO(*input_info)
    used_fns = set()
    num_fns = 0
    try:
        for fast5_fn in iterate_fast5_filenames(
                input_info.fast5s_dir, recursive=input_info.recursive,
                do_it_live=input_info.do_it_live):
            if not is_below_reads_limit():
                break
            fast5_fn_q.put(fast5_fn)
            num_fns += 1
            if input_info.do_it_live:
                used_fns.add(fast5_fn)
    except LiveDoneError:
        LOGGER.debug('LiveProcessingComplete')
        for fast5_fn in iterate_fast5_filenames(
                input_info.fast5s_dir, recursive=input_info.recursive):
            if fast5_fn not in used_fns:
                fast5_fn_q.put(fast5_fn)

    # add None to indicate to read enumeration workers that file enumeration
    # is complete
    for _ in range(input_info.num_read_enum_ts):
        fast5_fn_q.put(None)
    LOGGER.debug('Found {} total FAST5 files'.format(num_fns))
    LOGGER.debug('Closing')


def _fill_files_queue(input_info, fn_read_ids_q, num_reads_conn, aux_failed_q):
    """ Fill fn_read_ids_q with (fast5_fn, read_ids) tuples. When finished, send
    total number of reads to num_reads_conn.

    In order to maximize performance, read enumeration is spread over a
    multi-threading interface.
    """
    try:
        LOGGER.debug('Starting')
        valid_read_ids = mh.parse_read_ids(input_info.read_ids_fn)
        used_read_ids = set()
        fast5_fn_q = mp.Queue()
        file_reads_q = mp.Queue()
        below_reads_limit = True
        # spawn thread to enumerate files
        file_enum_t = threading.Thread(
            target=_fill_files_queue_worker,
            args=(fast5_fn_q, tuple(input_info), lambda: below_reads_limit),
            daemon=True, name='FileEnum')
        file_enum_t.start()
        # spawn threads to enumerate reads within files
        enum_ts = [file_enum_t]
        for re_i in range(input_info.num_read_enum_ts):
            enum_ts.append(threading.Thread(
                target=_fill_files_file_enum_worker,
                args=(fast5_fn_q, file_reads_q, lambda: below_reads_limit),
                daemon=True, name='ReadEnumThread{:03d}'.format(re_i)))
            enum_ts[-1].start()
        LOGGER.debug('InitComplete')
    except KeyboardInterrupt:
        raise
    except Exception as e:
        aux_failed_q.put((
            'FillFilesInitError', str(e), traceback.format_exc()))
        return

    try:
        # fill queue with read filename and read id tuples
        while below_reads_limit and (
                any(t.is_alive() for t in enum_ts) or
                not file_reads_q.empty()):
            try:
                fast5_fn, file_read_ids = file_reads_q.get(timeout=0.01)
            except queue.Empty:
                continue
            valid_file_read_ids = set(file_read_ids).difference(used_read_ids)
            if len(valid_file_read_ids) < len(file_read_ids):
                LOGGER.debug('RepeatedReadIDs {} {}'.format(
                    fast5_fn,
                    len(set(file_read_ids).intersection(used_read_ids))))
            if valid_read_ids is not None:
                valid_file_read_ids.intersection_update(valid_read_ids)
            if len(valid_file_read_ids) == 0:
                continue
            if input_info.num_reads is not None and \
               len(used_read_ids) + len(valid_file_read_ids) >= \
               input_info.num_reads:
                LOGGER.debug('NumReadsLimitReached')
                valid_file_read_ids = list(valid_file_read_ids)[
                    :input_info.num_reads - len(used_read_ids)]
                fn_read_ids_q.put((fast5_fn, valid_file_read_ids))
                used_read_ids.update(valid_file_read_ids)
                below_reads_limit = False
                break
            fn_read_ids_q.put((fast5_fn, valid_file_read_ids))
            used_read_ids.update(valid_file_read_ids)
    except KeyboardInterrupt:
        raise
    except Exception as e:
        aux_failed_q.put((
            'FillFilesProcessingError', str(e), traceback.format_exc()))
    finally:
        num_reads_conn.send(len(used_read_ids))
        # add None to indicate to each signal enumeration thread that
        # read enumeration is complete
        for _ in range(input_info.num_read_enum_ts *
                       input_info.num_extract_sig_proc):
            fn_read_ids_q.put(None)


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
