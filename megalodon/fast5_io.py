import os
from glob import glob

import numpy as np

from ont_fast5_api.fast5_interface import get_fast5_file

from megalodon import megalodon_helper as mh


def iterate_fast5_filenames(input_path, recursive=True):
    if recursive:
        for root, _, fns in os.walk(input_path, followlinks=True):
            for fn in fns:
                if not fn.endswith('.fast5'): continue
                yield os.path.join(root, fn)
    else:
        for fn in glob(os.path.join(input_path, '*.fast5')):
            yield fn
    return

def iterate_fast5_reads(fast5s_dir, limit=None, recursive=True):
    """Return iterator yielding reads in a directory of fast5 files or a
    single fast5 file.

    Each read is specified by a tuple (filepath, read_id)
    Files may be single or multi-read fast5s

    :param fast5s_dir: Directory containing fast5 files
    :param limit: Limit number of reads to consider
    :param recursive: Search path recursively for fast5 files.
    """
    nreads = 0
    for fast5_fn in iterate_fast5_filenames(fast5s_dir, recursive):
        with get_fast5_file(fast5_fn, 'r') as fast5_fp:
            for read_id in fast5_fp.get_read_ids():
                yield fast5_fn, read_id
                nreads += 1
                if limit is not None and nreads >= limit:
                    return
    return

def get_read(fast5_fn, read_id):
    return get_fast5_file(fast5_fn, mode="r").get_read(read_id)

def get_signal(read, scale=True):
    """ Get raw signal from read.
    """
    raw_sig = read.get_raw_data()

    if scale:
        med, mad = mh.med_mad(raw_sig)
        raw_sig = (raw_sig - med) / mad

    return raw_sig

def get_posteriors(read):
    # extract guppy StateData and calls
    latest_basecall = read.get_latest_analysis('Basecall_1D')
    state_data = read.get_analysis_dataset(
        latest_basecall + '/BaseCalled_template', 'StateData', proxy=False)
    state_attrs = read.get_analysis_attributes(
        latest_basecall + '/BaseCalled_template/StateData')
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
    latest_seg = read.get_latest_analysis('Segmentation')
    trim_start = read.get_summary_data(seg_name)['segmentation'][
        'first_sample_template']
    trim_len = read.get_summary_data(seg_name)['segmentation'][
        'duration_template']
    return trim_start, trim_len


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
