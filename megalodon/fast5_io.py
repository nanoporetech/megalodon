import os
from glob import glob

from ont_fast5_api.fast5_interface import get_fast5_file

from megalodon import megalodon_helper as mh


def iterate_fast5_filenames(input_path, recursive):
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

def get_signal(fast5_fn, read_id, scale=True):
    """ Get raw signal from read.
    """
    with get_fast5_file(fast5_fn, 'r') as fast5_fp:
        raw_sig = fast5_fp.get_read(read_id).get_raw_data()

    if scale:
        med, mad = mh.med_mad(raw_sig)
        raw_sig = (raw_sig - med) / mad

    return raw_sig


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
