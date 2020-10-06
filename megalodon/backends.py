import os
import re
import sys
import array
import subprocess
from abc import ABC
from time import sleep
from distutils.version import LooseVersion
from collections import defaultdict, namedtuple

import numpy as np

from megalodon import decode, fast5_io, logging, megalodon_helper as mh


LOGGER = logging.get_logger()

# model type specific information
TAI_NAME = 'taiyaki'
FAST5_NAME = 'fast5'
PYGUPPY_NAME = 'pyguppy'

# parameters for each backend run mode
TAI_PARAMS = namedtuple('TAI_PARAMS', (
    'available', 'taiyaki_model_fn', 'devices', 'chunk_size',
    'chunk_overlap', 'max_concur_chunks'))
TAI_PARAMS.__new__.__defaults__ = tuple([None, ] * 5)
FAST5_PARAMS = namedtuple('FAST5_PARAMS', (
    'available', 'fast5s_dir', 'num_startup_reads'))
FAST5_PARAMS.__new__.__defaults__ = tuple([None, ] * 2)
PYGUPPY_PARAMS = namedtuple('PYGUPPY_PARAMS', (
    'available', 'config', 'bin_path', 'port', 'timeout', 'devices', 'out_dir',
    'server_params'))
PYGUPPY_PARAMS.__new__.__defaults__ = tuple([None, ] * 7)
BACKEND_PARAMS = namedtuple('BACKEND_PARAMS', (
    TAI_NAME, FAST5_NAME, PYGUPPY_NAME))

COMPAT_GUPPY_MODEL_TYPES = set(('flipflop',))
GUPPY_HOST = 'localhost'
PYGUPPY_PER_TRY_TIMEOUT = 0.01
PYGUPPY_MAX_RECONNECT_ATTEMPTS = 5
GUPPY_LOG_BASE = 'guppy_log'
GUPPY_PORT_PAT = re.compile(r'Starting server on port:\W+(\d+)')
GUPPY_VERSION_PAT = re.compile(
    r'Oxford Nanopore Technologies, Limited. ' +
    r'Version\W+([0-9]+\.[0-9]+\.[0-9]+)\+[0-9a-z]+')
MIN_GUPPY_VERSION = LooseVersion('4.0')

PYGUPPY_CLIENT_KWARGS = {
    'move_and_trace_enabled': True,
    'state_data_enabled': True
}
PYGUPPY_CALLED_READ = namedtuple('PYGUPPY_CALLED_READ', (
    'model_type', 'model_stride', 'mod_long_names', 'output_alphabet',
    'state_size', 'trimmed_samples', 'scaling_shift', 'scaling_scale',
    'move', 'state', 'seq', 'qual'))

# maximum time (in seconds) to wait before assigning device
# over different processes. Actual time choosen randomly
# as many simultaneous calls sometimes conflict and kill
# some device assignments
MAX_DEVICE_WAIT = 1.0

SIGNAL_DATA = namedtuple('SIGNAL_DATA', (
    'fast5_fn', 'read_id', 'raw_len', 'dacs', 'raw_signal',
    'scale_params', 'stride', 'posteriors', 'channel_info'))
# set default value of None for ref, alts and ref_start
SIGNAL_DATA.__new__.__defaults__ = tuple([None, ] * 6)


def parse_device(device):
    if device is None or device == 'cpu':
        return 'cpu'
    # add colon for UGE/SGE device type
    if re.match('cuda[0-9]+', device) is not None:
        return 'cuda:{}'.format(device[4:])
    # if integer device is specified
    elif not device.startswith('cuda'):
        return 'cuda:{}'.format(device)
    return device


def parse_backend_params(args, num_fast5_startup_reads=5):
    # parse taiyaki backend params
    if any(not hasattr(args, k) for k in (
            'chunk_size', 'chunk_overlap', 'max_concurrent_chunks',
            'taiyaki_model_filename', 'devices')) or \
            args.taiyaki_model_filename is None:
        tai_params = TAI_PARAMS(False)
    else:
        tai_model_fn = mh.resolve_path(args.taiyaki_model_filename)
        if all(param is not None for param in (
                tai_model_fn, args.chunk_size, args.chunk_overlap,
                args.max_concurrent_chunks)):
            tai_params = TAI_PARAMS(
                True, tai_model_fn, args.devices, args.chunk_size,
                args.chunk_overlap, args.max_concurrent_chunks)
        else:
            tai_params = TAI_PARAMS(False)

    # parse fast5 post_out backend params
    if hasattr(args, 'fast5s_dir') and args.fast5s_dir is not None:
        fast5_params = FAST5_PARAMS(
            True, args.fast5s_dir, num_fast5_startup_reads)
    else:
        fast5_params = FAST5_PARAMS(False)

    # parse pyguppy backend params
    if any(not hasattr(args, k) for k in (
            'guppy_config', 'guppy_server_path', 'guppy_server_port',
            'guppy_timeout', 'devices', 'output_directory', 'guppy_params',
            'do_not_use_guppy_server')) or args.do_not_use_guppy_server:
        pyguppy_params = PYGUPPY_PARAMS(False)
    else:
        if args.guppy_server_port is None:
            args.guppy_server_port = 'auto'
        pyguppy_params = PYGUPPY_PARAMS(
            available=True, config=args.guppy_config,
            bin_path=args.guppy_server_path, port=args.guppy_server_port,
            timeout=args.guppy_timeout, devices=args.devices,
            out_dir=args.output_directory, server_params=args.guppy_params)

    if hasattr(args, 'basecalls_format') and \
       args.basecalls_format == 'fastq' and \
       not pyguppy_params.available:
        LOGGER.warning(
            'Quality score computation not enabled for taiyaki or FAST5 ' +
            'backends. Quality scores will be invalid.')

    return BACKEND_PARAMS(tai_params, fast5_params, pyguppy_params)


def extract_seq_summary_info(read, channel_info, na_str='NA'):
    """ Extract non-basecalling sequencing summary information from
    ont_fast5_api read object
    """
    try:
        fn = read.filename
        read_id = read.read_id
        samp_rate = channel_info[mh.CHAN_INFO_SAMP_RATE]
        try:
            raw_attrs = read.handle[read.raw_dataset_group_name].attrs
            mux = raw_attrs['start_mux']
            start_time = '{:.6f}'.format(raw_attrs['start_time'] / samp_rate)
            dur = '{:.6f}'.format(raw_attrs['duration'] / samp_rate)
        except AttributeError:
            mux = start_time = dur = na_str
        run_id = read.run_id
        try:
            run_id = run_id.decode()
        except AttributeError:
            pass
        batch_id = na_str
        chan = channel_info[mh.CHAN_INFO_CHANNEL_SLOT]
    except Exception:
        # if anything goes wrong set all values to na_str
        fn = read_id = run_id = batch_id = chan = mux = start_time = \
            dur = na_str
    return mh.SEQ_SUMM_INFO(
        filename=fn, read_id=read_id, run_id=run_id, batch_id=batch_id,
        channel=chan, mux=mux, start_time=start_time, duration=dur)


def _log_softmax_axis1(x):
    """ Compute log softmax over axis=1
    """
    e_x = np.exp((x.T - np.max(x, axis=1)).T)
    with np.errstate(divide='ignore'):
        return np.log((e_x.T / e_x.sum(axis=1)).T)


def get_pyguppy_read(read_id, raw_data, channel_info):
    return {
        'read_tag': np.random.randint(0, int(2**32 - 1)),
        'read_id': read_id,
        'raw_data': raw_data,
        'daq_offset': float(channel_info[mh.CHAN_INFO_OFFSET]),
        'daq_scaling': float(channel_info[mh.CHAN_INFO_RANGE]) / channel_info[
            mh.CHAN_INFO_DIGI]}


def parse_pyguppy_called_read(called_read):
    read_metadata = called_read['metadata']
    out_alphabet = read_metadata['base_mod_alphabet'] \
        if len(read_metadata['base_mod_alphabet']) > 0 \
        else mh.ALPHABET
    read_datasets = called_read['datasets']
    return PYGUPPY_CALLED_READ(
        model_type=read_metadata['basecall_type'],
        model_stride=read_metadata['model_stride'],
        mod_long_names=read_metadata['base_mod_long_names'].split(),
        output_alphabet=out_alphabet, state_size=read_metadata['state_size'],
        trimmed_samples=read_metadata['trimmed_samples'],
        scaling_shift=read_metadata['scaling_median'],
        scaling_scale=read_metadata['scaling_med_abs_dev'],
        move=read_datasets['movement'], state=read_datasets['state_data'],
        seq=read_datasets['sequence'], qual=read_datasets['qstring'])


class AbstractModelInfo(ABC):
    @property
    def n_can_state(self):
        ncan_base = len(self.can_alphabet)
        return (ncan_base + ncan_base) * (ncan_base + 1)

    def _compute_mod_alphabet_attrs(self):
        """ Parse alphabet attributes into more user-friendly data structures.

        Requires the following attributes to be set:
            - `can_nmods` (list): Number of modified bases associated with each
                canonical base.
            - `output_alphabet` (str)
            - `ordered_mod_long_names` (list)
        """
        self.can_alphabet = ''
        self.can_indices = []
        self.mod_long_names = []
        self.mod_base_to_can = {}
        self.str_to_int_mod_labels = {}
        self.can_base_mods = defaultdict(list)
        curr_can_offset = 0
        curr_nmods = 0
        for can_base_nmods in self.can_nmods:
            can_base = self.output_alphabet[curr_can_offset]
            self.can_alphabet += can_base
            self.can_indices.append(curr_can_offset)
            for mod_i, mod_base in enumerate(self.output_alphabet[
                    curr_can_offset + 1:
                    curr_can_offset + can_base_nmods + 1]):
                self.mod_long_names.append((
                    mod_base, self.ordered_mod_long_names[
                        curr_nmods + mod_i]))
                self.str_to_int_mod_labels[mod_base] = mod_i + 1
                self.can_base_mods[can_base].append(mod_base)
                self.mod_base_to_can[mod_base] = can_base

            curr_can_offset += can_base_nmods + 1
            curr_nmods += can_base_nmods

        self.can_indices.append(curr_can_offset)
        self.can_indices = np.array(self.can_indices).astype(np.uintp)
        self.can_base_mods = dict(self.can_base_mods)

    def _parse_minimal_alphabet_info(self):
        """ Parse minimal alphabet information pertaining to a model.

        The following attributes should be set:
            - `ordered_mod_long_names` (list)
            - `output_alphabet` (str)
        """
        # compute values required for standard model attribute extraction
        self.n_mods = len(self.ordered_mod_long_names)
        self.is_cat_mod = self.n_mods > 0
        self.can_alphabet = None
        for v_alphabet in mh.VALID_ALPHABETS:
            if all(b in self.output_alphabet for b in v_alphabet):
                self.can_alphabet = v_alphabet
                break
        if self.can_alphabet is None:
            LOGGER.error(
                'Model information contains invalid alphabet ({})'.format(
                    self.output_alphabet))
            raise mh.MegaError('Invalid alphabet.')
        # compute number of modified bases for each canonical base

        if self.is_cat_mod:
            self.can_nmods = np.diff(np.array(
                [self.output_alphabet.index(b) for b in self.can_alphabet] +
                [len(self.output_alphabet), ])) - 1
            self._compute_mod_alphabet_attrs()
            # compute indices over which to compute softmax for mod bases
            self.can_raw_mod_indices = []
            curr_n_mods = 0
            for bi_nmods in self.can_nmods:
                if bi_nmods == 0:
                    self.can_raw_mod_indices.append(None)
                else:
                    # global canonical category is index 0 then mod cat indices
                    self.can_raw_mod_indices.append(np.insert(
                        np.arange(curr_n_mods + 1, curr_n_mods + 1 + bi_nmods),
                        0, 0))
                curr_n_mods += bi_nmods
        else:
            self.str_to_int_mod_labels = {}
            self.mod_long_names = []
            self.can_nmods = None

    def format_mod_scores(self, bc_seq, mods_scores, min_prob):
        """ Convert 2D array of scores to hts-spec modified base output.
        See https://github.com/samtools/hts-specs/pull/418

        Args:
            bc_seq (str): basecall sequence
            mods_scores (np.ndarray): 2D array with basecall position rows and
                modbase columns.
            min_prob (float): Minimum probability to include modified base

        Returns:
            Mm string tag and Ml array tag
        """
        mm_tag, ml_tag = '', array.array('B')
        prev_bases = 0
        for can_base, can_nmods in zip(self.can_alphabet, self.can_nmods):
            if can_nmods == 0:
                prev_bases += 1
                continue
            mod_bases = self.can_base_mods[can_base]
            if len(mod_bases) != can_nmods:
                raise mh.MegaError((
                    'Number of modified bases ({}) associated with {} does ' +
                    'not match expected number of columns in mod scores: ' +
                    '{}.').format(','.join(mod_bases), can_base, can_nmods))
            can_bc_pos = np.array([b == can_base for b in bc_seq], dtype=bool)
            for mod_base, mod_index in zip(mod_bases, range(
                    prev_bases + 1, prev_bases + 1 + can_nmods)):
                probs = np.exp(mods_scores[can_bc_pos, mod_index])
                valid_prob_locs = np.where(probs > min_prob)[0]
                mm_tag += '{}+{}{};'.format(
                    can_base, mh.convert_legacy_mods(mod_base),
                    ''.join(',{}'.format(d) for d in np.diff(np.insert(
                        valid_prob_locs, 0, -1)) - 1))
                # extract mod scores and scale to 0-255 range
                scaled_probs = np.floor(probs[valid_prob_locs] * 256)
                # last interval includes prob=1
                scaled_probs[scaled_probs == 256] = 255
                ml_tag.extend(scaled_probs.astype(np.uint8))
            prev_bases += can_nmods + 1

        return mm_tag, ml_tag

    def get_alphabet_str(self):
        if self.is_cat_mod:
            return 'Using canonical alphabet {} and modified bases {}'.format(
                self.can_alphabet, '; '.join(
                    '{}={} (alt to {})'.format(
                        mod_b, mln, self.mod_base_to_can[mod_b])
                    for mod_b, mln in self.mod_long_names))
        return 'Using canonical alphabet {}'.format(self.can_alphabet)


class DetachedModelInfo(AbstractModelInfo):
    """ DetachedModelInfo represents a wrapper similar to ModelInfo, but allows
    manual setting of attributes instead of loading from a real model.
    """

    def __init__(self, alphabet=mh.ALPHABET, mod_long_names=None):
        self.output_alphabet = alphabet
        self.ordered_mod_long_names = mod_long_names
        if self.ordered_mod_long_names is None:
            self.ordered_mod_long_names = []
        self._parse_minimal_alphabet_info()
        self.output_size = (41 + len(self.mod_long_names)
                            if self.is_cat_mod else 40)


class ModelInfo(AbstractModelInfo):
    """ ModelInfo wraps the model backends supported by Megalodon. Currently,
    this includes (in preference order) guppy, taiyaki and fast5. Class
    initialization will load the first backend available as found in
    backend_params.

    Useful methods include:
        - `prep_model_worker`: Load model onto GPU device
        - `extract_signal_info`: Extract signal information
        - `iter_basecalled_reads`: Basecall a batch of reads yielding
            signal info, sequencing summary info, sequence, quality values,
            basecall positions within posterior array, canonical only posterior
            array, posterior array with modified base scores, and basecall
            anchored modified base scores.
    """

    def _load_taiyaki_model(self):
        LOGGER.info('Loading taiyaki basecalling backend')
        self.model_type = TAI_NAME

        devices = self.params.taiyaki.devices
        if devices is None:
            devices = ['cpu', ]
        self.process_devices = [
            parse_device(devices[i]) for i in np.tile(
                np.arange(len(devices)),
                (self.num_proc // len(devices)) + 1)][:self.num_proc]

        try:
            # import modules
            from taiyaki.helpers import (
                load_model as load_taiyaki_model, guess_model_stride)
            from taiyaki.basecall_helpers import run_model as tai_run_model
            from taiyaki.layers import GlobalNormFlipFlopCatMod
        except ImportError:
            LOGGER.error(
                'Failed to import taiyaki. Ensure working ' +
                'installations to run megalodon')
            sys.exit(1)
        try:
            import torch
        except ImportError:
            LOGGER.error(
                'Failed to import pytorch. Ensure working ' +
                'installations to run megalodon')
            sys.exit(1)

        # store modules in object
        self.load_taiyaki_model = load_taiyaki_model
        self.tai_run_model = tai_run_model
        self.torch = torch

        tmp_model = self.load_taiyaki_model(
            self.params.taiyaki.taiyaki_model_fn)
        ff_layer = tmp_model.sublayers[-1]
        self.is_cat_mod = (
            GlobalNormFlipFlopCatMod is not None and isinstance(
                ff_layer, GlobalNormFlipFlopCatMod))
        self.stride = guess_model_stride(tmp_model)
        self.output_size = ff_layer.size
        if self.is_cat_mod:
            # Modified base model is defined by 3 fixed fields in taiyaki
            # can_nmods, output_alphabet and modified_base_long_names
            self.output_alphabet = ff_layer.output_alphabet
            self.can_nmods = ff_layer.can_nmods
            self.ordered_mod_long_names = ff_layer.ordered_mod_long_names
            self._compute_mod_alphabet_attrs()
        else:
            if mh.nstate_to_nbase(self.output_size) != 4:
                raise NotImplementedError(
                    'Naive modified base flip-flop models are not ' +
                    'supported.')
            self.output_alphabet = mh.ALPHABET
            self.can_alphabet = mh.ALPHABET
            self.mod_long_names = []
            self.str_to_int_mod_labels = {}
            self.can_nmods = None
        self.n_mods = len(self.mod_long_names)

    def _load_fast5_post_out(self):
        def get_model_info_from_fast5(read):
            try:
                stride = fast5_io.get_stride(read)
                mod_long_names, out_alphabet = fast5_io.get_mod_base_info(read)
                out_size = fast5_io.get_posteriors(read).shape[1]
                mod_long_names = mod_long_names.split()
            except KeyError:
                LOGGER.error(
                    'Fast5 read does not contain required attributes.')
                raise mh.MegaError(
                    'Fast5 read does not contain required attributes.')
            return stride, mod_long_names, out_alphabet, out_size

        LOGGER.info('Loading FAST5 basecalling backend')
        self.model_type = FAST5_NAME
        self.process_devices = [None, ] * self.num_proc

        read_iter = fast5_io.iterate_fast5_reads(self.params.fast5.fast5s_dir)
        nreads = 0
        try:
            fast5_fn, read_id = next(read_iter)
            read = fast5_io.get_read(fast5_fn, read_id)
            (self.stride, self.ordered_mod_long_names, self.output_alphabet,
             self.output_size) = get_model_info_from_fast5(read)
        except StopIteration:
            LOGGER.error('No reads found.')
            raise mh.MegaError('No reads found.')

        for fast5_fn, read_id in read_iter:
            read = fast5_io.get_read(fast5_fn, read_id)
            r_s, r_omln, r_oa, r_os = get_model_info_from_fast5(read)
            if (
                    self.stride != r_s or
                    self.ordered_mod_long_names != r_omln or
                    self.output_alphabet != r_oa or
                    self.output_size != r_os):
                LOGGER.error(
                    'Model information from FAST5 files is inconsistent. ' +
                    'Assure all reads were called with the same model.')
                raise mh.MegaError(
                    'Model information from FAST5 files is inconsistent.')
            nreads += 1
            if nreads >= self.params.fast5.num_startup_reads:
                break

        self._parse_minimal_alphabet_info()

    def _load_pyguppy(self, init_sig_len=1000):
        def _check_guppy_version(pyguppy_version_str):
            try:
                check_vers_proc = subprocess.Popen(
                    [self.params.pyguppy.bin_path, '--version'], shell=False,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            except FileNotFoundError:
                raise mh.MegaError(
                    'Guppy server executable not found. Please specify path ' +
                    'via `--guppy-server-path` argument.')
            version_out = str(check_vers_proc.stdout.read())
            version_match = GUPPY_VERSION_PAT.search(version_out)
            if version_match is None:
                raise mh.MegaError(
                    'Guppy version string does not match expected ' +
                    'pattern: "{}"'.format(version_out))
            guppy_version_str = version_match.groups()[0]
            LOGGER.debug('Guppy version: "{}"'.format(guppy_version_str))
            LOGGER.debug('Pyguppy version: "{}"'.format(pyguppy_version_str))
            guppy_version = LooseVersion(guppy_version_str)
            if guppy_version < MIN_GUPPY_VERSION:
                raise mh.MegaError(
                    'Megalodon requires Guppy version>=4.0. Got: "{}"'.format(
                        guppy_version_str))
            pyguppy_version = LooseVersion(pyguppy_version_str)
            if len(pyguppy_version.version) < 2 or \
               len(guppy_version.version) < 2:
                LOGGER.warning('Invalid guppy or pyguppy versions.')
            elif guppy_version.version[0] != pyguppy_version.version[0]:
                LOGGER.warning(
                    'Guppy and pyguppy major versions do not match. This ' +
                    'will likely lead to a failure. Please install matching ' +
                    'Guppy and pyguppy versions.')
            elif guppy_version.version[1] != pyguppy_version.version[1]:
                LOGGER.warning(
                    'Guppy and pyguppy minor versions do not match. This ' +
                    'may lead to a failure. Please install matching ' +
                    'Guppy and pyguppy versions.')
            elif guppy_version.version[2] != pyguppy_version.version[2]:
                LOGGER.debug(
                    'Guppy and pyguppy point versions do not match. ' +
                    'Please install matching Guppy and pyguppy versions.')

        def start_guppy_server():
            def get_server_port():
                next_line = guppy_out_read_fp.readline()
                if next_line is None:
                    return None
                try:
                    return int(GUPPY_PORT_PAT.search(next_line).groups()[0])
                except AttributeError:
                    return None

            # set guppy logs output locations
            self.guppy_log = os.path.join(
                self.params.pyguppy.out_dir, GUPPY_LOG_BASE)
            self.guppy_out_fp = open(self.guppy_log + '.out', 'w')
            guppy_out_read_fp = open(self.guppy_log + '.out', 'r')
            self.guppy_err_fp = open(self.guppy_log + '.err', 'w')
            # prepare args to start guppy server
            server_args = [
                self.params.pyguppy.bin_path,
                '-p', str(self.params.pyguppy.port),
                '-l', self.guppy_log,
                '-c', self.params.pyguppy.config,
                '--post_out']
            if self.params.pyguppy.devices is not None and \
               len(self.params.pyguppy.devices) > 0 and \
               self.params.pyguppy.devices[0] != 'cpu':
                devices_str = ' '.join(
                    parse_device(device) for device in
                    self.params.pyguppy.devices)
                server_args.extend(('-x', devices_str))
            if self.params.pyguppy.server_params is not None:
                server_args.extend(self.params.pyguppy.server_params.split())
            LOGGER.debug('guppy server init command: "{}"'.format(
                ' '.join(server_args)))
            self.guppy_server_proc = subprocess.Popen(
                server_args, shell=False,
                stdout=self.guppy_out_fp, stderr=self.guppy_err_fp)
            # wait until server is successfully started or fails
            while True:
                used_port = get_server_port()
                if used_port is not None:
                    break
                if self.guppy_server_proc.poll() is not None:
                    raise mh.MegaError(
                        'Guppy server initialization failed. See guppy logs ' +
                        'in --output-directory for more details.')
                sleep(0.01)
            guppy_out_read_fp.close()
            self.params = self.params._replace(
                pyguppy=self.params.pyguppy._replace(port=used_port))

        def set_pyguppy_model_attributes():
            init_client = self.pyguppy_GuppyBasecallerClient(
                '{}:{}'.format(GUPPY_HOST, self.params.pyguppy.port),
                self.params.pyguppy.config, **PYGUPPY_CLIENT_KWARGS)
            try:
                init_client.connect()
                init_read = [(SIGNAL_DATA(
                    fast5_fn='init_test_read', read_id='init_test_read',
                    raw_len=0, dacs=np.zeros(init_sig_len, dtype=np.int16),
                    channel_info={mh.CHAN_INFO_OFFSET: 0,
                                  mh.CHAN_INFO_RANGE: 1,
                                  mh.CHAN_INFO_DIGI: 1}), None)]
                try:
                    init_called_read, _, _ = next(self.pyguppy_basecall(
                        init_client, init_read))
                except mh.MegaError:
                    raise mh.MegaError(
                        'Failed to run test read with Guppy. See Guppy logs ' +
                        'in --output-directory.')
            except ConnectionError:
                self.close()
                raise mh.MegaError(
                    'Error connecting to Guppy server. Server unavailable.')
            except ValueError:
                self.close()
                raise mh.MegaError(
                    'Error connecting to Guppy server. Missing barcode kit.')
            except RuntimeError as e:
                self.close()
                raise mh.MegaError(
                    'Error connecting to Guppy server. Undefined error: ' +
                    str(e))
            finally:
                init_client.disconnect()

            if init_called_read.model_type not in COMPAT_GUPPY_MODEL_TYPES:
                raise mh.MegaError((
                    'Megalodon is not compatible with guppy model type: ' +
                    '{}').format(init_called_read.model_type))

            self.stride = init_called_read.model_stride
            self.ordered_mod_long_names = init_called_read.mod_long_names
            self.output_alphabet = init_called_read.output_alphabet
            self.output_size = init_called_read.state_size
            if self.ordered_mod_long_names is None:
                self.ordered_mod_long_names = []

        LOGGER.info('Loading guppy basecalling backend')
        self.model_type = PYGUPPY_NAME
        self.process_devices = [None, ] * self.num_proc

        # load necessary packages and store in object attrs
        self.pyguppy_retries = max(
            1, int(self.params.pyguppy.timeout / PYGUPPY_PER_TRY_TIMEOUT))
        from pyguppy_client_lib.pyclient import PyGuppyClient
        from pyguppy_client_lib import __version__ as pyguppy_version_str
        self.pyguppy_GuppyBasecallerClient = PyGuppyClient

        _check_guppy_version(pyguppy_version_str)
        start_guppy_server()
        set_pyguppy_model_attributes()
        self._parse_minimal_alphabet_info()

    def __init__(self, backend_params, num_proc=1):
        self.num_proc = num_proc

        self.params = backend_params

        if self.params.pyguppy.available:
            self._load_pyguppy()
        elif self.params.taiyaki.available:
            self._load_taiyaki_model()
        elif self.params.fast5.available:
            self._load_fast5_post_out()
        else:
            raise mh.MegaError('No basecall model backend enabled.')

    def prep_model_worker(self, device=None):
        """ Load model onto device (when object is loaded into process to run
        basecaller).
        """
        if self.model_type == TAI_NAME:
            # setup for taiyaki model
            self.model = self.load_taiyaki_model(
                self.params.taiyaki.taiyaki_model_fn)
            if device is None or device == 'cpu':
                self.device = self.torch.device('cpu')
            else:
                sleep(np.random.uniform(0, MAX_DEVICE_WAIT))
                try:
                    self.device = self.torch.device(device)
                    self.torch.cuda.set_device(self.device)
                    self.model = self.model.to(self.device)
                except RuntimeError:
                    LOGGER.error('Invalid CUDA device: {}'.format(device))
                    raise mh.MegaError('Error setting CUDA GPU device.')
            self.model = self.model.eval()
        elif self.model_type == PYGUPPY_NAME:
            # open guppy client interface (None indicates using config
            # from server)
            self.client = self.pyguppy_GuppyBasecallerClient(
                '{}:{}'.format(GUPPY_HOST, self.params.pyguppy.port),
                self.params.pyguppy.config, **PYGUPPY_CLIENT_KWARGS)
            self.client.connect()

    def extract_signal_info(self, fast5_fp, read_id, extract_dacs=False):
        """ Extract signal information from fast5 file pointer.

        Args:
            fast5_fp (:ont_fast5_api.fast5_file:`Fast5File`): FAST5 file
                pointer object.
            read_id (str): Read identifier to extract.
            extract_dacs (bool): Extract raw DAC values.

        Returns:
            backends.SIGNAL_DATA and backends.SEQ_SUMM_INFO namedtuples
        """
        read = fast5_fp.get_read(read_id)
        channel_info = read.get_channel_info()
        seq_summ_info = extract_seq_summary_info(read, channel_info)
        dacs = scale_params = raw_sig = None
        if extract_dacs:
            # if not processing signal mappings, don't save dacs
            dacs = fast5_io.get_signal(read, scale=False)
            # scale parameters and trimming computed by guppy
            if self.model_type != PYGUPPY_NAME:
                med, mad = mh.med_mad(dacs)
                raw_sig = (dacs - med) / mad
                # scale_params are relative to current
                rd_factor = channel_info[mh.CHAN_INFO_RANGE] / \
                    channel_info[mh.CHAN_INFO_DIGI]
                scale_params = (
                    (med + channel_info[mh.CHAN_INFO_OFFSET]) * rd_factor,
                    mad * rd_factor)

        if self.model_type == TAI_NAME:
            if raw_sig is None:
                raw_sig = fast5_io.get_signal(read, scale=True)
            sig_data = SIGNAL_DATA(
                raw_signal=raw_sig, dacs=dacs, scale_params=scale_params,
                raw_len=raw_sig.shape[0], fast5_fn=fast5_fp.filename,
                read_id=read_id, stride=self.stride)
            return sig_data, seq_summ_info
        elif self.model_type == FAST5_NAME:
            bc_mod_post = fast5_io.get_posteriors(read)
            if extract_dacs:
                trim_start, trim_len = fast5_io.get_signal_trim_coordiates(
                    read)
                dacs = dacs[trim_start:trim_start + trim_len]
            sig_data = SIGNAL_DATA(
                raw_len=bc_mod_post.shape[0] * self.stride, dacs=dacs,
                scale_params=scale_params, fast5_fn=fast5_fp.filename,
                read_id=read_id,
                stride=self.stride, posteriors=bc_mod_post)
            return sig_data, seq_summ_info
        elif self.model_type == PYGUPPY_NAME:
            if dacs is None:
                dacs = fast5_io.get_signal(read, scale=False)
            sig_data = SIGNAL_DATA(
                dacs=dacs, raw_len=dacs.shape[0], fast5_fn=fast5_fp.filename,
                read_id=read_id, stride=self.stride,
                channel_info=channel_info)
            return sig_data, seq_summ_info

        raise mh.MegaError('Invalid model type')

    def _run_taiyaki_model(self, raw_sig, n_can_state=None):
        if self.model_type != TAI_NAME:
            raise mh.MegaError(
                'Attempted to run taiyaki model with non-taiyaki ' +
                'initialization.')
        try:
            trans_weights = self.tai_run_model(
                raw_sig, self.model, self.params.taiyaki.chunk_size,
                self.params.taiyaki.chunk_overlap,
                self.params.taiyaki.max_concur_chunks)
        except AttributeError:
            raise mh.MegaError('Out of date or incompatible model')
        except RuntimeError as e:
            LOGGER.debug('Likely out of memory error: {}'.format(str(e)))
            raise mh.MegaError(
                'Likely out of memory error. See log for details.')
        if self.device != self.torch.device('cpu'):
            self.torch.cuda.empty_cache()
        if n_can_state is not None:
            trans_weights = (
                np.ascontiguousarray(trans_weights[:, :n_can_state]),
                np.ascontiguousarray(trans_weights[:, n_can_state:]))

        return trans_weights

    def _softmax_mod_weights(self, raw_mod_weights):
        mod_layers = []
        for lab_indices in self.can_raw_mod_indices:
            if lab_indices is None:
                mod_layers.append(
                    np.ones((raw_mod_weights.shape[0], 1), dtype=np.float32))
            else:
                mod_layers.append(_log_softmax_axis1(
                    raw_mod_weights[:, lab_indices]))
        return np.concatenate(mod_layers, axis=1)

    def pyguppy_basecall(self, client, reads_batch, failed_reads_q=None):
        def do_sleep():
            # function for profiling purposes
            sleep(PYGUPPY_PER_TRY_TIMEOUT)

        def get_completed_reads():
            do_retry = True
            n_reconnect_attempts = 0
            comp_reads = err_str = None
            while do_retry:
                try:
                    comp_reads = client.get_completed_reads()
                    do_retry = False
                except ConnectionError as e:
                    if n_reconnect_attempts < PYGUPPY_MAX_RECONNECT_ATTEMPTS:
                        LOGGER.debug('Reconnecting to server (get reads)')
                        client.connect()
                        n_reconnect_attempts += 1
                    else:
                        err_str = ('Pyguppy pass read connection error ' +
                                   '"{}"').format(str(e))
                        do_retry = False
                    err_str = ('Pyguppy get completed reads connection ' +
                               'error "{}"').format(str(e))
                except RuntimeError as e:
                    err_str = ('Pyguppy get completed reads invalid error ' +
                               '"{}"').format(str(e))
                    do_retry = False
            return comp_reads, err_str

        saved_input_data = {}
        completed_reads = []
        for sig_info, seq_summ_info in reads_batch:
            err_str = None
            read_sent = False
            n_reconnect_attempts = 0
            pyguppy_read = get_pyguppy_read(
                sig_info.read_id, sig_info.dacs, sig_info.channel_info)
            while not read_sent:
                try:
                    read_sent = client.pass_read(pyguppy_read)
                except ValueError as e:
                    err_str = ('Pyguppy pass read malformed error ' +
                               '"{}"').format(str(e))
                    read_sent = True
                except ConnectionError as e:
                    # attempt to reconnect to server when connection error
                    # occurs
                    if n_reconnect_attempts < PYGUPPY_MAX_RECONNECT_ATTEMPTS:
                        LOGGER.debug('Reconnecting to server (pass read)')
                        client.connect()
                        n_reconnect_attempts += 1
                    else:
                        err_str = ('Pyguppy pass read connection error ' +
                                   '"{}"').format(str(e))
                        read_sent = True
                except RuntimeError as e:
                    err_str = ('Pyguppy pass read undefined error ' +
                               '"{}"').format(str(e))
                    read_sent = True
                # get completed reads while sending reads so server doesn't
                # back up indefinitely
                iter_comp_reads = get_completed_reads()[0]
                if iter_comp_reads is not None:
                    completed_reads.extend(iter_comp_reads)

            if err_str is None:
                saved_input_data[sig_info.read_id] = (sig_info, seq_summ_info)
            else:
                raw_len = sig_info.raw_len \
                    if hasattr(sig_info, 'raw_len') else 0
                fn_rid = '{}:::{}'.format(sig_info.fast5_fn, sig_info.read_id)
                if failed_reads_q is not None:
                    failed_reads_q.put(tuple(mh.READ_STATUS(
                        is_err=True, do_update_prog=True, err_type=err_str,
                        fast5_fn=fn_rid, n_sig=raw_len)))
                LOGGER.debug('{} BasecallingFailed "{}"'.format(
                    sig_info.read_id, err_str))

        # yield reads that have been called already
        for called_read in completed_reads:
            read_id = called_read["metadata"]["read_id"]
            try:
                sig_info, seq_summ_info = saved_input_data[read_id]
            except KeyError:
                # read submitted in last batch now finished
                LOGGER.debug('{} timeout read finished'.format(read_id))
                continue
            LOGGER.debug('{} BasecallingCompleted'.format(read_id))
            yield (parse_pyguppy_called_read(called_read), sig_info,
                   seq_summ_info)
            del saved_input_data[read_id]

        # process reads until timeout for batch
        err_str = 'Guppy server timeout (see --guppy-timeout argument)'
        for _ in range(self.pyguppy_retries):
            completed_reads, comp_reads_err = get_completed_reads()
            if comp_reads_err is not None:
                err_str = comp_reads_err
                break
            for called_read in completed_reads:
                read_id = called_read["metadata"]["read_id"]
                try:
                    sig_info, seq_summ_info = saved_input_data[read_id]
                except KeyError:
                    # read submitted in last batch now finished
                    LOGGER.debug('{} timeout read finished'.format(read_id))
                    continue
                LOGGER.debug('{} BasecallingCompleted'.format(read_id))
                yield (parse_pyguppy_called_read(called_read), sig_info,
                       seq_summ_info)
                del saved_input_data[read_id]
            # if all reads have been processed break from loop
            if len(saved_input_data) == 0:
                err_str = None
                break
            do_sleep()

        # if there are any left over reads report timeout errors
        for read_id, (sig_info, _) in saved_input_data.items():
            raw_len = sig_info.raw_len if hasattr(sig_info, 'raw_len') else 0
            fn_rid = '{}:::{}'.format(sig_info.fast5_fn, read_id)
            if failed_reads_q is not None:
                failed_reads_q.put(tuple(mh.READ_STATUS(
                    is_err=True, do_update_prog=True, fast5_fn=fn_rid,
                    n_sig=raw_len, err_type=err_str)))
            LOGGER.debug('{} BasecallingFailed {}'.format(read_id, err_str))

    def _postprocess_pyguppy_called_read(
            self, called_read, sig_info, seq_summ_info, return_post_w_mods,
            return_mod_scores, update_sig_info, signal_reversed,
            mod_bc_min_prob):
        # compute run length cumsum from move table
        rl_cumsum = np.where(called_read.move)[0]
        rl_cumsum = np.insert(rl_cumsum, rl_cumsum.shape[0],
                              called_read.move.shape[0])

        post_w_mods = mods_scores = None
        if self.is_cat_mod:
            # split canonical posteriors and mod transition weights
            can_post = np.ascontiguousarray(
                called_read.state[:, :self.n_can_state])
            if return_mod_scores or return_post_w_mods:
                mods_weights = self._softmax_mod_weights(
                    called_read.state[:, self.n_can_state:])
                if return_post_w_mods:
                    post_w_mods = np.concatenate(
                        [can_post, mods_weights], axis=1)
                if return_mod_scores:
                    mods_scores = mods_weights[rl_cumsum[:-1]]
                    if signal_reversed:
                        mods_scores = mods_scores[::-1]
                    mods_scores = self.format_mod_scores(
                        called_read.seq, mods_scores, mod_bc_min_prob)
        else:
            can_post = called_read.state

        # check validity of pyguppy results
        if len(called_read.seq) != len(called_read.qual) or \
           len(called_read.seq) != rl_cumsum.shape[0] - 1:
            LOGGER.debug((
                'Invalid results recieved from pyguppy backend: ' +
                '{}\t{}\t').format(len(called_read.seq), len(called_read.qual),
                                   rl_cumsum.shape[0]))
            raise mh.MegaError(
                'Invalid results recieved from pyguppy backend.')

        if update_sig_info:
            # add scale_params and trimmed dacs to sig_info
            trimmed_dacs = sig_info.dacs[called_read.trimmed_samples:]
            # guppy does not apply the med norm factor
            scale_params = (
                called_read.scaling_shift,
                called_read.scaling_scale * mh.MED_NORM_FACTOR)
            sig_info = sig_info._replace(
                raw_len=trimmed_dacs.shape[0], dacs=trimmed_dacs,
                scale_params=scale_params)

        if signal_reversed:
            called_read = called_read._replace(
                seq=called_read.seq[::-1], qual=called_read.qual[::-1])

        # update seq summary info with basecalling info
        try:
            samp_rate = sig_info.channel_info[mh.CHAN_INFO_SAMP_RATE]
            try:
                tmplt_start = '{:.6f}'.format(
                    float(seq_summ_info.start_time) + (
                        called_read.trimmed_samples / samp_rate))
            except ValueError:
                tmplt_start = seq_summ_info.start_time
            tmplt_dur = '{:.6f}'.format(
                (sig_info.dacs.shape[0] - called_read.trimmed_samples) /
                samp_rate)
            seq_len = len(called_read.seq)
            mean_q_score = '{:.6f}'.format(mh.get_mean_q_score(
                called_read.qual))
            med = '{:.6f}'.format(called_read.scaling_shift)
            mad = '{:.6f}'.format(called_read.scaling_scale)
            seq_summ_info = seq_summ_info._replace(
                template_start=tmplt_start, template_duration=tmplt_dur,
                sequence_length_template=seq_len,
                mean_qscore_template=mean_q_score, median_template=med,
                mad_template=mad)
        except Exception:
            # if anything goes wrong don't let it fail the read
            pass

        return (sig_info, seq_summ_info, called_read, rl_cumsum, can_post,
                post_w_mods, mods_scores)

    def _run_pyguppy_backend(
            self, reads_batch, return_post_w_mods, return_mod_scores,
            update_sig_info, signal_reversed, mod_bc_min_prob, failed_reads_q):
        if self.model_type != PYGUPPY_NAME:
            raise mh.MegaError(
                'Attempted to run pyguppy model with non-pyguppy ' +
                'initialization.')

        for called_read, sig_info, seq_summ_info in self.pyguppy_basecall(
                self.client, reads_batch, failed_reads_q):
            try:
                yield self._postprocess_pyguppy_called_read(
                    called_read, sig_info, seq_summ_info, return_post_w_mods,
                    return_mod_scores, update_sig_info, signal_reversed,
                    mod_bc_min_prob)
            # only catch Megalodon errors here, all others caught upstream
            except mh.MegaError as e:
                raw_len = sig_info.raw_len \
                    if hasattr(sig_info, 'raw_len') else 0
                fn_rid = '{}:::{}'.format(sig_info.fast5_fn, sig_info.read_id)
                failed_reads_q.put(tuple(mh.READ_STATUS(
                    is_err=True, do_update_prog=True, err_type=str(e),
                    fast5_fn=fn_rid, n_sig=raw_len)))
                LOGGER.debug('{} Failed {}'.format(sig_info.read_id, str(e)))

    def _run_non_pyguppy_backend(
            self, sig_info, seq_summ_info, return_post_w_mods=True,
            return_mod_scores=False, update_sig_info=False,
            signal_reversed=False, mod_bc_min_prob=mh.DEFAULT_MOD_MIN_PROB):
        post_w_mods = mod_weights = None
        if self.model_type == TAI_NAME:
            # run neural network with taiyaki
            if self.is_cat_mod:
                bc_weights, mod_weights = self._run_taiyaki_model(
                    sig_info.raw_signal, self.n_can_state)
            else:
                bc_weights = self._run_taiyaki_model(sig_info.raw_signal)
            # perform forward-backward algorithm on neural net output
            can_post = decode.crf_flipflop_trans_post(bc_weights, log=True)
            if return_post_w_mods and self.is_cat_mod:
                post_w_mods = np.concatenate([can_post, mod_weights], axis=1)
            # set mod_weights to None if mods_scores not requested to
            # avoid extra computation
            if not return_mod_scores:
                mod_weights = None
        else:
            # FAST5 stored posteriors backend
            if self.is_cat_mod:
                # split canonical posteriors and mod transition weights
                # producing desired return arrays
                can_post = np.ascontiguousarray(
                    sig_info.posteriors[:, :self.n_can_state])
                if return_mod_scores or return_post_w_mods:
                    # convert raw neural net mod weights to softmax weights
                    mod_weights = self._softmax_mod_weights(
                        sig_info.posteriors[:, self.n_can_state:])
                    if return_post_w_mods:
                        post_w_mods = np.concatenate(
                            [can_post, mod_weights], axis=1)
                    if not return_mod_scores:
                        mod_weights = None
            else:
                can_post = sig_info.posteriors

        # decode posteriors to sequence and per-base mod scores
        r_seq, _, rl_cumsum, mods_scores = decode.decode_post(
            can_post, self.can_alphabet, mod_weights, self.can_nmods)
        if return_mod_scores:
            mods_scores = self.format_mod_scores(
                r_seq, mods_scores, mod_bc_min_prob)
        # TODO implement quality extraction for taiyaki and fast5 modes
        # and add mean_qscore_template to seq summary
        r_qual = None

        if seq_summ_info is not None:
            try:
                # update seq summary info with basecalling info
                seq_summ_info = seq_summ_info._replace(
                    template_start=seq_summ_info.start_time,
                    template_duration='{:.6f}'.format(
                        sig_info.dacs.shape[0] /
                        sig_info.channel_info[mh.CHAN_INFO_SAMP_RATE]),
                    sequence_length_template=len(r_seq),
                    median_template='{:.4f}'.format(sig_info.scale_params[0]),
                    mad_template='{:.4f}'.format(sig_info.scale_params[1]))
            except Exception:
                pass

        return (sig_info, seq_summ_info, r_seq, r_qual, rl_cumsum, can_post,
                post_w_mods, mods_scores)

    def iter_basecalled_reads(
            self, reads_batch, return_post_w_mods=True,
            return_mod_scores=False, update_sig_info=False,
            signal_reversed=False, mod_bc_min_prob=mh.DEFAULT_MOD_MIN_PROB,
            failed_reads_q=None):
        if self.model_type not in (TAI_NAME, FAST5_NAME, PYGUPPY_NAME):
            raise mh.MegaError('Invalid model backend')

        if self.model_type == PYGUPPY_NAME:
            for bc_res in self._run_pyguppy_backend(
                    reads_batch, return_post_w_mods, return_mod_scores,
                    update_sig_info, signal_reversed, mod_bc_min_prob,
                    failed_reads_q):
                yield bc_res
        else:
            for sig_info, seq_summ_info in reads_batch:
                try:
                    yield self._run_non_pyguppy_backend(
                        sig_info, seq_summ_info, return_post_w_mods,
                        return_mod_scores, update_sig_info, signal_reversed,
                        mod_bc_min_prob)
                # only catch Megalodon errors here, all others caught upstream
                except mh.MegaError as e:
                    if failed_reads_q is not None:
                        raw_len = sig_info.raw_len \
                            if hasattr(sig_info, 'raw_len') else 0
                        fn_rid = '{}:::{}'.format(
                            sig_info.fast5_fn, sig_info.read_id)
                        failed_reads_q.put(tuple(mh.READ_STATUS(
                            is_err=True, do_update_prog=True, err_type=str(e),
                            fast5_fn=fn_rid, n_sig=raw_len)))
                    LOGGER.debug('{} Failed {}'.format(
                        sig_info.read_id, str(e)))

    def close(self):
        if self.model_type == PYGUPPY_NAME:
            self.guppy_server_proc.terminate()
            self.guppy_out_fp.close()
            self.guppy_err_fp.close()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
