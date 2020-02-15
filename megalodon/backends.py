import re
import sys
from time import sleep
from collections import defaultdict, namedtuple

import numpy as np

from megalodon import decode, fast5_io, logging, megalodon_helper as mh


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
    'available', 'host', 'port', 'timeout'))
PYGUPPY_PARAMS.__new__.__defaults__ = tuple([None, ] * 3)
BACKEND_PARAMS = namedtuple('BACKEND_PARAMS', (
    TAI_NAME, FAST5_NAME, PYGUPPY_NAME))

COMPAT_GUPPY_MODEL_TYPES = set(('flipflop',))
DEFAULT_GUPPY_HOST = 'localhost'
DEFAULT_GUPPY_PORT = 5555
DEFAULT_GUPPY_TIMEOUT = 5.0
PYGUPPY_PER_TRY_TIMEOUT = 0.1

# maximum time (in seconds) to wait before assigning device
# over different processes. Actual time choosen randomly
# as many simultaneous calls sometimes conflict and kill
# some device assignments
MAX_DEVICE_WAIT = 1.0

SIGNAL_DATA = namedtuple('SIGNAL_DATA', (
    'fast5_fn', 'read_id', 'raw_len', 'dacs', 'raw_signal',
    'scale_params', 'stride', 'posteriors'))
# set default value of None for ref, alts and ref_start
SIGNAL_DATA.__new__.__defaults__ = tuple([None, ] * 5)


def parse_backend_params(args, num_fast5_startup_reads=5):
    # TOOD try to import high level packages to ensure available
    # parse taiyaki backend params
    if any(not hasattr(args, k) for k in (
            'chunk_size', 'chunk_overlap', 'max_concurrent_chunks',
            'load_default_taiyaki_model', 'taiyaki_model_filename')):
        tai_params = TAI_PARAMS(False)
    else:
        tai_model_fn = mh.get_model_fn(
            args.taiyaki_model_filename, args.load_default_taiyaki_model)
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
            'guppy_server_host', 'guppy_server_port',
            'do_not_use_guppy_server')) or args.do_not_use_guppy_server:
        pyguppy_params = PYGUPPY_PARAMS(False)
    else:
        pyguppy_params = PYGUPPY_PARAMS(
            True, args.guppy_server_host, args.guppy_server_port,
            args.guppy_timeout)

    return BACKEND_PARAMS(tai_params, fast5_params, pyguppy_params)


def _log_softmax_axis1(x):
    """ Compute log softmax over axis=1
    """
    e_x = np.exp((x.T - np.max(x, axis=1)).T)
    return np.log((e_x.T / e_x.sum(axis=1)).T)


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


class ModelInfo(object):
    def compute_mod_alphabet_attrs(self):
        # parse these values to more user-friendly data structures
        self.can_alphabet = ''
        self.can_indices = []
        self.mod_long_names = []
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

            curr_can_offset += can_base_nmods + 1
            curr_nmods += can_base_nmods

        self.can_indices.append(curr_can_offset)
        self.can_indices = np.array(self.can_indices).astype(np.uintp)
        self.can_base_mods = dict(self.can_base_mods)

    def _load_taiyaki_model(self):
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
            logger = logging.get_logger()
            logger.error(
                'Failed to import taiyaki. Ensure working ' +
                'installations to run megalodon')
            sys.exit(1)
        try:
            import torch
        except ImportError:
            logger = logging.get_logger()
            logger.error(
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
            self.compute_mod_alphabet_attrs()
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

    def _parse_minimal_alphabet_info(self):
        # compute values required for standard model attribute extraction
        self.n_mods = len(self.ordered_mod_long_names)
        self.is_cat_mod = self.n_mods > 0
        self.can_alphabet = None
        for v_alphabet in mh.VALID_ALPHABETS:
            if all(b in self.output_alphabet for b in v_alphabet):
                self.can_alphabet = v_alphabet
                break
        if self.can_alphabet is None:
            logger = logging.get_logger()
            logger.error(
                'Model information from FAST5 files contains invalid ' +
                'alphabet ({})'.format(self.output_alphabet))
            raise mh.MegaError('Invalid alphabet.')
        # compute number of modified bases for each canonical base

        if self.is_cat_mod:
            self.can_nmods = np.diff(np.array(
                [self.output_alphabet.index(b) for b in self.can_alphabet] +
                [len(self.output_alphabet), ])) - 1
            self.compute_mod_alphabet_attrs()
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
            if mh.nstate_to_nbase(self.output_size) != 4:
                raise NotImplementedError(
                    'Naive modified base flip-flop models are not ' +
                    'supported.')
            self.str_to_int_mod_labels = {}
            self.mod_long_names = []
            self.can_nmods = None

    def _load_fast5_post_out(self):
        def get_model_info_from_fast5(read):
            try:
                stride = fast5_io.get_stride(read)
                mod_long_names, out_alphabet = fast5_io.get_mod_base_info(read)
                out_size = fast5_io.get_posteriors(read).shape[1]
                mod_long_names = mod_long_names.split()
            except KeyError:
                logger.error(
                    'Fast5 read does not contain required attributes.')
                raise mh.MegaError(
                    'Fast5 read does not contain required attributes.')
            return stride, mod_long_names, out_alphabet, out_size

        logger = logging.get_logger()
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
            logger.error('No reads found.')
            raise mh.MegaError('No reads found.')

        for fast5_fn, read_id in read_iter:
            read = fast5_io.get_read(fast5_fn, read_id)
            r_s, r_omln, r_oa, r_os = get_model_info_from_fast5(read)
            if (
                    self.stride != r_s or
                    self.ordered_mod_long_names != r_omln or
                    self.output_alphabet != r_oa or
                    self.output_size != r_os):
                logger.error(
                    'Model information from FAST5 files is inconsistent. ' +
                    'Assure all reads were called with the same model.')
                raise mh.MegaError(
                    'Model information from FAST5 files is inconsistent.')
            nreads += 1
            if nreads >= self.params.fast5.num_startup_reads:
                break

        self._parse_minimal_alphabet_info()

    def _load_pyguppy(self, init_sig_len=40):
        self.model_type = PYGUPPY_NAME
        self.process_devices = [None, ] * self.num_proc

        from zmq.error import Again
        from pyguppyclient.decode import ReadData
        from pyguppyclient.client import GuppyBasecallerClient
        self.zmq_Again = Again
        self.pyguppy_ReadData = ReadData
        self.pyguppy_GuppyBasecallerClient = GuppyBasecallerClient

        self.pyguppy_retries = max(
            1, int(self.params.pyguppy.timeout / PYGUPPY_PER_TRY_TIMEOUT))

        init_client = self.pyguppy_GuppyBasecallerClient(
            None, host=self.params.pyguppy.host,
            port=self.params.pyguppy.port,
            timeout=PYGUPPY_PER_TRY_TIMEOUT,
            retries=self.pyguppy_retries)
        try:
            self.config_name = init_client.get_configs()[0].ConfigName()
        except IndexError:
            raise mh.MegaError('No config found from guppy server.')
        except self.zmq_Again:
            raise mh.MegaError('Unable to contact guppy server.')
        try:
            self.config_name = self.config_name.decode()
        except AttributeError:
            # already decoded
            pass
        init_client.disconnect()
        init_client = self.pyguppy_GuppyBasecallerClient(
            self.config_name, host=self.params.pyguppy.host,
            port=self.params.pyguppy.port,
            timeout=PYGUPPY_PER_TRY_TIMEOUT,
            retries=self.pyguppy_retries)
        try:
            init_client.connect()
        except self.zmq_Again:
            raise mh.MegaError(
                'Cannot connect to guppy server running on port {}.'.format(
                    self.params.pyguppy.port))
        try:
            init_read = init_client.basecall(
                ReadData(np.zeros(init_sig_len, dtype=np.int16), 'a'),
                state=True, trace=True)
        except TimeoutError:
            raise mh.MegaError(
                'Cannot connect to guppy server running on port {}.'.format(
                    self.params.pyguppy.port))
        init_client.disconnect()
        if init_read.state is None:
            raise mh.MegaError(
                'Guppy results do not contain posterior state. Ensure ' +
                '--post_out is set during guppy server initialization.')
        if init_read.model_type not in COMPAT_GUPPY_MODEL_TYPES:
            raise mh.MegaError(
                'Megalodon is not compatible with guppy model type: {}'.format(
                    init_read.model_type))

        self.stride = init_read.model_stride
        self.ordered_mod_long_names = init_read.mod_long_names
        self.output_alphabet = init_read.mod_alphabet
        self.output_size = init_read.state_size
        if self.ordered_mod_long_names is None:
            self.ordered_mod_long_names = []
        if self.output_alphabet is None:
            self.output_alphabet = mh.ALPHABET

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

        return

    @property
    def n_can_state(self):
        ncan_base = len(self.can_alphabet)
        return (ncan_base + ncan_base) * (ncan_base + 1)

    def prep_model_worker(self, device):
        """ Load model onto a newly spawned process
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
                    logger = logging.get_logger()
                    logger.error('Invalid CUDA device: {}'.format(device))
                    raise mh.MegaError('Error setting CUDA GPU device.')
            self.model = self.model.eval()
        elif self.model_type == PYGUPPY_NAME:
            # open guppy client interface (None indicates using config
            # from server)
            self.client = self.pyguppy_GuppyBasecallerClient(
                self.config_name, host=self.params.pyguppy.host,
                port=self.params.pyguppy.port, timeout=PYGUPPY_PER_TRY_TIMEOUT,
                retries=self.pyguppy_retries)
            self.client.connect()

        return

    def extract_signal_info(self, fast5_fn, read_id, extract_sig_map_info):
        read = fast5_io.get_read(fast5_fn, read_id)
        dacs = scale_params = raw_sig = None
        if extract_sig_map_info:
            # if not processing signal mappings, don't save dacs
            dacs = fast5_io.get_signal(read, scale=False)
            # scale parameters and trimming computed by guppy
            if not self.model_type == PYGUPPY_NAME:
                scale_params = mh.med_mad(dacs)
                raw_sig = (dacs - scale_params[0]) / scale_params[1]

        if self.model_type == TAI_NAME:
            if raw_sig is None:
                raw_sig = fast5_io.get_signal(read, scale=True)
            return SIGNAL_DATA(
                raw_signal=raw_sig, dacs=dacs, scale_params=scale_params,
                raw_len=raw_sig.shape[0], fast5_fn=fast5_fn, read_id=read_id,
                stride=self.stride)
        elif self.model_type == FAST5_NAME:
            bc_mod_post = fast5_io.get_posteriors(read)
            if extract_sig_map_info:
                trim_start, trim_len = fast5_io.get_signal_trim_coordiates(
                    read)
                dacs = dacs[trim_start:trim_start + trim_len]
            return SIGNAL_DATA(
                raw_len=bc_mod_post.shape[0] * self.stride, dacs=dacs,
                fast5_fn=fast5_fn, read_id=read_id, stride=self.stride,
                posteriors=bc_mod_post)
        elif self.model_type == PYGUPPY_NAME:
            if dacs is None:
                dacs = fast5_io.get_signal(read, scale=False)
            return SIGNAL_DATA(
                dacs=dacs, raw_len=dacs.shape[0], fast5_fn=fast5_fn,
                read_id=read_id, stride=self.stride)

        raise mh.MegaError('Invalid model type')

    def run_taiyaki_model(self, raw_sig, n_can_state=None):
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
            logging.get_logger().debug(
                'Likely out of memory error: {}'.format(str(e)))
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

    def run_pyguppy_model(
            self, sig_info, return_post_w_mods, return_mod_scores,
            update_sig_info):
        if self.model_type != PYGUPPY_NAME:
            raise mh.MegaError(
                'Attempted to run pyguppy model with non-pyguppy ' +
                'initialization.')

        post_w_mods = mods_scores = None
        try:
            called_read = self.client.basecall(
                self.pyguppy_ReadData(sig_info.dacs, sig_info.read_id),
                state=True, trace=True)
        except TimeoutError:
            raise mh.MegaError(
                'Pyguppy server timeout. See --guppy-timeout option')

        # compute rl_cumsum from move table
        rl_cumsum = np.where(called_read.move)[0]
        rl_cumsum = np.insert(rl_cumsum, rl_cumsum.shape[0],
                              called_read.move.shape[0])

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
                    # TODO apply np.NAN mask to scores not applicable to
                    # canonical basecalls
                    mods_scores = np.ascontiguousarray(
                        mods_weights[rl_cumsum[:-1]])
        else:
            can_post = called_read.state

        if update_sig_info:
            # add scale_params and trimmed dacs to sig_info
            trimmed_dacs = sig_info.dacs[called_read.trimmed_samples:]
            scale_params = (called_read.scaling['median'],
                            called_read.scaling['med_abs_dev'])
            sig_info = sig_info._replace(
                raw_len=trimmed_dacs.shape[0], dacs=trimmed_dacs,
                raw_signal=((trimmed_dacs - scale_params[0]) /
                            scale_params[1]).asttype(np.float32))

        return (called_read.seq, rl_cumsum, can_post, sig_info, post_w_mods,
                mods_scores)

    def basecall_read(
            self, sig_info, return_post_w_mods=True, return_mod_scores=False,
            update_sig_info=False):
        if self.model_type not in (TAI_NAME, FAST5_NAME, PYGUPPY_NAME):
            raise mh.MegaError('Invalid model backend')

        # decoding is performed within pyguppy server, so shortcurcuit return
        # here as other methods require megalodon decoding.
        if self.model_type == PYGUPPY_NAME:
            return self.run_pyguppy_model(
                sig_info, return_post_w_mods, return_mod_scores,
                update_sig_info)

        post_w_mods = mod_weights = None
        if self.model_type == TAI_NAME:
            # run neural network with taiyaki
            if self.is_cat_mod:
                bc_weights, mod_weights = self.run_taiyaki_model(
                    sig_info.raw_signal, self.n_can_state)
            else:
                bc_weights = self.run_taiyaki_model(sig_info.raw_signal)
            # perform forward-backward algorithm on neural net output
            can_post = decode.crf_flipflop_trans_post(bc_weights, log=True)
            if return_post_w_mods and self.is_cat_mod:
                post_w_mods = np.concatenate([can_post, mod_weights], axis=1)
            # set mod_weights to None if mod_scores not requested to
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

        return r_seq, rl_cumsum, can_post, sig_info, post_w_mods, mods_scores


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
