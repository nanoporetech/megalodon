import re
import sys
import string
import numpy as np
from time import sleep
from collections import defaultdict, namedtuple

from megalodon import decode, fast5_io, logging, megalodon_helper as mh


# model type specific information
TAI_NAME = 'taiyaki'
FAST5_NAME = 'fast5'

# maximum time (in seconds) to wait before assigning device
# over different processes. Actual time choosen randomly
# as many simultaneous calls sometimes conflict and kill
# some device assignments
MAX_DEVICE_WAIT = 1.0

SIGNAL_DATA = namedtuple('SIGNAL_DATA', (
    'fast5_fn', 'read_id', 'raw_len', 'dacs', 'raw_signal',
    'scale_params', 'stride', 'posteriors'))
# set default value of None for ref, alts and ref_start
SIGNAL_DATA.__new__.__defaults__ = tuple([None,] * 5)


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
        return

    def _load_taiyaki_model(self):
        if any(arg is None for arg in (
                self.chunk_size, self.chunk_overlap, self.max_concur_chunks)):
            logger = logging.get_logger()
            logger.debug(
                'Must provide chunk_size, chunk_overlap, ' +
                'max_concur_chunks in order to run the taiyaki ' +
                'base calling backend.')
        self.model_type = TAI_NAME

        if self.devices is None:
            self.devices = ['cpu',]
        self.process_devices = [parse_device(self.devices[i]) for i in np.tile(
            np.arange(len(self.devices)),
            (self.num_proc // len(self.devices)) + 1)][:self.num_proc]

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

        tmp_model = self.load_taiyaki_model(self.taiyaki_model_fn)
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
        self.n_mods = len(self.mod_long_names)

        return

    def _load_fast5_post_out(self):
        def get_model_info_from_fast5(read):
            try:
                stride = fast5_io.get_stride(read)
                mod_long_names, out_alphabet = fast5_io.get_mod_base_info(read)
                out_size = fast5_io.get_posteriors(read).shape[1]
                mod_long_names = mod_long_names.split()
            except KeyError:
                logger.error('Fast5 read does not contain required attributes.')
                raise mh.MegaError(
                    'Fast5 read does not contain required attributes.')
            return stride, mod_long_names, out_alphabet, out_size

        logger = logging.get_logger()
        self.model_type = FAST5_NAME
        self.process_devices = [None,] * self.num_proc

        read_iter = fast5_io.iterate_fast5_reads(self.fast5s_dir)
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
            if (self.stride != r_s or
                self.ordered_mod_long_names != r_omln or
                self.output_alphabet != r_oa or
                self.output_size != r_os):
                logger.error(
                    'Model information from FAST5 files is inconsistent. ' +
                    'Assure all reads were called with the same model.')
                raise mh.MegaError(
                    'Model information from FAST5 files is inconsistent.')
            nreads += 1
            if nreads >= self.num_startup_reads:
                break

        # compute values required for standard model attribute extraction
        self.n_mods = len(self.ordered_mod_long_names)
        self.is_cat_mod = self.n_mods > 0
        self.can_alphabet = None
        for v_alphabet in mh.VALID_ALPHABETS:
            if all(b in self.output_alphabet for b in v_alphabet):
                self.can_alphabet = v_alphabet
                break
        if self.can_alphabet is None:
            logger.error(
                'Model information from FAST5 files contains invalid ' +
                'alphabet ({})'.format(self.output_alphabet))
            raise mh.MegaError('Invalid alphabet.')
        # compute number of modified bases for each canonical base
        self.can_nmods = np.diff(np.array(
            [self.output_alphabet.index(b) for b in self.can_alphabet] +
            [len(self.output_alphabet),])) - 1

        if self.is_cat_mod:
            self.compute_mod_alphabet_attrs()
        else:
            if mh.nstate_to_nbase(self.output_size) != 4:
                raise NotImplementedError(
                    'Naive modified base flip-flop models are not ' +
                    'supported.')
            self.str_to_int_mod_labels = {}
            self.mod_long_names = []

        return

    def __init__(
            self, num_proc=1, fast5s_dir=None, num_startup_reads=5,
            taiyaki_model_fn=None, devices=None, chunk_size=None,
            chunk_overlap=None, max_concur_chunks=None):
        self.num_proc = num_proc

        # guppy posterior backend args
        self.fast5s_dir = fast5s_dir
        # number of reads to read in to identify model attributes
        self.num_startup_reads = num_startup_reads

        # taiyaki backend args
        self.taiyaki_model_fn = taiyaki_model_fn
        self.devices = devices
        self.chunk_size = chunk_size
        self.chunk_overlap = chunk_overlap
        self.max_concur_chunks = max_concur_chunks

        if self.taiyaki_model_fn is not None:
            self._load_taiyaki_model()
        elif self.fast5s_dir is not None:
            self._load_fast5_post_out()
        else:
            raise mh.MegaError('Invalid model specification.')

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
            self.model = self.load_taiyaki_model(self.taiyaki_model_fn)
            if device is None or device == 'cpu':
                self.device = self.torch.device('cpu')
            else:
                sleep(np.random.uniform(0, MAX_DEVICE_WAIT))
                self.device = self.torch.device(device)
                self.torch.cuda.set_device(self.device)
                self.model = self.model.cuda()
            self.model = self.model.eval()

        return

    def extract_signal_info(self, fast5_fn, read_id, extract_sig_map_info):
        read = fast5_io.get_read(fast5_fn, read_id)
        dacs = scale_params = raw_sig = None
        if extract_sig_map_info:
            # if not processing signal mappings, don't save dacs
            dacs = fast5_io.get_signal(read, scale=False)
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
                trim_start, trim_len = fast5_io.get_signal_trim_coordiates(read)
                dacs = dacs[trim_start:trim_start + trim_len]
            return SIGNAL_DATA(
                raw_len=bc_mod_post.shape[0] * self.stride, dacs=dacs,
                fast5_fn=fast5_fn, read_id=read_id, stride=self.stride,
                posteriors=bc_mod_post)

        raise mh.MegaError('Invalid model type')
        return

    def run_model(self, raw_sig, n_can_state=None):
        if self.model_type == TAI_NAME:
            if any(arg is None for arg in (
                    self.chunk_size, self.chunk_overlap,
                    self.max_concur_chunks)):
                logger = logging.get_logger()
                logger.error(
                    'Must provide chunk_size, chunk_overlap, ' +
                    'max_concur_chunks in order to run the taiyaki ' +
                    'base calling backend.')
            try:
                trans_weights = self.tai_run_model(
                    raw_sig, self.model, self.chunk_size, self.chunk_overlap,
                    self.max_concur_chunks)
            except AttributeError:
                raise mh.MegaError('Out of date or incompatible model')
            except RuntimeError as e:
                raise mh.MegaError(
                    'Likely out of memory error: {}'.format(str(e)))
            if self.device != self.torch.device('cpu'):
                self.torch.cuda.empty_cache()
            if n_can_state is not None:
                trans_weights = (
                    np.ascontiguousarray(trans_weights[:,:n_can_state]),
                    np.ascontiguousarray(trans_weights[:,n_can_state:]))
        else:
            raise mh.MegaError('Invalid model type.')

        return trans_weights

    def get_posteriors(self, sig_info):
        mod_weights, can_nmods = None, None
        if self.model_type == TAI_NAME:
            if self.is_cat_mod:
                bc_weights, mod_weights = self.run_model(
                    sig_info.raw_signal, self.n_can_state)
                can_nmods = self.can_nmods
            else:
                bc_weights = self.run_model(sig_info.raw_signal)
            bc_post = decode.crf_flipflop_trans_post(bc_weights, log=True)
        elif self.model_type == FAST5_NAME:
            bc_mod_post = sig_info.posteriors
            if self.is_cat_mod:
                bc_post, mod_weights = (
                    np.ascontiguousarray(bc_mod_post[:,:self.n_can_state]),
                    np.ascontiguousarray(bc_mod_post[:,self.n_can_state:]))
                can_nmods = self.can_nmods
            else:
                bc_post = bc_mod_post
        else:
            raise mh.MegaError('Invalid model type')

        return bc_post, mod_weights, can_nmods


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
