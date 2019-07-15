import sys
import string
import numpy as np
from collections import defaultdict, namedtuple

from megalodon import logging, megalodon_helper as mh


# model type specific information
TAI_NAME = 'taiyaki'

CAN_ALPHABET = 'ACGT'
MOD_ALPHABET = CAN_ALPHABET + ''.join(
    b for b in string.ascii_uppercase[::-1] if b not in CAN_ALPHABET)


class ModelInfo(object):
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
        base_proc_per_device = int(np.ceil(self.num_proc / len(self.devices)))
        procs_per_device = np.repeat(
            base_proc_per_device, len(self.devices))
        if base_proc_per_device * len(self.devices) > self.num_proc:
            procs_per_device[
                -(base_proc_per_device * len(self.devices) -
                  self.num_proc):] -= 1
        assert sum(procs_per_device) == self.num_proc
        self.process_devices = [
            dv for dv, n_dv in zip(self.devices, procs_per_device)
            for _ in range(n_dv)]

        try:
            # import modules
            from taiyaki.helpers import load_model as load_taiyaki_model
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
        self.output_size = ff_layer.size
        if self.is_cat_mod:
            # Modified base model is defined by 3 fixed fields in taiyaki
            # can_nmods, output_alphabet and modified_base_long_names
            self.output_alphabet = ff_layer.output_alphabet
            self.can_nmods = ff_layer.can_nmods
            self.ordered_mod_long_names = ff_layer.ordered_mod_long_names

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
        else:
            if mh.nstate_to_nbase(ff_layer.size) != 4:
                raise NotImplementedError(
                    'Naive modified base flip-flop models are not ' +
                    'supported.')
            self.output_alphabet = mh.ALPHABET
            self.can_alphabet = mh.ALPHABET
            self.mod_long_names = []
            self.str_to_int_mod_labels = {}
        self.n_mods = len(self.mod_long_names)

        return

    def __init__(
            self, taiyaki_model_fn, devices=None,
            num_proc=1, chunk_size=None, chunk_overlap=None,
            max_concur_chunks=None):
        self.taiyaki_model_fn = taiyaki_model_fn
        self.devices = devices
        self.num_proc = num_proc
        self.chunk_size = chunk_size
        self.chunk_overlap = chunk_overlap
        self.max_concur_chunks = max_concur_chunks

        if self.taiyaki_model_fn is not None:
            self._load_taiyaki_model()
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
                self.device = self.torch.device(device)
                self.torch.cuda.set_device(self.device)
                self.model = self.model.cuda()
            self.model = self.model.eval()

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
            except RuntimeError:
                raise mh.MegaError('Likely out of memory error.')
            if self.device != self.torch.device('cpu'):
                self.torch.cuda.empty_cache()
            if n_can_state is not None:
                trans_weights = (
                    np.ascontiguousarray(trans_weights[:,:n_can_state]),
                    np.ascontiguousarray(trans_weights[:,n_can_state:]))
        else:
            raise mh.MegaError('Invalid model type.')

        return trans_weights


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
