import string
import numpy as np
from collections import namedtuple

from megalodon import megalodon_helper as mh


# model type specific information
TAI_NAME = 'taiyaki'
FLP_NAME = 'flappie'

CAN_ALPHABET = 'ACGT'
MOD_ALPHABET = CAN_ALPHABET + ''.join(
    b for b in string.ascii_uppercase[::-1] if b not in CAN_ALPHABET)


class ModelInfo(object):
    def __init__(self, flappie_model_name, taiyaki_model_fn, device):
        if flappie_model_name is not None:
            self.model_type = FLP_NAME
            self.name = flappie_model_name

            # import module
            global flappy
            import flappy

            # compiled flappie models require hardcoding or running fake model
            if flappie_model_name in ('r941_cat_mod', 'r941_5mC'):
                if flappie_model_name == 'r941_cat_mod':
                    self.is_cat_mod = True
                    can_nmods, output_alphabet = flappy.get_cat_mods()
                    self.can_nmods = can_nmods
                    self.n_mods = self.can_nmods.sum()
                    self.output_size = 41 + self.n_mods
                    self.can_offsets = np.cumsum(np.concatenate(
                        [[0], self.can_nmods[:-1] + 1]))
                    can_bases = ''.join(
                        output_alphabet[b_i] for b_i in self.can_offsets)
                    mod_bases = ''.join(
                        output_alphabet[b_i + 1:b_i + 1 + b_nmods]
                        for b_i, b_nmods in
                        zip(self.can_offsets, self.can_nmods))
                    self.alphabet = can_bases + mod_bases
                    self.collapse_alphabet = can_bases + ''.join(
                        b for can_b, b_nmods in zip(can_bases, self.can_nmods)
                        for b in can_b * b_nmods)
                else:
                    self.is_cat_mod = False
                    self.n_mods = 1
                    self.output_size = 60
                    self.alphabet = MOD_ALPHABET[:5]
                    self.collapse_alphabet = CAN_ALPHABET + 'C'
            else:
                self.is_cat_mod = False
                self.alphabet = mh.ALPHABET
                self.collapse_alphabet = mh.ALPHABET
                self.output_size = 40
                self.n_mods = 0
        elif taiyaki_model_fn is not None:
            self.model_type = TAI_NAME
            self.fn = taiyaki_model_fn
            self.device = device

            # import modules
            global load_taiyaki_model, GlobalNormFlipFlopCatMod, torch, \
                tai_run_model
            from taiyaki.helpers import load_model as load_taiyaki_model
            from taiyaki.basecall_helpers import run_model as tai_run_model
            try:
                from taiyaki.layers import GlobalNormFlipFlopCatMod
            except ImportError:
                GlobalNormFlipFlopCatMod = None
            import torch

            tmp_model = load_taiyaki_model(taiyaki_model_fn)
            self.is_cat_mod = (
                GlobalNormFlipFlopCatMod is not None and isinstance(
                    tmp_model.sublayers[-1], GlobalNormFlipFlopCatMod))
            self.output_size = tmp_model.sublayers[-1].size
            try:
                # TODO Add mod long names
                self.alphabet = tmp_model.alphabet
                self.collapse_alphabet = tmp_model.collapse_alphabet
            except AttributeError:
                self.alphabet = mh.ALPHABET
                self.collapse_alphabet = mh.ALPHABET
            ncan_base = len(set(self.collapse_alphabet))
            self.n_mods = len(self.alphabet) - ncan_base
        else:
            raise mh.MegaError('Invalid model specification.')

        return

    def prep_model_worker(self):
        if self.model_type == TAI_NAME:
            # setup for taiyaki model
            self.model = load_taiyaki_model(self.fn)
            if self.device is not None:
                self.device = torch.device(self.device)
                torch.cuda.set_device(self.device)
                self.model = self.model.cuda()
            else:
                self.device = torch.device('cpu')
            self.model = self.model.eval()

        return

    def run_model(
            self, raw_sig, n_can_state=None, chunk_size=10000, overlap=200):
        if self.model_type == FLP_NAME:
            rt = flappy.RawTable(raw_sig)
            # flappy will return split bc and mods based on model
            trans_weights = flappy.run_network(rt, self.name)
        elif self.model_type == TAI_NAME:
            raw_sig_t = torch.from_numpy(raw_sig).to(self.device)
            try:
                trans_weights = tai_run_model(
                    raw_sig_t, self.model, chunk_size, overlap)
            except AttributeError:
                raise MegaError('Out of date or incompatible model')
            except RuntimeError:
                raise MegaError('Likely out of memory error.')
            if self.device != torch.device('cpu'):
                torch.cuda.empty_cache()
            if n_can_state is not None:
                trans_weights = (
                    np.ascontiguousarray(trans_weights[:,:n_can_state]),
                    np.ascontiguousarray(trans_weights[:,n_can_state:]))
        else:
            raise MegaError('Invalid model type.')

        return trans_weights


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `megalodon -h`')
    sys.exit(1)
