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
                    self.cat_mods = flappy.get_cat_mods()
                    self.n_mods = self.cat_mods.sum()
                    self.is_cat_mod = True
                    self.output_size = 41 + self.n_mods
                    self.alphabet = MOD_ALPHABET[:4 + self.n_mods]
                    self.collapse_alphabet = CAN_ALPHABET + ''.join(
                        b for b_i, c in enumerate(self.cat_mods)
                        for b in [CAN_ALPHABET[b_i]] * c)
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
        elif taiyaki_model_fn is not None:
            self.model_type = TAI_NAME
            self.fn = taiyaki_model_fn
            self.device = device

            # import modules
            global load_taiyaki_model, GlobalNormFlipFlopCatMod, torch
            from taiyaki.helpers import load_model as load_taiyaki_model
            from taiyaki.layers import (
                GlobalNormFlipFlopCatMod, GlobalNormFlipFlopCatMod)
            import torch

            tmp_model = load_taiyaki_model(taiyaki_model_fn)
            self.is_cat_mod = isinstance(
                tmp_model.sublayers[-1], GlobalNormFlipFlopCatMod)
            self.output_size = tmp_model.sublayers[-1].size
            try:
                self.alphabet = tmp_model.alphabet
                self.collapse_alphabet = tmp_model.collapse_alphabet
            except AttributeError:
                self.alphabet = mh.ALPHABET
                self.collapse_alphabet = mh.ALPHABET
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

    def run_model(self, raw_sig, n_can_state=None):
        if self.model_type == FLP_NAME:
            rt = flappy.RawTable(raw_sig)
            # flappy will return split bc and mods based on model
            trans_weights = flappy.run_network(rt, self.name)
        elif self.model_type == TAI_NAME:
            with torch.no_grad():
                raw_sig_t = torch.from_numpy(
                    raw_sig[:,np.newaxis,np.newaxis]).to(self.device)
                try:
                    trans_weights = self.model(
                        raw_sig_t).squeeze().cpu().numpy()
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
