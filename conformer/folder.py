import torch
import torch.nn as nn
import torch.nn.functional as F

import numpy as np
import os
import sys
from conformer.base import get_base_model
from conformer.msa import get_msa_subsampler
from conformer.mask import get_masker
from conformer.data.load_data import load_fasta_sequences
import hydra


def yaml_replace(item):
    return (
        item.replace("None", "null").replace("False", "false").replace("True", "true")
    )


class Folder(nn.Module):
    def __init__(self, config, override_kwargs={}):
        super().__init__()
        if isinstance(config, str):
            overrides = [f"{key}={value}" for key, value in override_kwargs.items()]
            hydra.core.global_hydra.GlobalHydra.instance().clear()
            hydra.initialize("./configs/")
            self.config = hydra.compose(config, overrides=overrides)
        else:
            self.config = config
        print(config)
        self.base_model = get_base_model(self.config)
        if "msa" in self.config and self.config.msa is not None:
            msa_subsampler = get_msa_subsampler(self.config)
            self.base_model = msa_subsampler.wrap_base_model(self.base_model)
        if "mask" in self.config and self.config.mask is not None:
            masker = get_masker(self.config)
            self.base_model = masker.wrap_base_model(self.base_model)

    def run(
        self,
        fasta_dir,
        alignment_dir,
        template_dir,
        output_dir,
        structure_dir,
        input_dict=None,
    ):
        print(fasta_dir, alignment_dir, template_dir, output_dir, structure_dir)
        data_dirs = {
            "fasta_dir": fasta_dir,
            "alignment_dir": alignment_dir,
            "template_dir": template_dir,
            "output_dir": output_dir,
            "structure_dir": structure_dir,
        }
        if input_dict is None:
            input_dict = load_fasta_sequences(fasta_dir)
        output_paths = self.base_model.run(input_dict=input_dict, data_dirs=data_dirs)
        return output_paths
