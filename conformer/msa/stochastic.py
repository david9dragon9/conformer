import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import os
import sys
import glob
import random
from conformer.utils.af_cluster_utils import load_fasta
from conformer.msa.base import MSASubsampler


class StochasticMSASubsampler(MSASubsampler):
    def __init__(self, config):
        self.config = config

    def modify(self, input_dict, data_dirs):
        path_results = {}
        for tag in input_dict.keys():
            globbed_msas = glob.glob(
                os.path.join(data_dirs["alignment_dir"], tag) + "/*"
            )
            assert len(globbed_msas) > 0
            all_msas = []
            for msa_path in globbed_msas:
                loaded_ids, loaded_msa = load_fasta(msa_path)
                all_msas.extend(loaded_msa)
            paths = self.modify_single(tag, all_msas, data_dirs)
            random.shuffle(paths)
            path_results[tag] = paths
        return path_results

    def modify_single(self, tag, input_msa, data_dirs):
        all_paths = []
        for i in range(self.config.num_samples):
            sampled_indices = np.random.choice(
                len(input_msa),
                size=[self.config.num_msa_per_sample],
                replace=self.config.replace_msa,
            )
            sampled_msa = [input_msa[si] for si in sampled_indices]
            curr_path = os.path.join(data_dirs["output_dir"], f"{tag}_{i}")
            all_paths.append(curr_path)
            os.makedirs(os.path.join(curr_path, tag), exist_ok=True)
            with open(os.path.join(curr_path, tag, f"{tag}_{i}.a3m"), "w") as f:
                for j in range(len(sampled_msa)):
                    f.write(f">s{sampled_indices[j]}_{j}\n")
                    f.write(sampled_msa[j] + "\n")
        return all_paths
