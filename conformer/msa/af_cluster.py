import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import os
import sys
import glob
import random
from conformer.utils.af_cluster_utils import load_fasta
from conformer.utils.cluster_msa import cluster_msa
from conformer.msa.base import MSASubsampler


class AFClusterMSASubsampler(MSASubsampler):
    def modify(self, input_dict, data_dirs):
        path_results = {}
        for tag in input_dict.keys():
            globbed_msas = glob.glob(
                os.path.join(data_dirs["alignment_dir"], tag) + "/*"
            )
            assert len(globbed_msas) > 0
            all_ids = []
            all_msas = []
            for msa_path in globbed_msas:
                loaded_ids, loaded_msa = load_fasta(msa_path)
                all_ids.extend(loaded_ids)
                all_msas.extend(loaded_msa)
            paths = self.modify_single(tag, all_ids, all_msas, data_dirs)
            random.shuffle(paths)
            path_results[tag] = paths
        return path_results

    def modify_single(self, tag, ids, msas, data_dirs):
        self.config.o = data_dirs["output_dir"]
        self.config.keyword = tag
        (
            all_output_paths,
            uniform_output_paths_10,
            uniform_output_paths_100,
        ) = cluster_msa(self.config, loaded_msa=(ids, msas))
        final_paths = []
        if self.config.use_all_output_paths:
            final_paths.extend(all_output_paths)
        if self.config.use_uniform_output_paths_10:
            final_paths.extend(uniform_output_paths_10)
        if self.config.use_uniform_output_paths_100:
            final_paths.extend(uniform_output_paths_100)
        return final_paths
