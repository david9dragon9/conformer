import os
import glob
import random
import numpy as np
from conformer.mask.base import MSAMasker
from conformer.utils.utils import dedup
from conformer.utils.af_cluster_utils import load_fasta
import torch
import torch.nn as nn
import torch.nn.functional as F
import copy


class RandomSelectedMSAMasker(MSAMasker):
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
            paths = self.modify_single(
                tag, input_dict[tag], all_ids, all_msas, data_dirs
            )
            random.shuffle(paths)
            path_results[tag] = paths
        return path_results

    def modify_single(self, tag, seq, input_ids, input_msa, data_dirs):
        input_msa = [
            "".join([x for x in s if x.isupper() or x == "-"]) for s in input_msa
        ]
        if self.config.residue_selection == "chosen":
            chosen_residues = self.config.chosen_residues
        elif self.config.residue_selection is None:
            chosen_residues = []
        else:
            raise ValueError(
                f"Unsupported residue selection: {self.config.residue_selection}"
            )

        if self.config.mask_mode == "preserve_residues":
            mask = 1 - torch.sum(
                F.one_hot(torch.as_tensor(chosen_residues), len(seq)), dim=0
            )
        elif self.config.mask_mode == "ignore_residues":
            mask = torch.sum(
                F.one_hot(torch.as_tensor(chosen_residues), len(seq)), dim=0
            )
        elif self.config.mask_mode == "ignore_all":
            mask = torch.ones(
                [
                    len(seq),
                ]
            )
        else:
            raise ValueError(f"Unsupported mask_mode: {self.config.mask_mode}")
        mask = mask.type(torch.float32).numpy()
        all_paths = []
        for i in range(self.config.num_samples):
            new_msa = copy.deepcopy(input_msa)
            for k in range(len(new_msa)):
                new_msa[k] = list(new_msa[k])
            random_selection = (
                np.random.random([len(input_msa), len(seq)]) < self.config.p
            ).astype(np.float32)
            final_mask = mask[None, :] * random_selection
            for pos in np.stack(np.where(final_mask), axis=-1).tolist():
                new_msa[pos[0]][pos[1]] = "-"
            for k in range(len(new_msa)):
                new_msa[k] = "".join(new_msa[k])
            curr_path = os.path.join(data_dirs["output_dir"], f"{tag}_{i}")
            all_paths.append(curr_path)
            os.makedirs(os.path.join(curr_path, tag), exist_ok=True)
            with open(os.path.join(curr_path, tag, f"{tag}_{i}.a3m"), "w") as f:
                for j in range(len(new_msa)):
                    f.write(f">{input_ids[j]}\n")
                    f.write(new_msa[j] + "\n")
        return all_paths
