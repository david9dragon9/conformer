import sys

sys.path.append("EigenFold/")
sys.path.insert(0, "OmegaFoldE/")
from EigenFold.model import get_model
from EigenFold.utils.logging import get_logger
from EigenFold.utils.inference import inference_epoch
from EigenFold.utils.dataset import get_loader, get_args_suffix
from OmegaFoldE.omegafold.__main__ import OmegaFoldModel
import wandb
import datetime
from tools import smart_rmsd

logger = get_logger(__name__)

import pandas as pd
import os
import torch
import numpy as np
import tqdm
from argparse import Namespace
from conformer.eval.download_structures import download_structure

str_to_dtype = {
    "float16": torch.float16,
    "bfloat16": torch.bfloat16,
    "float32": torch.float32,
}


def make_omegafold_embeddings(config, input_dict, data_dirs):
    tags, sequences = [], []
    for key, value in input_dict.items():
        tags.append(key)
        sequences.append(value)
    splits = {
        "tags": tags,
        "name": tags,
        "seq": sequences,
        "seqlen": [len(seq) for seq in sequences],
        "seqres": sequences,
    }
    suffix = (
        ".".join(map(str, ["omegafold_num_recycling", config.omegafold_num_recycling]))
        + ".npz"
    )
    precision = str_to_dtype[config.precision]

    # load OmegaFold model
    omegafold = OmegaFoldModel(config.lm_weights_path, device=config.device)
    omegafold.model.to(precision)
    skipping = 0
    doing = 0
    for index, path in tqdm.tqdm(enumerate(splits["tags"])):
        embeddings_dir = os.path.join(data_dirs["output_dir"], "embeddings", path[:2])
        if not os.path.exists(embeddings_dir):
            os.makedirs(embeddings_dir)
        embeddings_path = os.path.join(embeddings_dir, path) + "." + suffix

        if os.path.exists(embeddings_path):
            skipping += 1
            continue

        doing += 1
        fasta_lines = [f">{path}", splits["seqres"][index]]

        edge_results, node_results = omegafold.inference(
            fasta_lines, config.omegafold_num_recycling  # , precision=precision
        )
        np.savez(
            embeddings_path,
            node_repr=node_results[0].to(torch.float32),
            edge_repr=edge_results[0].to(torch.float32),
        )
    print(splits, "DONE")
    print("Skipped", skipping)
    print("Done", doing)


def unwrap(config):
    x = config.__dict__["_content"]
    x = {k: v._val for k, v in x.items()}
    return x


class EigenFold:
    def __init__(self, config, full_config=None):
        self.config = config
        self.device = (
            self.config.device
            if self.config.device is not None
            else torch.device("cuda" if torch.cuda.is_available() else "cpu")
        )

        logger.info("Constructing model")
        model = get_model(self.config).to(self.device)
        ckpt = os.path.join(self.config.model_dir, self.config.ckpt)

        logger.info(f"Loading weights from {ckpt}")
        state_dict = torch.load(ckpt, map_location=torch.device("cpu"))
        model.load_state_dict(state_dict["model"], strict=True)
        self.model = model
        self.ep = state_dict["epoch"]

    def run(self, input_dict, data_dirs):
        if self.config.wandb:
            wandb.init(
                project="conformer-EigenFold",
                group="-".join(input_dict.keys()),
                name="-".join(input_dict.keys())
                + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                config=unwrap(self.config),
            )
        if self.config.embeddings_dir is None:
            self.config.embeddings_dir = os.path.join(
                data_dirs["output_dir"], "embeddings"
            )
        make_omegafold_embeddings(self.config, input_dict, data_dirs)
        tags, sequences = [], []
        for key, value in input_dict.items():
            tags.append(key)
            sequences.append(value)
        splits = {
            "name": tags,
            "seq": sequences,
            "seqlen": [len(seq) for seq in sequences],
            "seqres": sequences,
        }
        splits = pd.DataFrame(splits).set_index("name")
        val_loader = get_loader(
            Namespace(**unwrap(self.config)),
            None,
            splits,
            mode=self.config.split_key,
            shuffle=False,
        )
        samples, log = inference_epoch(
            Namespace(**unwrap(self.config)),
            self.model,
            val_loader.dataset,
            device=self.device,
            pdbs=True,
            elbo=self.config.elbo,
        )

        means = {key: np.mean(log[key]) for key in log if key != "path"}
        logger.info(f"Inference epoch {self.ep}: len {len(log['rmsd'])} MEANS {means}")

        # inf_name = f"{self.config.splits.split('/')[-1]}.ep{self.ep}.num{self.config.num_samples}.step{self.config.inf_step}.alpha{self.config.alpha}.beta{self.config.beta}"
        inf_name = f"{'.'.join(tags)}.ep{self.ep}.num{self.config.num_samples}.step{self.config.inf_step}.alpha{self.config.alpha}.beta{self.config.beta}"
        if self.config.inf_step != self.config.elbo_step:
            inf_name += f".elbo{self.config.elbo_step}"
        csv_path = os.path.join(data_dirs["output_dir"], f"{inf_name}.csv")
        pd.DataFrame(log).set_index("path").to_csv(csv_path)
        logger.info(f"Saved inf csv {csv_path}")

        if not os.path.exists(os.path.join(data_dirs["output_dir"], inf_name)):
            os.mkdir(os.path.join(data_dirs["output_dir"], inf_name))
        abrmsd, a1rmsd, a2rmsd = [], [], []
        all_output_paths = {}
        for i, samp in enumerate(samples):
            samp.pdb.write(
                os.path.join(
                    data_dirs["output_dir"],
                    inf_name,
                    samp.path.split("/")[-1] + f".{samp.copy}.anim.pdb",
                ),
                reverse=True,
            )
            curr_path = os.path.join(
                data_dirs["output_dir"],
                inf_name,
                samp.path.split("/")[-1] + f".{samp.copy}.pdb",
            )
            if samp.info.name not in all_output_paths:
                all_output_paths[samp.info.name] = []
            all_output_paths[samp.info.name].append(curr_path)
            samp.pdb.clear().add(samp.Y).write(
                os.path.join(
                    data_dirs["output_dir"],
                    inf_name,
                    samp.path.split("/")[-1] + f".{samp.copy}.pdb",
                )
            )
            if self.config.wandb:
                wandb.save(curr_path)
                os.makedirs(data_dirs["structure_dir"], exist_ok=True)
                p1_path = (
                    self.config.p1
                    if os.path.exists(self.config.p1)
                    else os.path.join(
                        data_dirs["structure_dir"], f"{self.config.p1[:4]}.pdb"
                    )
                )
                p2_path = (
                    self.config.p2
                    if os.path.exists(self.config.p2)
                    else os.path.join(
                        data_dirs["structure_dir"], f"{self.config.p2[:4]}.pdb"
                    )
                )
                if not os.path.exists(p1_path):
                    download_structure(self.config.p1, data_dirs["structure_dir"])
                if not os.path.exists(p2_path):
                    download_structure(self.config.p2, data_dirs["structure_dir"])
                between_rmsd = smart_rmsd(
                    p1_path,
                    p2_path,
                    chain1=self.config.p1[5:],
                    chain2=self.config.p2[5:],
                )
                p1_rmsd = smart_rmsd(
                    curr_path, p1_path, chain1="A", chain2=self.config.p1[5:]
                )
                p2_rmsd = smart_rmsd(
                    curr_path, p2_path, chain1="A", chain2=self.config.p2[5:]
                )
                wandb.log(
                    {
                        "e_between_rmsd": between_rmsd,
                        "e_p1_rmsd": p1_rmsd,
                        "e_p2_rmsd": p2_rmsd,
                        "sample_number": i,
                    }
                )
                abrmsd.append(between_rmsd)
                a1rmsd.append(p1_rmsd)
                a2rmsd.append(p2_rmsd)
                wandb.save(curr_path)
        if self.config.wandb:
            wandb.log(
                {
                    "between_rmsd": np.mean(abrmsd),
                    "p1_rmsd": np.mean(a1rmsd),
                    "p2_rmsd": np.mean(a2rmsd),
                }
            )
        wandb.finish()
        return all_output_paths
