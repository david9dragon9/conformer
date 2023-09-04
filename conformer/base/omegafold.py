import gc
import os
import sys
import argparse
import omegafold as of
import time
import torch
import logging
from omegafold import pipeline
from omegafold.pipeline import _load_weights
from omegafold.utils.protein_utils import residue_constants as rc
from omegafold import utils

# wget https://helixon.s3.amazonaws.com/release1.pt
def create_inputs(config, input_dict, data_dirs, deterministic=True):
    for i, (ch, fas) in enumerate(input_dict.items()):
        fas = fas.replace("Z", "E").replace("B", "D").replace("U", "C")
        aatype = torch.LongTensor(
            [rc.restypes_with_x.index(aa) if aa != "-" else 21 for aa in fas]
        )
        mask = torch.ones_like(aatype).float()
        assert torch.all(aatype.ge(0)) and torch.all(aatype.le(21)), (
            f"Only take 0-20 amino acids as inputs with unknown amino acid "
            f"indexed as 20"
        )
        name_max = 32
        if len(ch) < name_max:
            out_fname = ch.replace(os.path.sep, "-")
        else:
            out_fname = f"{i}th chain"
        out_fname = os.path.join(data_dirs["output_dir"], out_fname + ".pdb")

        num_res = len(aatype)
        data = list()
        g = None
        if deterministic:
            g = torch.Generator()
            g.manual_seed(num_res)
        for _ in range(config.num_cycle):
            p_msa = aatype[None, :].repeat(config.num_pseudo_msa, 1)
            p_msa_mask = torch.rand([config.num_pseudo_msa, num_res], generator=g).gt(
                config.pseudo_msa_mask_rate
            )
            p_msa_mask = torch.cat((mask[None, :], p_msa_mask), dim=0)
            p_msa = torch.cat((aatype[None, :], p_msa), dim=0)
            p_msa[~p_msa_mask.bool()] = 21
            data.append({"p_msa": p_msa, "p_msa_mask": p_msa_mask})
        yield utils.recursive_to(data, device=config.device), out_fname, ch


class OmegaFold:
    def __init__(self, config, full_config=None):
        self.config = config
        self.full_config = full_config
        self.forward_config = argparse.Namespace(
            subbatch_size=config.subbatch_size,
            num_recycle=config.num_cycle,
        )
        logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
        # create the output directory
        # get the model
        logging.info(f"Constructing OmegaFold")
        self.model = of.OmegaFold(of.make_config(config.model))
        weights_file = config.weights_file
        # if the output directory is not provided, we will create one alongside the
        # input fasta file
        if weights_file or config.weights_url:
            state_dict = _load_weights(config.weights_url, weights_file)
            state_dict = state_dict.pop("model", state_dict)
        else:
            state_dict = None
        if state_dict is None:
            logging.warning("Inferencing without loading weight")
        else:
            if "model" in state_dict:
                state_dict = state_dict.pop("model")
            self.model.load_state_dict(state_dict)
        self.model.eval()
        self.model.to(config.device)

    def run(self, input_dict, data_dirs):
        os.makedirs(data_dirs["output_dir"], exist_ok=True)
        # logging.info(f"Reading {args.input_file}")
        output_paths = {}
        for i, (input_data, save_path, tag) in enumerate(
            create_inputs(
                self.config,
                input_dict,
                data_dirs,
            )
        ):
            logging.info(f"Predicting {i + 1}th chain!")
            logging.info(f"{len(input_data[0]['p_msa'][0])} residues in this chain.")
            ts = time.time()
            try:
                output = self.model(
                    input_data,
                    predict_with_confidence=True,
                    fwd_cfg=self.forward_config,
                )
            except RuntimeError as e:
                logging.info(f"Failed to generate {save_path} due to {e}")
                logging.info(f"Skipping...")
                continue
            logging.info(f"Finished prediction in {time.time() - ts:.2f} seconds.")

            logging.info(f"Saving prediction to {save_path}")
            pipeline.save_pdb(
                pos14=output["final_atom_positions"],
                b_factors=output["confidence"] * 100,
                sequence=input_data[0]["p_msa"][0],
                mask=input_data[0]["p_msa_mask"][0],
                save_path=save_path,
                model=0,
            )
            output_paths[tag] = save_path
            logging.info(f"Saved")
            del output
            torch.cuda.empty_cache()
            gc.collect()
        logging.info("Done!")
        return output_paths
