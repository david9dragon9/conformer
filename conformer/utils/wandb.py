import wandb
import datetime
import os
import sys
from tools import smart_rmsd
from conformer.eval.download_structures import download_structure


def wandb_init(wandb_group, wandb_name, args):
    wandb_group = (
        wandb_group if wandb_group is not None else f"{args['p1']}-{args['p2']}"
    )
    wandb_name = (
        wandb_name
        if wandb_name is not None
        else datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    )
    wandb.init(
        project="conformer-OpenFold", group=wandb_group, name=wandb_name, config=args
    )


def wandb_log_results(new_path, args, structure_dir=None):
    os.makedirs(structure_dir, exist_ok=True)
    p1_path = (
        args["p1"]
        if os.path.exists(args["p1"])
        else os.path.join(structure_dir, f"{args['p1'][:4]}.pdb")
    )
    p2_path = (
        args["p2"]
        if os.path.exists(args["p2"])
        else os.path.join(structure_dir, f"{args['p2'][:4]}.pdb")
    )
    if not os.path.exists(p1_path):
        download_structure(args["p1"], structure_dir)
    if not os.path.exists(p2_path):
        download_structure(args["p2"], structure_dir)
    between_rmsd = smart_rmsd(
        p1_path, p2_path, chain1=args["p1"][5:], chain2=args["p2"][5:]
    )
    p1_rmsd = smart_rmsd(new_path, p1_path, chain1="A", chain2=args["p1"][5:])
    p2_rmsd = smart_rmsd(new_path, p2_path, chain1="A", chain2=args["p2"][5:])
    results = {"between_rmsd": between_rmsd, "p1_rmsd": p1_rmsd, "p2_rmsd": p2_rmsd}
    wandb.log(results)
    wandb.save(new_path)
    return results
