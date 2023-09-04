from tools import smart_rmsd
from colorama import Fore
import os
import json


def evaluate_and_filter_pairs(eval_inputs):
    filtered_pairs = []
    all_leans = {}
    for pair in eval_inputs.pairs:
        a = smart_rmsd(
            os.path.join(eval_inputs.baseline_outputs_dir, f"{pair[0]}.pdb"),
            os.path.join(eval_inputs.structure_dir, f"{pair[0][:4]}.pdb"),
            chain1="A",
            chain2=pair[0][5:],
        )
        b = smart_rmsd(
            os.path.join(eval_inputs.baseline_outputs_dir, f"{pair[1]}.pdb"),
            os.path.join(eval_inputs.structure_dir, f"{pair[0][:4]}.pdb"),
            chain1="A",
            chain2=pair[0][5:],
        )

        c = smart_rmsd(
            os.path.join(eval_inputs.baseline_outputs_dir, f"{pair[0]}.pdb"),
            os.path.join(eval_inputs.structure_dir, f"{pair[1][:4]}.pdb"),
            chain1="A",
            chain2=pair[1][5:],
        )
        d = smart_rmsd(
            os.path.join(eval_inputs.baseline_outputs_dir, f"{pair[1]}.pdb"),
            os.path.join(eval_inputs.structure_dir, f"{pair[1][:4]}.pdb"),
            chain1="A",
            chain2=pair[1][5:],
        )

        if max(a, b) < min(c, d):
            all_leans[pair[0]] = 0
            filtered_pairs.append(pair)
        elif max(c, d) < min(a, b):
            all_leans[pair[0]] = 1
            filtered_pairs.append(pair)
    with open(eval_inputs.leans_path, "w") as f:
        json.dump(all_leans, f)
    with open(eval_inputs.filtered_pair_path, "w") as f:
        for pair in filtered_pairs:
            f.write(f"{pair}\n")
    eval_inputs.filtered_pairs = filtered_pairs
    eval_inputs.leans = all_leans
    return eval_inputs
