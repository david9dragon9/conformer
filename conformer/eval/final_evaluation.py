from tools import smart_rmsd
import os
from tabulate import tabulate


def final_evaluation(eval_inputs):
    total = 0
    correct = 0
    for pair in eval_inputs.filtered_pairs:
        if eval_inputs.method_output_paths is None:
            output_paths = [
                os.path.join(eval_inputs.method_outputs_dir, f"{pair[0]}.pdb")
            ]
        else:
            output_paths = eval_inputs.method_output_paths[pair[0]]
        leans, others = [], []
        curr_correct = 0
        curr_total = 0
        for output_path in output_paths:
            a = smart_rmsd(
                output_path,
                os.path.join(eval_inputs.structure_dir, f"{pair[0][:4]}.pdb"),
                chain1="A",
                chain2=pair[0][5:],
            )
            b = smart_rmsd(
                output_path,
                os.path.join(eval_inputs.structure_dir, f"{pair[1][:4]}.pdb"),
                chain1="A",
                chain2=pair[1][5:],
            )
            print(
                tabulate(
                    [[a, b]],
                    headers=["LEANS", "OTHER"]
                    if eval_inputs.leans[pair[0]] == 0
                    else ["OTHER", "LEANS"],
                    tablefmt="simple_grid",
                )
            )
            # print(a, b)
            if eval_inputs.leans[pair[0]] == 0:
                leans.append(a)
                others.append(b)
                if b < a and b < 3.5:
                    curr_correct += 1
            elif eval_inputs.leans[pair[0]] == 1:
                leans.append(b)
                others.append(a)
                if a < b and a < 3.5:
                    curr_correct += 1
            else:
                assert False
            curr_total += 1
        if eval_inputs.config.eval_config == "best":
            total += int(curr_total > 0)
            correct += int(curr_correct > 0)
        elif eval_inputs.config.eval_config == "average":
            total += int(curr_total > 0)
            correct += curr_correct / curr_total
    print(f"FINAL ACCURACY: {correct / total}")
