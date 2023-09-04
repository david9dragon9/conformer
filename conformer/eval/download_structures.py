import os


def download_structure(tag, structure_dir):
    os.system(
        f"wget -q https://files.rcsb.org/download/{tag[:4]}.pdb -P {structure_dir}"
    )


def download_structures(eval_inputs):
    os.makedirs(eval_inputs.structure_dir, exist_ok=True)
    for tag in eval_inputs.all_tags:
        download_structure(tag, eval_inputs.structure_dir)
    return eval_inputs
