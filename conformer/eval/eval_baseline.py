import json
import os


def eval_baseline(eval_inputs):
    input_dict = {
        pdb_id: eval_inputs.sequence_dict[pdb_id] for pdb_id in eval_inputs.all_tags
    }
    eval_inputs.baseline_folder.run(
        fasta_dir=None,
        alignment_dir=eval_inputs.alignment_dir,
        template_dir=eval_inputs.template_dir,
        output_dir=eval_inputs.baseline_outputs_dir,
        structure_dir=eval_inputs.structure_dir,
        input_dict=input_dict,
    )
    return eval_inputs
