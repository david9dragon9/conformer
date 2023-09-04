import os
import ast


def run_method(eval_inputs):
    all_method_output_paths = {}
    for pair in eval_inputs.filtered_pairs:
        method_output_paths = eval_inputs.method_folder.run(
            fasta_dir=None,
            alignment_dir=eval_inputs.alignment_dir,
            template_dir=eval_inputs.template_dir,
            output_dir=eval_inputs.method_outputs_dir,
            structure_dir=eval_inputs.structure_dir,
            input_dict={pair[0]: eval_inputs.sequence_dict[pair[0]]},
        )
        for key, value in method_output_paths.items():
            if key not in all_method_output_paths:
                all_method_output_paths[key] = []
            all_method_output_paths[key].extend(value)
    eval_inputs.method_output_paths = all_method_output_paths
    with open(
        os.path.join(eval_inputs.method_outputs_dir, "method_output_paths.txt"), "w"
    ) as f:
        f.write(str(all_method_output_paths))
    return eval_inputs
