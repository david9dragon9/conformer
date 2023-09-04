import os
from conformer.utils.utils import load_json, load_pairs, dedup
from conformer.folder import Folder
import ast


class EvalInputs:
    def __init__(
        self,
        config,
        pair_path,
        output_dir,
        sequence_dict_path,
        alignment_dir=None,
        template_dir=None,
        structure_dir=None,
        baseline_outputs_dir=None,
        method_outputs_dir=None,
        leans_path=None,
        filtered_pair_path=None,
        method_output_paths=None,
    ):
        self.config = config
        if config.eval_method:
            self.method_folder = Folder(config.method)
        if config.eval_baseline:
            self.baseline_folder = Folder(config.baseline)
        self.pair_path = pair_path
        self.output_dir = output_dir
        self.sequence_dict_path = sequence_dict_path
        self.alignment_dir = (
            alignment_dir
            if alignment_dir is not None
            else os.path.join(output_dir, "alignments")
        )
        self.template_dir = (
            template_dir
            if template_dir is not None
            else os.path.join(output_dir, "templates")
        )
        self.structure_dir = (
            structure_dir
            if structure_dir is not None
            else os.path.join(output_dir, "structures")
        )
        self.baseline_outputs_dir = (
            baseline_outputs_dir
            if baseline_outputs_dir is not None
            else os.path.join(output_dir, "baseline_outputs")
        )
        self.method_outputs_dir = (
            method_outputs_dir
            if method_outputs_dir is not None
            else os.path.join(output_dir, "method_outputs")
        )
        self.leans_path = (
            leans_path
            if leans_path is not None
            else os.path.join(output_dir, "leans.json")
        )
        self.filtered_pair_path = (
            filtered_pair_path
            if filtered_pair_path is not None
            else os.path.join(output_dir, "filtered_pairs.txt")
        )
        if isinstance(method_output_paths, str):
            with open(method_output_paths, "r") as f:
                self.method_output_paths = ast.literal_eval(f.read())
        else:
            self.method_output_paths = method_output_paths

        self.sequence_dict = load_json(self.sequence_dict_path)
        self.leans = load_json(self.leans_path)
        self.pairs = load_pairs(self.pair_path)
        print(self.pair_path, self.pairs)
        self.all_tags = dedup(self.pairs) if self.pairs is not None else None
        self.filtered_pairs = load_pairs(self.filtered_pair_path)
