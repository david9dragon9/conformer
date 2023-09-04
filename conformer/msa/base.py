import numpy as np


class MSASubsampler:
    def __init__(self, config):
        self.config = config
        self.base_model = None

    def wrap_base_model(self, base_model):
        self.base_model = base_model
        return self

    def run(self, input_dict, data_dirs):
        path_results = self.modify(input_dict, data_dirs)
        all_output_paths = {}
        for tag, paths in path_results.items():
            print(tag, paths)
            curr_input_dict = {tag: input_dict[tag]}
            paths = np.random.choice(
                paths,
                size=[
                    self.config.num_samples,
                ],
                replace=False,
            ).tolist()
            for path in paths:
                data_dirs["output_dir"] = path
                data_dirs["alignment_dir"] = path
                print(curr_input_dict, data_dirs)
                output_paths = self.base_model.run(curr_input_dict, data_dirs)
                for key, value in output_paths.items():
                    if key not in all_output_paths:
                        all_output_paths[key] = []
                    all_output_paths[key].extend(value)
        return all_output_paths
