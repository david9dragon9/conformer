import random
import torch
import numpy as np
import json
import ast
import os
from typing import Iterable
import omegaconf


def load_json(json_path):
    if os.path.exists(json_path):
        with open(json_path, "r") as f:
            data = json.load(f)
        return data
    else:
        return None


def load_pairs(pair_path):
    if os.path.exists(pair_path):
        with open(pair_path, "r") as f:
            lines = f.readlines()
        pairs = [ast.literal_eval(x) for x in lines]
        return pairs
    else:
        return None


def random_seed(seed):
    if seed is None:
        seed = 0
    np.random.seed(seed)
    torch.manual_seed(seed)
    random.seed(seed)


def dedup(recursive_list):
    final_set = set()
    for x in recursive_list:
        if (
            isinstance(x, list)
            or isinstance(x, omegaconf.ListConfig)
            or isinstance(x, omegaconf.listconfig.ListConfig)
        ):
            curr_result = dedup(x)
            for y in curr_result:
                final_set.add(y)
        else:
            final_set.add(x)
    return list(final_set)
