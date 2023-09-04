from conformer.msa.af_cluster import AFClusterMSASubsampler
from conformer.msa.stochastic import StochasticMSASubsampler


def get_msa_subsampler(config):
    if config.msa.name == "af_cluster":
        return AFClusterMSASubsampler(config.msa)
    elif config.msa.name == "stochastic":
        return StochasticMSASubsampler(config.msa)
    else:
        raise ValueError(f"Unsupported subsampler: {config.msa.name}")
