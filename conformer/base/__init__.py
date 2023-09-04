def get_base_model(config):
    if config.base.name == "openfold":
        from conformer.base.openfold import OpenFold

        return OpenFold(config.base, full_config=config)
    elif config.base.name == "eigenfold":
        from conformer.base.eigenfold import EigenFold

        return EigenFold(config.base, full_config=config)
    elif config.base.name == "omegafold":
        from conformer.base.omegafold import OmegaFold

        return OmegaFold(config.base, full_config=config)
    else:
        raise ValueError(f"Unsupported base model: {config.base}")
