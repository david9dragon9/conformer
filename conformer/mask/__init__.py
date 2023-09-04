from conformer.mask.random_selected import RandomSelectedMSAMasker


def get_masker(config):
    if config.mask.name == "random_selected":
        return RandomSelectedMSAMasker(config.mask)
    else:
        raise ValueError(f"Unsupported masker: {config.name}")
