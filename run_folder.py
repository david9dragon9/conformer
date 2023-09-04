import hydra
from omegaconf import DictConfig, OmegaConf
from conformer.folder import Folder


@hydra.main("./conformer/configs/")
def main(cfg: DictConfig):
    folder = Folder(cfg)
    print(cfg)
    output_paths = folder.run(
        fasta_dir=cfg.fasta_dir,
        alignment_dir=cfg.alignment_dir,
        template_dir=cfg.template_dir,
        output_dir=cfg.output_dir,
        structure_dir=cfg.structure_dir,
    )
    print("--------------------------------------------------")
    print("--------------------------------------------------")
    print(output_paths)
    print("--------------------------------------------------")
    print("--------------------------------------------------")


if __name__ == "__main__":
    main()
