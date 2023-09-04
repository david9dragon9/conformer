import hydra
from conformer.eval.eval_inputs import EvalInputs
from conformer.eval.download_msa import download_msa
from conformer.eval.download_structures import download_structures
from conformer.eval.eval_method import run_method
from conformer.eval.eval_baseline import eval_baseline
from conformer.eval.evaluate_rmsd import evaluate_and_filter_pairs
from conformer.eval.final_evaluation import final_evaluation


@hydra.main("./conformer/configs/")
def main(cfg):
    eval_inputs = EvalInputs(
        config=cfg,
        pair_path=cfg.pair_path,
        output_dir=cfg.output_dir,
        sequence_dict_path=cfg.sequence_dict_path,
        alignment_dir=cfg.alignment_dir,
        template_dir=cfg.template_dir,
        structure_dir=cfg.structure_dir,
        baseline_outputs_dir=cfg.baseline_outputs_dir,
        method_outputs_dir=cfg.method_outputs_dir,
        leans_path=cfg.leans_path,
        filtered_pair_path=cfg.filtered_pair_path,
        method_output_paths=cfg.method_output_paths,
    )
    if cfg.download_msa:
        eval_inputs = download_msa(eval_inputs)
    if cfg.download_structures:
        eval_inputs = download_structures(eval_inputs)
    if cfg.eval_baseline:
        eval_inputs = eval_baseline(eval_inputs)
    if cfg.evaluate_rmsd:
        eval_inputs = evaluate_and_filter_pairs(eval_inputs)
    if cfg.eval_method:
        eval_inputs = run_method(eval_inputs)
    if cfg.final_evaluation:
        eval_inputs = final_evaluation(eval_inputs)


if __name__ == "__main__":
    main()
