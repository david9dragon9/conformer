# Conformer: A Library for Multi-Conformational Protein Folding
Machine-learning based protein structure prediction methods (e.g. AlphaFold, RoseTTAFold, etc...) have revolutionized the field of computational biology with highly accurate protein structure predictions without the cost of experiments.

However, many of these protein folding methods are unable to predict multiple conformational states of proteins, instead only predicting one of the two conformational states. Many important proteins (e.g. GPCRs, transport proteins, protein pumps and channels) have more than one conformational state, making this problem critical.

Conformer provides the tools necessary to for predicting multiple conformational states of given proteins, with implementations of multiple existing methods for multi-conformational protein folding, and seamless composability between different base models and multi-conformational protein folding specialized techniques.

# Capabilities
Conformer provides several different base models and techniques for multi-conformational protein folding that can be easily composed with each other.

**More methods will be coming soon!!!**

## Base Models
- OpenFold (AlphaFold weights can be used with the same code): OpenFold is a PyTorch reproduction of the original AlphaFold2 result.
- OmegaFold: OmegaFold is a highly accurate protein structure prediction method that uses a pre-trained protein language model and no MSA sequences.
- EigenFold: EigenFold is a conditional diffusion-based protein structure prediction method that predicts protein structure from only sequence and uses OmegaFold embeddings, gradually denoising the predicted structure using diffusion.

## Multi-Conformational Protein Folding Techniques
- AF-Cluster (MSA clustering): AF-Cluster attempts to predict multiple conformations of a given protein by clustering the MSA sequences into different clusters that have biases towards different conformational states.
- MSA Stochastic Subsampling: MSA Stochastic Subsampling randomly selects shallow subsets of the full MSA that randomly bias the model towards particular conformational states.
- MSA Masking: MSA Masking works by masking out particular determining factors that bias the model towards particular conformational states in order to coax the model into predicting alternative conformational states.
- Diffusion: Diffusion randomly samples protein structures by starting with random noise initially and randomly taking steps until a final structure is reached.

## Other
- MSA and structure download: Generate and download MSAs and structures for arbitrary protein sequences.
- Smart RMSD that uses alignment with gap penalty in order to find accurate RMSD to evaluate difference between protein structures.
- End-to-end evaluation code from PDB pairs to final accuracy.
- Hydra for easy configuration.
- wandb for result logging.

# Design
Conformer is centered around the `Folder` class: The `Folder` allows you to predict structure using any combination of base model and multi-conformational technique by passing in the inputs. The `Folder` class can be easily constructed using a config loaded from a YAML file. Base models, MSA selection, and MSA masking methods are housed in separate folders, and the EigenFold diffusion method is a separate base model.

# Setup
## General setup
For all methods, run:
```
pip install -r requirements.txt
chmod +x setup.sh && ./setup.sh
```

For each specific base model (e.g OpenFold, EigenFold, OmegaFold), run the corresponding setup file (e.g. `setup_eigenfold.sh` for EigenFold, etc...)

# Fold sequences
To fold a sequence, you can either use the `Folder` class interactively (e.g. in an IPython notebook), or by running `run_folder.py`.

You will need the following inputs:
- fasta_dir: A directory with .fasta files each containing a single ID and sequence. Example FASTA file:
```
>ID
ARNDCEQ
```
- alignment_dir: A directory with folders named **with the IDs in the FASTA files**. Each folder should have a single .a3m file which has the MSA. Example MSA file:
```
>information about MSA seq #1
ARNDCEQ
>information about MSA seq #2
ARND--Q
```
- template_dir: A directory with templates. **Conformer has not been tested with templates.**
- output_dir: The output directory. Modified MSAs will also be saved here.
- structure_dir: Optionally set the directory where PDB files will be downloaded and loaded from if you are using wandb to perform RMSD logging.

Set the following variables to use different methods:
- base: Set to either "openfold", "omegafold", or "eigenfold".
- msa: Set to either "af_cluster", "stochastic", or null to use AF-Cluster, Stochastic MSA Subsampling, or nothing respectively.
- mask: Set to either "random_selected" or null to use MSA masking or nothing.

Parameters that **must be set** for each method:
- AF-Cluster: None. You can use the defaults or tweak the parameters.
- MSA Subsampling: None. You can use the defaults or tweak the parameters.
- MSA Masking:
- - chosen_residues: which residues to ignore or preserve when masking.
- - mask_mode: whether to ignore the chosen residues (mask them out with probability p) or preserve the chosen residues (keep these, and mask everything else with probability p).

Parameters that **must be set** for each base model:
- OpenFold:
- - config_preset: Set to an OpenFold config preset (e.g. "model_3", "model_1", etc...)
- - openfold_checkpoint_path: Set to the OpenFold checkpoint path. To download the OpenFold weights, please refer to https://github.com/aqlaboratory/openfold/blob/main/notebooks/OpenFold.ipynb
- OmegaFold:
- - device: Set to "cuda" or "cpu" depending on what device you are running on
- - weights_file: Set to the weights file for OmegaFold.
- EigenFold:
- - device: Set to "cuda" or "cpu" depending on what device you are running on
- - model_dir: The directory in which the pretrained model is stored
- - ckpt: The path to the model checkpoint inside model_dir
- - lm_weights_path: The path to the OmegaFold weights.

For OpenFold and EigenFold, you can use `use_wandb` and `wandb` respectively to use wandb for result logging.

To run `run_folder.py`:
```
python3 run_folder.py --config-name basic_fold \
                      fasta_dir=/fasta/dir/ \
                      template_dir=/template/dir \
                      alignment_dir=/alignment/dir \
                      output_dir=/output/dir \
                      base.config_preset=model_3 \
                      base.openfold_checkpoint_path=/path/to/openfold.pt
```

To use `Folder` interactively:
```
from conformer.folder import Folder

# Load the config
import hydra
hydra.core.global_hydra.GlobalHydra.instance().clear()
hydra.initialize_config_dir("/content/conformer/conformer/configs/")
config = hydra.compose("basic_fold", overrides=["base=openfold"]) # Put other command line overrides here

# Create Folder
folder = Folder(config)

# Run Folder
folder.run(fasta_dir="/path/to/fasta/dir",
           alignment_dir="/path/to/alignment/dir",
           template_dir="/path/to/template/dir",
           output_dir="/path/to/output/dir",
           input_dict={"A": "ARNDCEQ"}) # You can pass an input dict to override fasta_dir instead as well
```

# End-to-end evaluation
To run an end-to-end evaluation:
1. Create a working folder for your evaluation. All intermediate outputs and final outputs will go here.
2. Create a file called pairs.txt inside your folder and put the pairs that you would like to evaluate on into this folder. Conformer will use the first ID in each pair as the sequence. Example pairs.txt file:
```
["1ABC_A", "1BCD_A"]
["2ABC_B", "2DCB_G"]
```
3. Create a separate config file in the conformer/configs folder with your method (e.g. with AF-Cluster, MSA masking). **You must do this if your method differs from the baseline. Save the config file for future reference.** Then, change eval.yaml to use your method config file in the method field (e.g. `my_config@method`)
Example method config file:
```
defaults:
 - base: openfold
 - msa: af_cluster
 - mask: null

fasta_dir: null
template_dir: null
alignment_dir: null 
output_dir: null
structure_dir: structures/
```


Then, run `run_eval.py`:

Important parameters to keep in mind:
- eval_config: Chosen out of ["best", "average"]. Indicates whether to choose the best output out of the sampled predicted structures or to calculate the percentage of predictions that successfully predicted the target state
- pair_path: This is where the input pairs are placed, typically set to {output_dir}/pairs.txt, default=output_dir/pairs.txt
- output_dir: This is the output directory where all intermediate and final outputs are stored, default=null
- sequence_dict_path: This is a path to a JSON with all PDB sequences (can be downloaded from https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz).
- alignment_dir: Where MSAs will be stored (or loaded from, if you already have precomputed MSAs), default=output_dir/alignments
- template_dir: Where templates will be stored (or loaded from, if you have precomputed templates). Optional if you are not using templates, default=output_dir/templates
- structure_dir: Where downloaded structures will be stored for the evaluation, default=output_dir/structures
- baseline_outputs_dir: Where baseline outputs will be stored, default=output_dir/baseline_outputs
- method_outputs_dir: Where method outputs will be stored, default=output_dir/method_outputs
- leans_path: Where the calculated leans JSON will be stored, default=output_dir/leans.json
- filtered_pair_path: Where the calculated filtered pairs will be stored, default=output_dir/filtered_pairs.txt
- method_output_paths: Where the calculated method output paths will be stored, default=output_dir/method_output_paths.txt

Stages of the evaluation:
- download_msa: Download MSAs from the ColabFold MSA API for each PDB ID in the pairs list (both PDB IDs in each pair).
- download_structures: Download structures from the PDB for each PDB ID in the pairs list (both PDB IDs in each pair). Used for evaluation and lean calculation.
- eval_baseline: Evaluate baseline for both PDB IDs in each pair to determine lean of the model.
- evaluate_rmsd: Calculate the lean of each predicted protein structure to both PDB IDs in each pair. Determines whether there is "bias" in the structure prediction.
- eval_method: Evaluates the method once for each PDB **pair**. As the outputs can be scattered within method_outputs_dir, this also keeps track of where each method prediction is stored, saves it to method_output_paths.txt, and prints it.
- final_evaluation: Evaluates whether the method was able to predict the alternative state of each pair (i.e. the conformation that the baseline was not able to predict.). Prints final accuracy.

While this evaluation pipeline is specifically aimed at identifying biases in the model and testing whether methods can predict the alternative states, you can run only certain stages (e.g. just evaluate your baseline) as well, simply by setting the boolean setting for each 

Example command:
```
python3 run_eval.py --config-name eval \
                     pair_path=/path/to/pairs/file.txt \
                     sequence_dict_path=/path/to/sequence/dict \
                     output_dir=/output/dir/ \
                     baseline.base.config_preset=model_3 \
                     baseline.base.openfold_checkpoint_path=/path/to/openfold/params \
                     method.base.config_preset=model_3 \
                     method.base.openfold_checkpoint_path=/path/to/openfold/params
```

# Roadmap
Coming soon:

Base Models:
- RoseTTAFold
- ESMFold

Multi-conformational techniques:
- Sequence Masking
- Template selection
- Residue selection methods for sequence and MSA masking

# FAQ
## What if my protein has more than 2 conformational states?
You can still run Folder and evaluate different methods on your protein, though the Conformer evaluation pipeline currently expects a list of pairs.

## Does Conformer support relaxation?
No, Conformer does not support relaxation at the moment.

## Does Conformer support templates?
Yes, Conformer does support templates, though it has not been extensively tested with templates.

# References
### AlphaFold
Code: https://github.com/deepmind/alphafold
Paper: Jumper, J., Evans, R., Pritzel, A. et al. Highly accurate protein structure prediction with AlphaFold. Nature 596, 583–589 (2021). https://doi.org/10.1038/s41586-021-03819-2

### ColabFold
Code: https://github.com/sokrypton/ColabFold
Paper: Mirdita M, Schütze K, Moriwaki Y, Heo L, Ovchinnikov S and Steinegger M. ColabFold: Making protein folding accessible to all.
Nature Methods (2022) doi: 10.1038/s41592-022-01488-1

### OpenFold
Code: https://github.com/aqlaboratory/openfold
Paper: Ahdritz, G., Bouatta, N., Kadyan, S., Xia, Q., Gerecke, W., O\textquoterightDonnell, T., Berenberg, D., Fisk, I., Zanichelli, N., Zhang, B., Nowaczynski, A., Wang, B., Stepniewska-Dziubinska, M., Zhang, S., Ojewole, A., Guney, M., Biderman, S., Watkins, A., Ra, S., Lorenzo, P., Nivon, L., Weitzner, B., Ban, Y.E., Sorger, P., Mostaque, E., Zhang, Z., Bonneau, R., & AlQuraishi, M. (2022). OpenFold: Retraining AlphaFold2 yields new insights into its learning mechanisms and capacity for generalization. bioRxiv.

### OmegaFold
Code: https://github.com/HeliXonProtein/omegafold
Paper: Wu, R., Ding, F., Wang, R., Shen, R., Zhang, X., Luo, S., Su, C., Wu, Z., Xie, Q., Berger, B., Ma, J., & Peng, J. (2022). High-resolution de novo structure prediction from primary sequence. bioRxiv.

### EigenFold
Code: https://github.com/bjing2016/EigenFold
Paper: Bowen Jing, Ezra Erives, Peter Pao-Huang, Gabriele Corso, Bonnie Berger, & Tommi Jaakkola. (2023). EigenFold: Generative Protein Structure Prediction with Diffusion Models.

### AF2 Conformations
Code: https://github.com/delalamo/af2_conformations
Paper: Alamo, D., Sala, D., Mchaourab, H., & Meiler, J. (2022). Sampling alternative conformational states of transporters and receptors with AlphaFold2. eLife, 11, e75751.

### AF-Cluster
Code: https://github.com/HWaymentSteele/AF_Cluster
Paper: Hannah K. Wayment-Steele, Sergey Ovchinnikov, Lucy Colwell, & Dorothee Kern (2022). Prediction of multiple conformational states by combining sequence clustering with AlphaFold2. bioRxiv.

### AF2 GPCR Kinase
Code: https://github.com/meilerlab/AF2_GPCR_Kinase
Paper: Sala, D., & Meiler, J. (2022). Biasing AlphaFold2 to predict GPCRs and Kinases with user-defined functional or structural properties. bioRxiv.

### SPEACH AF
Code: https://github.com/RSvan/SPEACH_AF
Paper: Stein, SPEACH_AF: Sampling protein ensembles and conformational heterogeneity with Alphafold2.
PLoS Comput Biol (2022) 18(8), e1010483 doi: 10.1371/journal.pcbi.1010483

### RMSD Computation
Code: https://github.com/Bernhard10/py_qcprot
Papers:
Liu P, Agrafiotis DK, & Theobald DL (2011) Reply to comment on: "Fast determination of the optimal rotation matrix for macromolecular superpositions." Journal of Computational Chemistry 32(1):185-186. [Open Access], doi:10.1002/jcc.21606

Liu P, Agrafiotis DK, & Theobald DL (2010) "Fast determination of the optimal rotation matrix for macromolecular superpositions." Journal of Computational Chemistry 31(7):1561-1563. [Open Access] doi:10.1002/jcc.21439

Douglas L Theobald (2005) "Rapid calculation of RMSDs using a quaternion-based characteristic polynomial." Acta Crystallogr A 61(4):478-480. [Open Access] doi:10.1107/S0108767305015266