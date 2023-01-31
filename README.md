# Exploring AlphaFold2’s capability to distinguish structural different variants’ conformation
This repository contains the code and data used in the master thesis.

### Abstract
AlphaFold2 (AF2) is a deep learning model developed by Google DeepMind, that accurately predicts the three-dimensional structure of proteins with atomic precision. This master’s thesis evaluates AF2’s capability to predict variant protein structures. AF2 was compared to a dataset of 116 variant crystal structures and four wild-type (WT) crystal structures from four different proteins. Furthermore, for each of the four WT structures, similar WT crystal structures were obtained through the Protein Data Bank (PDB). This was done to estimate the experimental inaccuracies in protein structure determination and was used as a control dataset. For every crystal variant structure, an AF2 structure was generated. The Root Mean Square Deviation (RMSD) was calculated for each structure pair and for each of the four protein variant datasets for both the crystal structures and the AF2 generated structures. The RMSD was calculated for each structure pair and for each of the four protein WT control datasets as well. The results showed that AF2 predicts protein variant structures of high quality but was not able to clearly and reliably distinguish between variant and WT structures as the calculated differences were within the range of experimental inaccuracy. The AF2 predicted χ1 angles were evaluated as well showing again that the predicted χ1 angles resemble both variant and WT χ1 angles. The results from AF2 were assessed further by comparing them to ESMFold, another protein structure prediction model. Furthermore, ablation of the input Multiple Sequence Alignment (MSA) to AF2 was shown to have little effect on producing structures that more closely resemble crystal variant structures compared to using the entire MSA.
### Dependencies
- Clone the repository.
- Setup the environment with the .yaml-file that has all the dependencies.
- Then open jupyter lab and run the different notebooks.

```bash
    conda env create -f  master_thesis.yml
    conda activate master_thesis
    jupyter-lab
```

The notebook folder contains all the code for the data analysis.

The scripts folder contains all the scripts to calculate RMSDs in ChimeraX. For this to work ChimeraX must be installed locally (https://www.cgl.ucsf.edu/chimerax/download.html)
To reproduce a specific calculation simply type:
```bash
    open example.py
```

The scripts_binf folder contains code from the groups server system and was used to create the WT control datasets. It is not executable and provided for documentation only.
