# Exploring AlphaFold2’s capability to distinguish structural different variants’ conformation
This repository contains the code and data used in the master thesis.

### Abstract
AlphaFold2 (AF2) is a deep learning model developed by Google DeepMind, that accurately predicts the three-dimensional structure of proteins with atomic precision. This master’s thesis evaluates AF2’s capability to predict variant protein structures. AF2 was compared to a dataset of 116 variant crystal structures and four wild-type (WT) crystal structures from four different proteins. For each of the 116 variant crystal structures and for each of the four WT structures the corresponding AF2 structures were generated. Furthermore, a WT crystal structure control dataset was obtained and used to assess the inherent inaccuracies in protein structure determination. The Root Mean Square Deviation (RMSD) was calculated for each structure pair and for each of the four protein variant datasets. The RMSD was calculated for each structure pair and for each of the four protein WT datasets as well. The results showed that AF2 predicts protein variant structures of high quality but was not able to clearly and reliably distinguish between variant and WT structures as the calculated differences were within the range of experimental inaccuracy. The AF2 predicted χ1 angles were evaluated as well but it was again shown that the predicted χ1 angles resembles both variant and WT χ1 angles. The results from AF2 were assessed further by comparing to ESMFold, another protein structure prediction model. Furthermore, ablation of the input Multiple Sequence Align- ment (MSA) to AF2 was shown to have little effect on producing structures that more closely resemble crystal variant structures compared to using the entire MSA.
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
