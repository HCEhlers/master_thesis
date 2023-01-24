import os, sys
import pandas as pd
import numpy as np
import chimerax
from chimerax.core.commands import run
import urllib.request

# Define consts
try:
    uniprot_id = 'P00698'
    esmfold_path = "/Users/holger/Desktop/master_thesis/data/results_esmfold/structures/P00698_rec1/"
    pdb_id_ref_list = ['1hem', '1heo', '1her']
    rmsd_all_list = list()
except (IOError, FileNotFoundError) as err:
    print("Cant open file:", str(err));
    #sys.exit(1)

all_rmsds_list = list()
for pdb in range(len(pdb_id_ref_list)):
    # Run ChimeraX part of script
    run(session, "open {}".format(esmfold_path + pdb_id_ref_list[pdb] + '.pdb'))
    run(session, "open 4lyz")

for model in range(0,len(pdb_id_ref_list)*2,2):
    rmsd_list = list()
    rmsd = run(session, "mm #{}/A to #{}/A".format(int(model+1),int(model+2)))
    rmsd_list.append(rmsd[0]['full RMSD'])
    rmsd_all_list.append(rmsd_list)

for rmsd in range(len(rmsd_all_list)):
    filepath = '/Users/holger/Desktop/master_thesis/data/results_rmsd/results_rmsd_esmfold/P00698/P00698_rec1_WT/rmsd_esmfold_to_mut_{}_{}.csv'.format(uniprot_id, pdb_id_ref_list[rmsd])
    df = pd.DataFrame(rmsd_all_list[rmsd], columns=["rmsd"])
    df.to_csv(filepath, index=False) 
    #russion, "close")
