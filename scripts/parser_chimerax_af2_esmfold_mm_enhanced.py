import os, sys
import pandas as pd
import numpy as np
import chimerax
from chimerax.core.commands import run
import urllib.request

# Define consts
cwd = os.getcwd()
DATA_PATH_ROOT = cwd + '/master_thesis/'

try:
    uniprot_id = 'P61626'
    esmfold_path = DATA_PATH_ROOT + "/data/results_esmfold/structures/P61626_rec1/"
    pdb_id_ref_list = ['1gb5', '1gb6', '1gb7', '1gb8', '1gb9', '1gbo', '1gbw', '1gbx', '1gby', '1gbz', '1gfh', '1gfj', '1gfk', '1gfr', '1gft', '1gfu','1gfv', '1inu', '1lhh', '1lhi', '1lhj', '1lhk', '1lhl', '1oub','1ouc', '1oud', '1oue', '1ouh', '1oui', '1ouj', '1wqm', '1wqn','1wqo', '1wqp', '1wqq', '1wqr', '1yam', '1yan', '1yao', '1yap','1yaq', '2hea', '2heb', '2hec', '2hed', '2hee', '2hef']

    rmsd_all_list = list()
except (IOError, FileNotFoundError) as err:
    print("Cant open file:", str(err));
    #sys.exit(1)

all_rmsds_list = list()
for pdb in range(len(pdb_id_ref_list)):
    # Run ChimeraX part of script
    run(session, "open {}".format(esmfold_path + pdb_id_ref_list[pdb] + '.pdb'))
    run(session, "open {}".format(pdb_id_ref_list[pdb]))

for model in range(0,len(pdb_id_ref_list)*2,2):
    rmsd_list = list()
    rmsd = run(session, "mm #{}/A to #{}/A".format(int(model+1),int(model+2)))
    rmsd_list.append(rmsd[0]['full RMSD'])
    rmsd_all_list.append(rmsd_list)

for rmsd in range(len(rmsd_all_list)):
    filepath = DATA_PATH_ROOT + 'data/results_rmsd/results_rmsd_esmfold/P61626/P61626_rec1_MUT/rmsd_esmfold_to_mut_{}_{}.csv'.format(uniprot_id, pdb_id_ref_list[rmsd])
    df = pd.DataFrame(rmsd_all_list[rmsd], columns=["rmsd"])
    df.to_csv(filepath, index=False) 
    #russion, "close")
