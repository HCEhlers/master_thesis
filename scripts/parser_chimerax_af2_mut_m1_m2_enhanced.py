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
    uniprot_id = 'P00698'
    af2_struc_list = ['1hem', '1heo', '1her']
except (IOError, FileNotFoundError) as err:
    print("Cant open file:", str(err));
    #sys.exit(1)

rmsd_list = []
for af2_struc in range(len(af2_struc_list)):
    # Define data path & change dir
    DATA_PATH = DATA_PATH_ROOT + 'data/results_af2/' + uniprot_id + '/' + uniprot_id + '_' + af2_struc_list[af2_struc] + '.result'
    os.chdir(DATA_PATH) 
    
    # Load data
    dir_list = os.listdir(DATA_PATH)
    dir_list = sorted(dir_list)
    
    # Run ChimeraX part of script
    models = [run(session, "open {}".format(pdb)) for pdb in dir_list if ('.pdb' in pdb and 'model_5' in pdb or 'model_1' in pdb and '.json' not in pdb)]

    
    # Run matchmaker to obtain RMSD for each mutant to mutant
    try:
        rmsd = run(session, "mm #{}/A to #{}/A".format(int(af2_struc+1),int(af2_struc+2)))[0]['final RMSD']
        rmsd_list.append(rmsd)
    except (chimerax.core.errors.UserError) as err:
        print('ChimeraX error')

# Save file
df = pd.DataFrame(rmsd_list, columns=["rmsd"])
df.to_csv(DATA_PATH_ROOT + 'data/results_rmsd/results_rmsd_af2/results_rmsd_P00698_m1_to_m5.csv', index=False)

# Close session
#run(session, "close")

