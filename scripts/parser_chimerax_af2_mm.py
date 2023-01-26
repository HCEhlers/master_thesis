import os, sys
import pandas as pd
import numpy as np
import chimerax
from chimerax.core.commands import run

# Define consts
cwd = os.getcwd()
DATA_PATH_ROOT = cwd + '/master_thesis/'

try:
    uniprot_id = 'P61626'
    pdb_id_ref = DATA_PATH_ROOT + 'data/results_af2/P61626/P61626_default/P61626_1lz1.result/P61626_1lz1_ca876_unrelaxed_rank_1_model_5.pdb'
    af2_struc_list =  ['2nwd',
 '5lsh',
 '1jsf',
 '7xf6',
 '1iwt',
 '1iwu',
 '1iwv',
 '1iww',
 '1iwx',
 '1iwy',
 '1jwr',
 '6lfh',
 '1iwz',
 '1lzr',
 '3fe0',
 '1rex',
 '1lz1',
 '7xf7',
 '1lzs',
 '7xf8',
 '1rey',
 '1rez',
 '2zil',
 '4ml7',
 '2zik',
 '2zij',
 '4i0c',
 '1rem',
 '1re2',
 '5lvk']


except (IOError, FileNotFoundError) as err:
    print("Cant open file:", str(err));
    #sys.exit(1)

for af2_struc in range(len(af2_struc_list)):
    # Define data path & change dir
    DATA_PATH = DATA_PATH_ROOT + 'data/results_af2/' + uniprot_id + '/' + uniprot_id + '_default/' + uniprot_id + '_' + af2_struc_list[af2_struc] + '.result'
    os.chdir(DATA_PATH) 
    
    # Load data
    dir_list = os.listdir(DATA_PATH)
    dir_list = sorted(dir_list)
    
    # Run ChimeraX part of script
    run(session, "open {}".format(pdb_id_ref))
    models = [run(session, "open {}".format(pdb)) for pdb in dir_list if '.pdb' in pdb]
    model_names = [pdb for pdb in dir_list if '.pdb' in pdb]
    
    # Run matchmaker to obtain RMSD for each mutant to the wild type
    rmsd_list = list()
    for i in range(2,len(models)+2):
        try:
            rmsd = run(session, "mm #1/A to #{}/A".format(int(i)))
            rmsd_list.append([model_names[i-2],rmsd[0]['full RMSD']])
        except (chimerax.core.errors.UserError) as err:
            with open('error.log', 'a') as f:
                print("Can't calculate RMSD for {}:".format(i, models[i]), str(err), file=f);
                f.close()
    # Save file
    df = pd.DataFrame(rmsd_list, columns=["name", "rmsd"])
    df.to_csv(DATA_PATH_ROOT + 'data/results_rmsd/results_rmsd_af2/results_rmsd_P61626_wt/rmsd_wt_{}_{}.csv'.format(uniprot_id, af2_struc_list[af2_struc]), index=False)

    # Close session
    run(session, "close")

