import os, sys
import pandas as pd
import numpy as np
import chimerax
from chimerax.core.commands import run

# Define consts
try:
    uniprot_id = 'P0A7Y4'
    pdb_id_ref = '2rn2'
    af2_struc_list = ['1gob', '1kva', '1kvb', '1lav', '1law', '1rbr', '1rbs', '1rbt',
       '1rbu', '1rbv', '1rda', '1rdb']
except (IOError, FileNotFoundError) as err:
    print("Cant open file:", str(err));
    #sys.exit(1)

for count, pdb in enumerate(af2_struc_list):
    # Define data path & change dir
    DATA_PATH = '/Users/holger/Desktop/master_thesis/data/results_pdb/output_' + uniprot_id + '/pdbs/' + pdb + '.pdb'
    
    # Run ChimeraX part of script
    if count == 0:
        run(session, "open {}".format(pdb_id_ref))
    run(session, "open {}".format(pdb))

    # Run matchmaker to obtain RMSD for each mutant to the wild type
    rmsd_list = list()
    try:
        rmsd = run(session, "mm #1/A to #{}/A".format(int(count+2)))
        rmsd_list.append([pdb,rmsd[0]['full RMSD']])
        rmsd = run(session, "mm #1/B to #{}/B".format(int(count+2)))
        rmsd_list.append([pdb,rmsd[0]['full RMSD']])
    except (chimerax.core.errors.UserError) as err:
        with open('/Users/holger/Desktop/master_thesis/data/results_rmsd/results_rmsd_pdb/' + 'error.log', 'w') as f:
            print("Can't calculate RMSD for {}:".format(pdb, str(err), file=f));
            f.close()
    # Save file
    df = pd.DataFrame(rmsd_list, columns=["name","rmsd"])
    df.to_csv('/Users/holger/Desktop/master_thesis/data/results_rmsd/results_rmsd_pdb/results_rmsd_{}_mut/rmsd_mut_a_b_{}_{}.csv'.format(uniprot_id,uniprot_id, pdb), index=False)

# Close session
#run(session, "close")

