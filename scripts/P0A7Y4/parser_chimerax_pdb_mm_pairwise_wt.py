import os, sys
import pandas as pd
import numpy as np
import chimerax
from chimerax.core.commands import run

# Define consts
try:
    uniprot_id = 'P0A7Y4'
    structures = ['2rn2', '7vsd', '7vsa', '7vsc', '7vsb', '7vse']
except (IOError, FileNotFoundError) as err:
    print("Cant open file:", str(err));
    #sys.exit(1)

shape = (len(structures), len(structures))
pw_matrix = np.zeros(shape)
for pdb in structures:
    # Define data path & change dir
    DATA_PATH = '/Users/holger/Desktop/master_thesis/data/results_pdb/output_' + uniprot_id + '/pdbs/' + pdb + '.pdb'
    
    # Run ChimeraX part of script
    run(session, "open {}".format(pdb))

for i in range(pw_matrix.shape[0]):
    for j in range(pw_matrix.shape[1]):
        if i == j:
            pw_matrix[i,j] = 0
        else:
            try:
                rmsd = run(session, "mm #{}/A to #{}/A".format(int(i+1), int(j+1)))
                rmsd = rmsd[0]['full RMSD']
                pw_matrix[i,j] = rmsd
            except (chimerax.core.errors.UserError) as err:
                print("Can't calculate RMSD for {} and {}".format(structures[i], structures[j]));
np.save( '/Users/holger/Desktop/master_thesis/data/arodz12/pw_dist_rmsd/M_{}_wt'.format(uniprot_id), pw_matrix)   

# Close session
run(session, "close")

