import os, sys
import pandas as pd
import numpy as np
import chimerax
from chimerax.core.commands import run

# Define consts
try:
    uniprot_id = 'P00720'
    structures = ['3lzm', '2lzm', '1lyd', '4lzm', '6lzm', '7lzm', '5lzm']
except (IOError, FileNotFoundError) as err:
    print("Cant open file:", str(err));
    #sys.exit(1)

cwd = os.getcwd()
DATA_PATH_ROOT = cwd + '/master_thesis/'

shape = (len(structures), len(structures))
pw_matrix = np.zeros(shape)
for pdb in structures:
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
np.save(DATA_PATH_ROOT + '/data/arodz12/pw_dist_rmsd/M_P00720_wt', pw_matrix)   

# Close session
run(session, "close")

