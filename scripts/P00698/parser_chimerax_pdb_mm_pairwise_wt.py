import os, sys
import pandas as pd
import numpy as np
import chimerax
from chimerax.core.commands import run

# Define consts
try:
    uniprot_id = 'P00698'
    structures = ['2vb1', '3lzt', '1iee', '4lzt', '7bcu', '7br5', '5lxw', '4b4e', '7avf', '7der', '7deq', '3ajn', '7bd0', '7bcx', '7bdm', '6agn', '6adf', '4ng8', '5k2q', '1lks', '4b0d', '6f1r', '4ngj', '5myy', '7lpm', '7r1q', '6ybi', '6ybf', '7llp', '7lp6', '1v7t', '8aj3', '6a10', '7lpl', '1v7s', '6g8a', '5f14', '5iel', '2z19', '7ol8', '2d4k', '4b49', '2z18', '2d4j', '2d4i', '3ato', '6abn', '6k5q', '2zyp', '7q0u', '5kxm', '5kxn', '5kxl', '4h9a', '4h8z', '4h8y', '4h8x', '4h9e', '4h90', '5kxt', '5v4g', '5kxo', '5ky1', '7jmu', '5kxw', '5kxs', '7vgo', '5f16', '5kxp', '4htk', '5ebh', '4h94', '1qio', '4n8z', '5kxy', '6fsj', '5kj9', '5fst', '5kxz', '6pl9', '6yjx', '5kxx', '4h9h', '4h91', '4h9i', '4h9f', '4h93', '3d9a', '2hub', '4bap', '4lyb', '6o2h', '5v4h', '5dla', '6hy8', '6agr', '8aj4', '7pnh', '4rds', '6bo1']
except (IOError, FileNotFoundError) as err:
    print("Cant open file:", str(err));
    #sys.exit(1)

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
np.save( '/Users/holger/Desktop/master_thesis/data/arodz12/pw_dist_rmsd/M_P00698_wt', pw_matrix)   

# Close session
run(session, "close")

