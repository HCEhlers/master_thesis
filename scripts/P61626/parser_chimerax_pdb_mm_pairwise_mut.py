import os, sys
import pandas as pd
import numpy as np
import chimerax
from chimerax.core.commands import run

# Define consts
try:
    uniprot_id = 'P61626'
    structures = ['1gb5', '1gb6', '1gb7', '1gb8', '1gb9', '1gbo', '1gbw', '1gbx', '1gby', '1gbz', '1gfh', '1gfj', '1gfk', '1gfr', '1gft', '1gfu','1gfv', '1inu', '1lhh', '1lhi', '1lhj', '1lhk', '1lhl', '1oub','1ouc', '1oud', '1oue', '1ouh', '1oui', '1ouj', '1wqm', '1wqn','1wqo', '1wqp', '1wqq', '1wqr', '1yam', '1yan', '1yao', '1yap','1yaq', '2hea', '2heb', '2hec', '2hed', '2hee', '2hef']

except (IOError, FileNotFoundError) as err:
    print("Cant open file:", str(err));
    #sys.exit(1)

shape = (len(structures), len(structures))
pw_matrix = np.zeros(shape)
muts_rmsd = []
pdb1 = []
pdb2 = []
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
                if rmsd > 2*0.171:
                    muts_rmsd.append(rmsd)
                    pdb1.append(structures[int(i)])
                    pdb2.append(structures[int(j)])
            except (chimerax.core.errors.UserError) as err:
                print("Can't calculate RMSD for {} and {}".format(structures[i], structures[j]));
np.save( '/Users/holger/Desktop/master_thesis/data/arodz12/pw_dist_rmsd/M_pair_mut', pw_matrix)   
# Save most significant variants to dataframe
index = [i for i in range(len(muts_rmsd))]
columns = ['rmsd', 'pdb1', 'pdb2']
df = pd.DataFrame({columns[0]: muts_rmsd, columns[1]: pdb1, columns[2]: pdb2})
print(df)
df.to_csv('/Users/holger/Desktop/master_thesis/data/arodz12/pw_dist_rmsd/signi_rmsd_pairs_pdb_P61626.csv')
# Close session
#run(session, "close")

