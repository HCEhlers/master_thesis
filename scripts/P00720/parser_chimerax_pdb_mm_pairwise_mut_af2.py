import os, sys
import pandas as pd
import numpy as np
import chimerax
from chimerax.core.commands import run

# Define consts
try:
    uniprot_id = 'P00720'
    structures = ['1dya', '1dyb', '1dyc', '1dyd', '1dye', '1dyf', '1dyg', '1l00',
       '1l02', '1l03', '1l05', '1l06', '1l07', '1l08', '1l09', '1l10',
       '1l11', '1l12', '1l13', '1l14', '1l15', '1l16', '1l19', '1l20',
       '1l21', '1l22', '1l23', '1l24', '1l25', '1l26', '1l27', '1l28',
       '1l29', '1l30', '1l31', '1l32', '1l33', '1l34', '1l37', '1l38',
       '1l44', '1l45', '1l46', '1l47', '1l48', '1l52', '1l53', '1l56',
       '1l57', '1l58', '1l60', '1l69', '1l98', '1l99']

except (IOError, FileNotFoundError) as err:
    print("Cant open file:", str(err));
    #sys.exit(1)

shape = (len(structures), len(structures))
pw_matrix = np.zeros(shape)
muts_rmsd = []
pdb1 = []
pdb2 = []
for pdb in structures:
    # Define data path
    DATA_PATH = '/Users/holger/Desktop/master_thesis/data/results_af2/' + uniprot_id + '/' + uniprot_id + '_default/'  + uniprot_id + '_' + pdb + '.result' 

    # Load data
    dir_list = os.listdir(DATA_PATH)
    dir_list = sorted(dir_list)

    # Get rank1 model
    model = [pdb for pdb in dir_list if '.pdb' in pdb and 'rank_1' in pdb]
    print(model)
    # Run ChimeraX part of script
    run(session, "open {}/{}".format(DATA_PATH, model[0]))

for i in range(pw_matrix.shape[0]):
    for j in range(pw_matrix.shape[1]):
        if i == j:
            pw_matrix[i,j] = 0
        else:
            try:
                rmsd = run(session, "mm #{}/A to #{}/A".format(int(i+1), int(j+1)))
                rmsd = rmsd[0]['full RMSD']
                pw_matrix[i,j] = rmsd
                if rmsd > 0.126 + 2*0.061:
                    muts_rmsd.append(rmsd)
                    pdb1.append(structures[int(i)])
                    pdb2.append(structures[int(j)])
            except (chimerax.core.errors.UserError) as err:
                print("Can't calculate RMSD for {} and {}".format(structures[i], structures[j]));
np.save( '/Users/holger/Desktop/master_thesis/data/arodz12/pw_dist_rmsd/M_pair_P00720_mut_af2', pw_matrix)   

# Save most significant variants to dataframe
index = [i for i in range(len(muts_rmsd))]
columns = ['rmsd', 'pdb1', 'pdb2']
df = pd.DataFrame({columns[0]: muts_rmsd, columns[1]: pdb1, columns[2]: pdb2})
print(df)
df.to_csv('/Users/holger/Desktop/master_thesis/data/arodz12/pw_dist_rmsd/signi_rmsd_paris.csv')
# Close session
#run(session, "close")
