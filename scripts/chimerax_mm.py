import chimerax, os, sys
from chimerax.core.commands import run
import pandas as pd

# Define consts
pdb_list = ['1gob', '1kva', '1kvb', '1lav', '1law', '1rbr', '1rbs', '1rbt',
       '1rbu', '1rbv', '1rda', '1rdb']
# Load models
models = [run(session, "open {}".format(pdb)) for pdb in pdb_list]

# Run matchmaker to obtain RMSD for each mutant to the wild type
rmsd_list = list()
for i in range(len(pdb_list)):
    try:
        rmsd = run(session, "mm #1/A to #{}/A".format(int(i+1)))
        rmsd_list.append(rmsd[0]['full RMSD'])
        rmsd = run(session, "mm #1/B to #{}/B".format(int(i+1)))
        rmsd_list.append(rmsd[0]['full RMSD'])
    except (chimerax.core.errors.UserError) as err:
        print("Can't calculate RMSD for {}:".format(pdb_list[i]))

# Save file
df = pd.DataFrame(rmsd_list, columns=["rmsd"])
df.to_csv('/Users/holger/Desktop/master_thesis/data/results_rmsd/results_rmsd_pdb/results_rmsd_P0A7Y4_mut/rmsd_wt_P0A7Y4_a&b.csv', index=False)

# Close session
run(session, "close")

