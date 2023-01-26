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
    uniprot_id = 'P0A7Y4'
    pdb_id_ref = DATA_PATH_ROOT + "data/results_af2/P0A7Y4/P0A7Y4_2rn2.result/P0A7Y4_2rn2_98129_unrelaxed_rank_1_model_2.pdb"
    af2_struc_list = ['1gob', '1kva', '1kvb', '1lav', '1law', '1rbr', '1rbs', '1rbt',
       '1rbu', '1rbv', '1rda', '1rdb']
except (IOError, FileNotFoundError) as err:
    print("Cant open file:", str(err));
    #sys.exit(1)

# Get mutation type
def get_pdb_seq(pdb, url):
    
    criteria = ['pdb_sequence']
    wt_seq = 'MLKQVEIFTDGSCLGNPGPGGYGAILRYRGREKTFSAGYTRTTNNRMELMAAIVALEALKEHCEVILSTDSQYVRQGITQWIHNWKKRGWKTADKKPVKNVDLWQRLDAALGQHQIKWEWVKGHAGHPENERCDELARAAAMNPTLEDTGYQVEV'
    df = pd.read_json(url)
    mutation_list = list()
    
    # Get mutation sites 
    try:
        df = pd.read_json(url)
        seq_struc = df[pdb][0][criteria[0]]
            
        # Find mismatches
        point_mutations = [wt_seq[mut]+str(mut+1)+seq_struc[mut] for mut in range(len(wt_seq)) if wt_seq[mut] != seq_struc[mut]]
    except TypeError:
        print("No PDBs found")

    return point_mutations


for af2_struc in range(len(af2_struc_list)):
    # Define data path & change dir
    DATA_PATH = DATA_PATH_ROOT + 'data/results_af2/' + uniprot_id + '/' + uniprot_id + '_' + af2_struc_list[af2_struc] + '.result'
    os.chdir(DATA_PATH) 
    
    # Load data
    dir_list = os.listdir(DATA_PATH)
    dir_list = sorted(dir_list)
    
    # Get mutation site
    url = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/%s/' % af2_struc_list[af2_struc]
    mutation = get_pdb_seq(af2_struc_list[af2_struc],url)
    
    # Run ChimeraX part of script
    run(session, "open {}".format(pdb_id_ref))
    models = [run(session, "open {}".format(pdb)) for pdb in dir_list if '.pdb' in pdb]
    model_names = [pdb for pdb in dir_list if '.pdb' in pdb]
    
    # Run matchmaker to obtain RMSD for each mutant to the wild type
    rmsd_list = list()
    for i in range(2,len(models)+2):
        try:
            rmsd = run(session, "mm #1/A to #{}/A".format(int(i)))
            rmsd_list.append([model_names[i-2],rmsd[0]['full RMSD'], mutation[:]])
        except (chimerax.core.errors.UserError) as err:
            with open('error.log', 'a') as f:
                print("Can't calculate RMSD for {}:".format(i, models[i]), str(err), file=f);
                f.close()
    # Save file
    df = pd.DataFrame(rmsd_list, columns=["name", "rmsd", "mut"])
    df.to_csv(DATA_PATH_ROOT + 'data/results_rmsd/results_rmsd_af2/results_rmsd_P0A7Y4_mut/rmsd_mut_{}_{}.csv'.format(uniprot_id, af2_struc_list[af2_struc]), index=False)

    # Close session
   # run(session, "close")

