import os, sys
import pandas as pd
import numpy as np
import chimerax
from chimerax.core.commands import run
import urllib.request

# Define consts
try:
    uniprot_id = 'P0A7Y4'
    af2_struc_list = ['1kvb' for i in range(10)]
    pdb_id_ref = ['1kvb' for i in range(len(af2_struc_list))]
except (IOError, FileNotFoundError) as err:
    print("Cant open file:", str(err));
    #sys.exit(1)

cwd = os.getcwd()
DATA_PATH_ROOT = cwd + '/master_thesis/'

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
    DATA_PATH = DATA_PATH_ROOT + 'data/results_af2/' + uniprot_id + '/' + uniprot_id + '_1kvb_32/' + uniprot_id + '_1kvb_32_' + str(af2_struc+1) + '.result'
    os.chdir(DATA_PATH) 
    
    # Load data
    dir_list = os.listdir(DATA_PATH)
    dir_list = sorted(dir_list)
    
    # Get mutation site
    url = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/%s/' % af2_struc_list[af2_struc]
    mutation = get_pdb_seq(af2_struc_list[af2_struc],url)
    
    # Run ChimeraX part of script
    run(session, "open {}".format(pdb_id_ref[af2_struc]))
    models = [run(session, "open {}".format(pdb)) for pdb in dir_list if '.pdb' in pdb]
    model_names = [pdb for pdb in dir_list if '.pdb' in pdb]
    
    # Run matchmaker to obtain RMSD for each mutant to mutant
    rmsd_list = list()
    if af2_struc == 0:
        offset = 0
    elif af2_struc == 1:
        offset = 5
    elif af2_struc >= 2:
        offset += 5
    for i in range(2,len(models)+2):
        try:
            rmsd = run(session, "mm #{}/A to #{}/A".format(int(af2_struc+1+offset),int(af2_struc+offset+i)))
            rmsd_list.append([model_names[i-2],rmsd[0]['full RMSD'], mutation[:]])
        except (chimerax.core.errors.UserError) as err:
            with open('error.log', 'a') as f:
                print("Can't calculate RMSD for {}:".format(i, models[i]), str(err), file=f);
                f.close()
    
    # Save file
    df = pd.DataFrame(rmsd_list, columns=["name", "rmsd", "mut"])
    df.to_csv(DATA_PATH_ROOT + 'data/results_rmsd/results_rmsd_af2/results_rmsd_1kvb_32_mut/rmsd_1kvb2mut_{}_{}.csv'.format(uniprot_id, str(af2_struc+1)), index=False)

    # Close session
    #run(session, "close")

