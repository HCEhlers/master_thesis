import os, sys
import pandas as pd
import numpy as np
import chimerax
from chimerax.core.commands import run
import urllib.request

# Define consts
try:
    uniprot_id = 'P61626'
    pdb_id_ref = ['1gb5', '1gb6', '1gb7', '1gb8', '1gb9', '1gbo', '1gbw', '1gbx',
       '1gby', '1gbz', '1gfh', '1gfj', '1gfk', '1gfr', '1gft', '1gfu',
       '1gfv', '1inu', '1lhh', '1lhi', '1lhj', '1lhk', '1lhl', '1oub',
       '1ouc', '1oud', '1oue', '1ouh', '1oui', '1ouj', '1wqm', '1wqn',
       '1wqo', '1wqp', '1wqq', '1wqr', '1yam', '1yan', '1yao', '1yap',
       '1yaq', '2hea', '2heb', '2hec', '2hed', '2hee', '2hef']
    af2_struc_list = ['1gb5', '1gb6', '1gb7', '1gb8', '1gb9', '1gbo', '1gbw', '1gbx',
       '1gby', '1gbz', '1gfh', '1gfj', '1gfk', '1gfr', '1gft', '1gfu',
       '1gfv', '1inu', '1lhh', '1lhi', '1lhj', '1lhk', '1lhl', '1oub',
       '1ouc', '1oud', '1oue', '1ouh', '1oui', '1ouj', '1wqm', '1wqn',
       '1wqo', '1wqp', '1wqq', '1wqr', '1yam', '1yan', '1yao', '1yap',
       '1yaq', '2hea', '2heb', '2hec', '2hed', '2hee', '2hef']
except (IOError, FileNotFoundError) as err:
    print("Cant open file:", str(err));
    #sys.exit(1)

# Get mutation type
def get_pdb_seq(pdb, url):
    
    criteria = ['pdb_sequence']
    wt_seq = 'KVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGYNTRATNYNAGDRSTDYGIFQINSRYWCNDGKTPGAVNACHLSCSALLQDNIADAVACAKRVVRDPQGIRAWVAWRNRCQNRDVRQYVQGCGV'
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


for af2_struc in range(1,len(af2_struc_list)):
    # Define data path & change dir
    DATA_PATH = '/Users/holger/Desktop/master_thesis/data/results_af2/' + uniprot_id + '/' + uniprot_id + '_' + af2_struc_list[af2_struc] + '.result'
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
    if af2_struc == 1:
        offset = 0
    elif af2_struc == 2:
        offset = 5
    elif af2_struc > 2:
        offset += 5
    for i in range(2,len(models)+2):
        try:
            rmsd = run(session, "mm #{}/A to #{}/A".format(int(af2_struc+offset),int(af2_struc-1+offset+i)))
            rmsd_list.append([model_names[i-2],rmsd[0]['full RMSD'], mutation[:]])
        except (chimerax.core.errors.UserError) as err:
            with open('error.log', 'a') as f:
                print("Can't calculate RMSD for {}:".format(i, models[i]), str(err), file=f);
                f.close()
    # Save file
    df = pd.DataFrame(rmsd_list, columns=["name", "rmsd", "mut"])
    df.to_csv('/Users/holger/Desktop/master_thesis/data/results_rmsd/results_rmsd_af2/results_rmsd_P61626_mut2mut/rmsd_mut2mut_{}_{}.csv'.format(uniprot_id, af2_struc_list[af2_struc]), index=False)

    # Close session
    #run(session, "close")

