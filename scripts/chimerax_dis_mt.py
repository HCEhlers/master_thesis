import os, sys
import pandas as pd
import numpy as np
import chimerax
from chimerax.core.commands import run
import urllib.request

# Define consts
try:
    uniprot_id = 'P61626'
    af2_struc_list = ['1gb5', '1gb6', '1gb7', '1gb8', '1gb9', '1gbo', '1gbw', '1gbx',
       '1gby', '1gbz', '1gfh', '1gfj', '1gfk', '1gfr', '1gft', '1gfu',
       '1gfv', '1inu', '1lhh', '1lhi', '1lhj', '1lhk', '1lhl', '1oub',
       '1ouc', '1oud', '1oue', '1ouh', '1oui', '1ouj', '1wqm', '1wqn',
       '1wqo', '1wqp', '1wqq', '1wqr', '1yam', '1yan', '1yao', '1yap',
       '1yaq', '2hea', '2heb', '2hec', '2hed', '2hee', '2hef']
    rsn_N = ["['V74I']", "['V74L']", "['V74M']", "['V74F']", "['V110G']", "['V110I']", "['V110L']", "['V110M']", "['V110F']", "['V74Y']", "['V74D']", "['V74N']", "['V74R']", "['V110Y']", "['V110D']", "['V110N']", "['V110R']", "['V110P']", "['D91P']", "['Y124F']", "['Y20F']", "['Y38F']", "['Y45F']", "['Y54F']", "['Y63F']", "['I106V']", "['I23V']", "['I56V']", "['I59V']", "['I89V']"]
except (IOError, FileNotFoundError) as err:
    print("Cant open file:", str(err));
    #sys.exit(1)

# Get c beta distances
for af2_struc in range(len(af2_struc_list)):
    # Define data path & change dir
    DATA_PATH = '/Users/holger/Desktop/master_thesis/data/results_af2/' + uniprot_id + '/' + uniprot_id + '_' + af2_struc_list[af2_struc] + '.result'
    os.chdir(DATA_PATH) 
    
    # Load data
    dir_list = os.listdir(DATA_PATH)
    dir_list = sorted(dir_list)
    pdb_file = [file for file in  dir_list if 'rank_1' in file and '.pdb' in file] 

    # Run ChimeraX part of script
    run(session, "open {}".format(af2_struc_list[af2_struc]))
    run(session, "open {}".format(pdb_file[0]))
    
    dis_list = list()
    #if rsn_N[2] != 'G' and rsn_N[-3] != 'G' and rsn_N[2] != 'A' and rsn_N[-3] != 'A':
    #    c_beta_dis = run(session,"distance #1:{}@cb  #{}:{}@cb".format(int(rsn_N[af2_struc][3:-3]), int(af2_struc), int(rsn_N[af2_struc][3:-3])))
    #    dis_list.append(c_beta_dis)

    # Close session
    #run(session, "close")

