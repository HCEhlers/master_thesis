import os, sys
import pandas as pd
import numpy as np
import chimerax
from chimerax.core.commands import run
import urllib.request

# Define consts
try:
    uniprot_id = 'P61626'
    pdb_id_ref = '/Users/holger/Desktop/master_thesis/data/results_af2/P61626/P61626_1lz1.result/P61626_1lz1_ca876_unrelaxed_rank_1_model_5.pdb'
    wt_seq = 'KVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGYNTRATNYNAGDRSTDYGIFQINSRYWCNDGKTPGAVNACHLSCSALLQDNIADAVACAKRVVRDPQGIRAWVAWRNRCQNRDVRQYVQGCGV'
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
def get_pdb_seq(pdb, url, wt_seq):
    criteria = ['pdb_sequence']
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
   
    return point_mutations, seq_struc


for af2_struc in range(len(af2_struc_list)):
    # Define data path & change dir
    DATA_PATH = '/Users/holger/Desktop/master_thesis/data/results_af2/' + uniprot_id + '/' + uniprot_id + '_' + af2_struc_list[af2_struc] + '.result'
    os.chdir(DATA_PATH) 
    
    # Load data
    dir_list = os.listdir(DATA_PATH)
    dir_list = sorted(dir_list)
    
    # Get mutation site
    url = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/%s/' % af2_struc_list[af2_struc]
    mutation, seq_struc = get_pdb_seq(af2_struc_list[af2_struc],url, wt_seq)
    
    # Run ChimeraX part of script
    run(session, "open {}".format(pdb_id_ref))
    models = [run(session, "open {}".format(pdb)) for pdb in dir_list if '.pdb' in pdb]
    model_names = [pdb for pdb in dir_list if '.pdb' in pdb]
    
    # Run matchmaker to obtain RMSD for each mutant to the wild type
    rmsd_list = list()
    torsion_angle_list = list()
    for i in range(2,len(models)+2):
        try:
            # Superimpose structures and calculate RMSD
            rmsd = run(session, "mm #1/A to #{}/A".format(int(i)))
            rmsd_list.append([model_names[i-2],rmsd[0]['full RMSD'], mutation[0:]])
            
            # Calculate chi 1 angle for all positions for rank 1 structure only
            if i == 2:
                for res in range(len(seq_struc)):
                    # Arginine
                    if wt_seq[res] == 'R':
                        torison_angle = run(session, "torsion #{}/A:{}@n,ca,cb,cg".format(int(af2_struc), int(res)+1))
                        print(torison_angle)
                        torsion_angle_list.append(torison_angle)
                    # Asparagine
                    if wt_seq[res] == 'N':
                        torison_angle = run(session, "torsion #{}/A:{}@n,ca,cb,cg".format(int(af2_struc), int(res)+1))
                        torsion_angle_list.append(torison_angle)
                    # Aspartic acid
                    if wt_seq[res] == 'D':
                        torison_angle = run(session, "torsion #{}/A:{}@n,ca,cb,cg".format(int(af2_struc), int(res)+1))
                        torsion_angle_list.append(torison_angle)
                    # Cysteine
                    if wt_seq[res] == 'C':
                        torison_angle = run(session, "torsion #{}/A:{}@n,ca,cb,sg".format(int(af2_struc), int(res)+1))
                        torsion_angle_list.append(torison_angle)
                    # Glutamine
                    if wt_seq[res] == 'Q':
                        torison_angle = run(session, "torsion #{}/A:{}@n,ca,cb,cg".format(int(af2_struc), int(res)+1))
                        torsion_angle_list.append(torison_angle)
                    # Glutamic acid
                    if wt_seq[res] == 'E':
                        torison_angle = run(session, "torsion #{}/A:{}@n,ca,cb,cg".format(int(af2_struc), int(res)+1))
                        torsion_angle_list.append(torison_angle)
                    # Histidine
                    if wt_seq[res] == 'H':
                        torison_angle = run(session, "torsion #{}/A:{}@n,ca,cb,cg".format(int(af2_struc), int(res)+1))
                        torsion_angle_list.append(torison_angle)
                    # Isoleucine
                    if wt_seq[res] == 'I':
                        torison_angle = run(session, "torsion #{}/A:{}@n,ca,cb,cg1".format(int(af2_struc), int(res)+1))
                        torsion_angle_list.append(torison_angle)
                    # Leucine
                    if wt_seq[res] == 'L':
                        torison_angle = run(session, "torsion #{}/A:{}@n,ca,cb,cg".format(int(af2_struc), int(res)+1))
                        torsion_angle_list.append(torison_angle)
                    # Lysine
                    if wt_seq[res] == 'K':
                        torison_angle = run(session, "torsion #{}/A:{}@n,ca,cb,cg".format(int(af2_struc), int(res)+1))
                        torsion_angle_list.append(torison_angle)
                    # Methionine
                    if wt_seq[res] == 'M':
                        torison_angle = run(session, "torsion #{}/A:{}@n,ca,cb,cg".format(int(af2_struc), int(res)+1))
                        torsion_angle_list.append(torison_angle)
                    # Phenylalanine
                    if wt_seq[res] == 'F':
                        torison_angle = run(session, "torsion #{}/A:{}@n,ca,cb,cg".format(int(af2_struc), int(res)+1))
                        torsion_angle_list.append(torison_angle)
                    # Phenylalanine
                    if wt_seq[res] == 'P':
                        torison_angle = run(session, "torsion #{}/A:{}@n,ca,cb,cg".format(int(af2_struc), int(res)+1))
                        torsion_angle_list.append(torison_angle)
                    # Serine
                    if wt_seq[res] == 'S':
                        torison_angle = run(session, "torsion #{}/A:{}@n,ca,cb,og".format(int(af2_struc), int(res)+1))
                        torsion_angle_list.append(torison_angle)
                    # Threonine
                    if wt_seq[res] == 'T':
                        torison_angle = run(session, "torsion #{}/A:{}@n,ca,cb,og1".format(int(af2_struc), int(res)+1))
                        torsion_angle_list.append(torison_angle)
                    # Tryptophan
                    if wt_seq[res] == 'W':
                        torison_angle = run(session, "torsion #{}/A:{}@n,ca,cb,cg".format(int(af2_struc), int(res)+1))
                        torsion_angle_list.append(torison_angle)
                    # Tyrosine
                    if wt_seq[res] == 'Y':
                        torison_angle = run(session, "torsion #{}/A:{}@n,ca,cb,cg".format(int(af2_struc), int(res)+1))
                        torsion_angle_list.append(torison_angle)
                    # Valine
                    if wt_seq[res] == 'V':
                        torison_angle = run(session, "torsion #{}/A:{}@n,ca,cb,cg1".format(int(af2_struc), int(res)+1))
                        torsion_angle_list.append(torison_angle)
                    # Cannot calculate Chi1 for Alainine and Glycine
                    if wt_seq[res] == 'A' or 'G':
                        torsion_angle = np.nan
                        torsion_angle_list.append(torison_angle)
        except (chimerax.core.errors.UserError) as err:
            with open('error.log', 'a') as f:
                print("Can't calculate RMSD for {}:".format(i, models[i]), str(err), file=f);
                f.close()
        print(torsion_angle_list)
    # Save RMSD file
    df = pd.DataFrame(rmsd_list, columns=["name", "rmsd", "mut"])
    df.to_csv('/Users/holger/Desktop/master_thesis/data/results_rmsd/results_rmsd_af2/results_rmsd_P61626_mut/rmsd_mut_{}_{}.csv'.format(uniprot_id, af2_struc_list[af2_struc]), index=False)
    
    # Save torsion angle file
    df1 = pd.DataFrame(torsion_angle_list, columns=["chi_1",])
    df1.to_csv('/Users/holger/Desktop/master_thesis/data/results_torsion_angle/results_torsion_angle_af2/torsion_angle_mut_{}_{}.csv'.format(uniprot_id, af2_struc_list[af2_struc]), index=False)
    
    # Close session
    run(session, "close")

