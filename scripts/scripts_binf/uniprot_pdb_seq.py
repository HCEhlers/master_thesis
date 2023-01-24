import pandas as pd
import numpy as np
import os, sys, argparse, json, time
import urllib.request
import requests


if len(sys.argv) != 3:
    print("Usage: uniprot_pdb_seq.py <UniProtID> <output.csv>" )
    sys.exit(1)

try:
    DATA_PATH = '/storage1/holger/master_thesis/data/cytosol_data/'
    uniprot_id = sys.argv[1]
    ALIGNMENT_FILE = open(DATA_PATH + 'alignments_' + uniprot_id + '.fasta', 'r')
    outfile = sys.argv[2]

except (IOError, FileNotFoundError) as err:
    print("Cant open file:", str(err));
    #sys.exit(1)

def return_strucmap_df():
   """
   Return a Pandas DataFrame from strucmap_*.csv file
   """
   DF_STRUCMAP = pd.read_csv(DATA_PATH + 'strucmap_' + uniprot_id + '.csv', delimiter=';')
   return DF_STRUCMAP

def get_pdbs_with_coverage(df):
    """
    Gets PDBs from strucmap files that have coverage > 0.7 into a list.
    Parameters
    ----------
    df1 : Pandas DataFrame
    """
    pdb_ids = list()
    # Get entries that have coverage > 0.7
    for i in range(len(df['pdb'].unique())):
        if df['cover_exp'][i] >= 0.7:
            pdb_ids.append(df['pdb'].unique()[i])
    return pdb_ids

def get_pdb_criteria_base(pdb_ids):
    """
    This function calls the PDB URL and checks the base criteria (not as restrictive for now...)
    for a PDB as defined in DOI 10.1002/prot.24073. Return list of pdbs to be REMOVED.
    Parameters
    ----------
    pdb_ids : list of strings
    """
    criteria = ['experimental_method', 'resolution']
    pdbs_rm = list()

    for pdb in range(len(pdb_ids)):
        url = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/experiment/%s/' % pdb_ids[pdb]
        df = pd.read_json(url)
        
        # Is the experimental method X-ray diffraction?
        try:
            if df[pdb_ids[pdb]][0][criteria[0]] != 'X-ray diffraction':
                print("The {} entry is not determined by X-ray diffraction ".format(pdb_ids[pdb]), 'and is excluded.')
                pdbs_rm.append(pdb_ids[pdb])
        except TypeError:
            print("No PDBs found")

        # Is resolution > 3.0? 
        try:
            if df[pdb_ids[pdb]][0][criteria[1]] > 3.0:
                print("The {} entry is above 3.0 Ã in resolution ".format(pdb_ids[pdb]), 'and is excluded.')
                pdbs_rm.append(pdb_ids[pdb])
        except (TypeError, KeyError) as e:
            print('Resolution for: {}'.format(pdb_ids[pdb]), "is missing")
            pdbs_rm.append(pdb_ids[pdb])

    return pdbs_rm

def get_mutations(ALIGNMENT_FILE, pdbs):
    """
    This functions looks at pdbs that fulfill the criteria and then
    opens the alignment files to find mutations.
     ----------
    alignment_file : fasta file with sequences for each pdb
    pdb_ids: list of pdb ids
    """
    first_loop = True
    entire_seqlist = list()
    seqlist = list()
    pdblist = list()
    chain_list = list()
    mutations = list()
    for line in ALIGNMENT_FILE:
        if line.startswith('>'):
            uniprotid = line.rstrip()[1:7]
            break
    
    for line in ALIGNMENT_FILE:
        if line.startswith('>'):
            if first_loop:
                wt_seq = ''.join(seqlist).upper()
                seqlist = []
                header = line.rstrip()[1:5]
                chain = line.rstrip()[6] 
                pdblist.append(header)
                chain_list.append(chain)
                first_loop = False
            else:
                seq = ''.join(seqlist).upper()
                entire_seqlist.append(seq)
                seqlist = []
                header = line.rstrip()[1:5]
                chain = line.rstrip()[6]
                pdblist.append(header)
                chain_list.append(chain)
        else:
            seqlist.append(''.join(line.split()))
   
    # Do last sequence
    seq = ''.join(seqlist).upper()
    seqlist = []

    # Find mutations
    for i in range(len(entire_seqlist)):
        if pdblist[i] in pdbs:
            # Find mismatches
            point_mutations = [mut for mut in range(len(wt_seq)) if (wt_seq[mut] != entire_seqlist[i][mut] and entire_seqlist[i][mut] != '-')]
            for mut in range(len(wt_seq)):
                if mut in point_mutations:
                    if len(point_mutations) == 1:
                        mutations.append([pdblist[i], chain_list[i], "".join([wt_seq[mut].strip(), str(mut+1),entire_seqlist[i][mut].strip()]), entire_seqlist[i]])
                    # We do not look at consecutive mutants in beginning , middle  or end  of the sequence for chrystallographic reasons
                    elif len(point_mutations) > 1 and point_mutations != list(range(min(point_mutations), max(point_mutations)+1)) and entire_seqlist[i][mut:mut+2] != 'XX' and entire_seqlist[i][mut-1] != 'X':
                        mutations.append([pdblist[i], chain_list[i], "".join([wt_seq[mut].strip(), str(mut+1),entire_seqlist[i][mut].strip()]), entire_seqlist[i]])

    print(uniprotid)      
    ALIGNMENT_FILE.close()    
    return mutations, wt_seq

def write_csv(mutations, outfile):
    """
    Takes input a list of mutations and different info from the associated UniProtID
    and writes to .csv file with PDB ID, chain, seq and mutation.
    ----------
    pdbs : list of strings
    sys.argv[2]: csv file to write
    """
    index = [i for i in range(len(mutations))]
    columns = ['pdb_id','chain', 'mutation', 'sequence']
    df = pd.DataFrame(mutations, columns=columns, index=index)
    df = pd.DataFrame(df.groupby(['pdb_id', 'chain', 'sequence'])['mutation'].apply(lambda x:'\n'.join(x)))
    return df.to_csv(outfile)

if __name__ == '__main__':
   df = return_strucmap_df()
   pdb_ids = get_pdbs_with_coverage(df)
   pdbs_rm = get_pdb_criteria_base(pdb_ids)
   pdbs = [pdb for pdb in pdb_ids if pdb not in pdbs_rm]
   mutations, wt_seq = get_mutations(ALIGNMENT_FILE, pdbs)
   print(mutations)
   write_csv(mutations, outfile)
   # Write potential good proteins for analysis to file.
   if len(mutations) > 60:
       with open('output.txt', 'a') as f:
           print(outfile, file=f)
