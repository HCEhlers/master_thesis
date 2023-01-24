import pandas as pd
import numpy as np
import os, sys, argparse, json, time
import urllib.request
import requests


if len(sys.argv) != 4:
    print("Usage: find_wt_structures.py <UniProtID> <wt_seq> <output.csv>" )
    sys.exit(1)

try:
    DATA_PATH = '/storage1/holger/master_thesis/data/output_P00698/'
    uniprot_id = sys.argv[1]
    wt_seq = str(sys.argv[2])
    outfile = sys.argv[3]

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
    Gets PDBs from strucmap files that have coverage = 1.0 into a list.
    Parameters
    ----------
    df1 : Pandas DataFrame
    """
    pdb_ids = list()
    # Get entries that have coverage = 1.0
    for i in range(len(df['pdb'].unique())):
        if df['cover_exp'][i] == 1.0:
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
    criteria = ['experimental_method', 'resolution', 'pdb_sequence']
    pdb_list = list()
    seq_list =list()
    
    for pdb in range(len(pdb_ids)):
        url = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/experiment/%s/' % pdb_ids[pdb]
        df = pd.read_json(url)
        fulfill = True
        
        # Is the experimental method X-ray diffraction?
        try:
            if df[pdb_ids[pdb]][0][criteria[0]] != 'X-ray diffraction':
                print("The {} entry is not determined by X-ray diffraction ".format(pdb_ids[pdb]), 'and is excluded.')
                fulfill = False
        except (TypeError, KeyError) as e:
            print("No PDBs found")
            fulfill = False
        # Is resolution < 2.5, then it is excluded? 
        try:
            if df[pdb_ids[pdb]][0][criteria[1]] > 2.5:
                print("The {} entry is above 2.5 in resolution ".format(pdb_ids[pdb]), 'and is excluded.')
                fulfill = False
        except (TypeError, KeyError) as e:
            print('Resolution for: {}'.format(pdb_ids[pdb]), "is missing")
            fulfill = False
        # Is the structure wild type 
        try:
            url = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/%s/' % pdb_ids[pdb]
            df = pd.read_json(url)
            seq_struc = df[pdb_ids[pdb]][0][criteria[2]] 
            if wt_seq != seq_struc:
                print("The {} entry is not identical to wt seq ".format(pdb_ids[pdb]), 'and is excluded.')
                fulfill = False
        except (TypeError, KeyError) as e:
            print("No PDBs found")
            fulfill = False

        # Get the sequences of structures
        if fulfill:
            pdb_list.append(pdb_ids[pdb])
            seq_list.append(seq_struc)
    return pdb_list, seq_list


def write_csv(pdbs,seqs, outfile):
    """
    Takes input a list of mutations and different info from the associated UniProtID
    and writes to .csv file with PDB ID, chain, seq and mutation.
    ----------
    pdbs : list of strings
    sys.argv[2]: csv file to write
    """
    index = [i for i in range(len(pdbs))]
    columns = ['pdb_id', 'seq']
    df = pd.DataFrame({columns[0]: pdbs,
        columns[1]: seqs})
    return df.to_csv(outfile)

if __name__ == '__main__':
   df = return_strucmap_df()
   pdb_ids = get_pdbs_with_coverage(df)
   pdb_list, seq_list = get_pdb_criteria_base(pdb_ids)
   write_csv(pdb_list, seq_list, outfile)

