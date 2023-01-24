import pandas as pd
import numpy as np
import os, sys, argparse, json, time
import urllib.request
import requests


if len(sys.argv) !=3:
    print("Usage: find_wt_structures.py <UniProtID> <uniprot.fasta>" )
    sys.exit(1)

try:
    uniprot_id = sys.argv[1]
    outfile = sys.argv[2]

except (IOError, FileNotFoundError) as err:
    print("Cant open file:", str(err));
    #sys.exit(1)

def get_pdb_seqs(pdb_ids):
    seq_list = [] 
    for pdb in range(len(pdb_ids)):
        try:
            url = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/%s/' % pdb_ids[pdb]
            df = pd.read_json(url)
            seq_struc = df[pdb_ids[pdb]][0]['pdb_sequence'] 
            seq_list.append(seq_struc)
        except (TypeError, KeyError) as e:
            print("No PDBs found")
    return seq_list


def write_fasta(pdb,seq):
    with open(outfile, 'a') as f:
        f.write('>' + pdb + '\n')
        f.write(seq + '\n')

if __name__ == '__main__':
    pdbs= ['1hem', '1heo', '1her']
    seq_list = get_pdb_seqs(pdbs)
    
    for i in range(len(seq_list)):
        print(pdbs[i],str(seq_list[i]))
        write_fasta(pdbs[i],str(seq_list[i]))

