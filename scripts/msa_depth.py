import pandas as pd
import numpy as np
import os, sys, argparse, json, time

if len(sys.argv) != 3:
    print("Usage: msa_depth.py <pdb_id> <n_lines>" )
    sys.exit(1)

try:
    pdb_id = str(sys.argv[1])
    n_lines = int(sys.argv[2])

except ValueError:
    print("Oops!  That was no valid number.  Try again...")

cwd = os.getcwd()
os.chdir('..')
DATA_PATH_ROOT = os.getcwd() + '/'
os.chdir(cwd)

def msa_depth(n_lines):
    n_seqs = n_lines / 2
    return n_seqs

def write_csv(n_seqs):
    #index = [i for i in range(n_dirs)]  
    d = {'pdb_id': [pdb_id], 'msa_depth': [n_seqs]}
    df = pd.DataFrame(data=d)
    return df.to_csv(DATA_PATH_ROOT + 'data/results_msa/P0A7Y4/msa_128_{}.csv'.format(pdb_id))

if __name__ == '__main__':
    n_seqs =  msa_depth(n_lines)
    write_csv(n_seqs)

