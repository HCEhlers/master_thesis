#! /bin/bash -f

## Define path to code directory
RDIR="/Users/holger/Desktop/master_thesis/scripts"

# Define PDBs to calculate RMSDs for
PDBS_array=('1gb5', '1gb6', '1gb7', '1gb8', '1gb9', '1gbo', '1gbw', '1gbx',
       '1gby', '1gbz', '1gfh', '1gfj', '1gfk', '1gfr', '1gft', '1gfu',
       '1gfv', '1inu', '1lhh', '1lhi', '1lhj', '1lhk', '1lhl', '1oub',
       '1ouc', '1oud', '1oue', '1ouh', '1oui', '1ouj', '1wqm', '1wqn',
       '1wqo', '1wqp', '1wqq', '1wqr', '1yam', '1yan', '1yao', '1yap',
       '1yaq', '2hea', '2heb', '2hec', '2hed', '2hee', '2hef')

for pdb in ${PDBS_array[@]}; do
    python ${RDIR}/exe_chimerax.py P61626 1lz1 ${pdb}
done
