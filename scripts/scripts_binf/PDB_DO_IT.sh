#! /bin/bash -f

## Define path to code directory
RDIR="/storage1/holger/master_thesis/scripts"

## Define path to data set
DDIR="/storage1/holger/master_thesis/data/cytosol_data"

## Define main data path
DDIR_MAIN="/storage1/holger/master_thesis/data"

## Define file counter
let file_counter=1

## Create and go to results folder
mkdir ${DDIR_MAIN}/results_3000_cytosol
cd ${DDIR_MAIN}/results_3000_cytosol

# Get UniProt ID from cytosol.csv and loop through each UNIPROT ID that have an associated strucmap_*.csv (file with information about structures)
for file in ${DDIR}/strucmap_*.csv
do
UNIPROT_ID=$(echo ${file} | awk -F'[_.]' '{print $4}')

# Create and go to output folder
mkdir output_${UNIPROT_ID}
cd output_${UNIPROT_ID}

# Copy strucmaps and alignments to folder (not all UniProt IDs have alignment - if statement catches this)
cp ${DDIR}/strucmap_${UNIPROT_ID}.csv ${DDIR_MAIN}/results_3000_cytosol/output_${UNIPROT_ID}/strucmap_${UNIPROT_ID}.csv
if test -f "${DDIR}/alignments_${UNIPROT_ID}.fasta";	

then
cp ${DDIR}/alignments_${UNIPROT_ID}.fasta ${DDIR_MAIN}/results_3000_cytosol/output_${UNIPROT_ID}/alignments_${UNIPROT_ID}.fasta

# Run python script to get PDB IDs with mutations for first 1000 UniProt IDs ("Usage of uniprot_pdb_seq.py: uniprot_pdb_seq.py <UniProtID> <output.csv>" )
/storage1/cagiada/opt/anaconda3/bin/python ${RDIR}/uniprot_pdb_seq.py ${UNIPROT_ID} ${UNIPROT_ID}_mutations.csv
else 
echo "${UNIPROT_ID}.fasta was not found"
fi

cd ..

# Only run first 100 entries
let file_counter++
if [[ "${file_counter}" == '3000' ]]
then
echo "Have finished  ${file_counter} files!"
break
fi

# Print number of output
ls -1 | wc -l
done
