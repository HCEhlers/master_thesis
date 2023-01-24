## Define path to data set
DIRPATH="/Users/holger/Desktop/master_thesis/data/results_af2/P0A7Y4/"

for dir in ${DIRPATH}/*
do 

PDB_ID=$(echo ${dir} | awk -F'[_.]' '{print $4}')

echo "Count sequences for: " ${PDB_ID} 

seq_count=$(wc -l < ${dir}/*.a3m | xargs)

echo ${seq_count}

python msa_depth.py $PDB_ID $dir_count $seq_count

done
