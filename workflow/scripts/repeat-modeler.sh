#!/bin/bash

threads=$1 
echo "Using" $threads "threads"
fasta=$2
echo "Fasta:" $fasta

# Check if the correct number of arguments were given
if [ $# -ne 3 ]; then
    echo "Usage: $0 <threads> <fasta> <outdir>"
    exit 1
fi

# Check if the threads argument is an integer number
if ! [[ $threads =~ ^[0-9]+$ ]]; then
    echo "Threads argument must be an integer"
    exit 1
fi

# Check if the fasta file exists
if [ ! -f "$fasta" ]; then
    echo "Fasta file does not exist"
    exit 1
fi 

lineage=$(basename "$fasta" .fasta)
dir=$(dirname "$fasta")
outdir=$3
echo "Making directory:" ${dir}/${outdir}
mkdir -p ${dir}/${outdir}


#### RepeatModeler ####
echo "Starting RepeatModeler"
echo "Building lineage database"
mkdir -p ${dir}/${outdir}/${lineage}_db/
BuildDatabase -name ${dir}/${outdir}/${lineage}_db/${lineage} -engine ncbi ${fasta}

echo "Running RepeatModeler"
RepeatModeler -pa ${threads} -engine ncbi -database ${dir}/${outdir}/${lineage}_db/${lineage} & PID=$!
wait ${PID}
echo $PID
echo "Separating known and unknown families"
cat ${dir}/${outdir}/${lineage}_db/${lineage}-families.fa | seqkit fx2tab | grep -v "Unknown" | seqkit tab2fx > ${dir}/${outdir}/${lineage}_known.fa
cat ${dir}/${outdir}/${lineage}_db/${lineage}-families.fa | seqkit fx2tab | grep "Unknown" | seqkit tab2fx > ${dir}/${outdir}/${lineage}_unknown.fa

echo "Cleaning up"
rm -r RM_${PID}* | true