#!/bin/bash

threads=$1
echo $threads
lineage=$2
echo $lineage
repeat_dir=$3
echo "Making directory:" $repeat_dir
mkdir -p $repeat_dir
cd $repeat_dir
ln -s -r -f ../${lineage}.fasta .
echo "Working directory:" 
pwd

#### RepeatModeler ####
echo "Starting RepeatModeler"
mkdir ${lineage}_db 
echo "Building lineage database"
BuildDatabase -name ${lineage} -engine ncbi ../${lineage}.fasta
echo "Running RepeatModeler"
RepeatModeler -pa ${threads} -engine ncbi -database ${lineage}
echo "Moving database files"
mv ${lineage}.* ${lineage}_db 
echo "Separating known and unknown families"
cat ${lineage}-families.fa | seqkit fx2tab | grep -v "Unknown" | seqkit tab2fx > ${lineage}_known.fa
cat ${lineage}-families.fa | seqkit fx2tab | grep "Unknown" | seqkit tab2fx > ${lineage}_unknown.fa
