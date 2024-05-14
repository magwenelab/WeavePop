#!/bin/bash

threads=$1 
echo "Using" $threads "threads"
fasta=$2
echo "Fasta:" $fasta

# Check if the correct number of arguments were given
if [ $# -ne 3 ]; then
    echo "Usage: $0 <threads> <fasta> <repeats_dir>"
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
fasta_abs=$(realpath $fasta)
work_dir=$(dirname "$fasta")
echo "Moving to working directory:" $work_dir
cd $work_dir
repeats_dir=$3
echo "Repeats directory:" $repeats_dir
mkdir -p ${repeats_dir}
repeats_dir_abs=$(realpath ${repeats_dir})
#### RepeatModeler ####
echo "Starting RepeatModeler"
echo "Building lineage database"
mkdir -p ${repeats_dir}/${lineage}_db/
BuildDatabase -name ${repeats_dir}/${lineage}_db/${lineage} -engine ncbi ${fasta_abs}
echo "Finished Building lineage database"

rmodeler_dir=RModeler
mkdir -p ${repeats_dir}/${rmodeler_dir}
echo "Getting absolute path of ${repeats_dir}/${rmodeler_dir}"
rmodeler_dir_abs=$(realpath ${repeats_dir}/${rmodeler_dir})
echo "Absolute path of ${repeats_dir}/${rmodeler_dir}: ${rmodeler_dir_abs}"
echo "Running RepeatModeler"
RepeatModeler -pa ${threads} -engine ncbi -database ${repeats_dir}/${lineage}_db/${lineage} -dir ${rmodeler_dir_abs}

echo "Separating known and unknown families"
cat ${repeats_dir}/${lineage}_db/${lineage}-families.fa | seqkit fx2tab | grep -v "Unknown" | seqkit tab2fx > ${repeats_dir}/${lineage}_known.fa
cat ${repeats_dir}/${lineage}_db/${lineage}-families.fa | seqkit fx2tab | grep "Unknown" | seqkit tab2fx > ${repeats_dir}/${lineage}_unknown.fa
