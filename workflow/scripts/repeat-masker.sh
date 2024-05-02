#!/bin/bash
set -e

if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <threads> <database_path> <fasta> <known> <unknown> <output>"
    exit 1
fi

threads=$1 
echo "Using" $threads "threads"

# Check if $threads is an integer
if ! [[ $threads =~ ^[0-9]+$ ]]; then
    echo "Error: $threads is not a valid integer"
    exit 1
fi

database_path=$2
echo "Using database in:" $database_path

# Check if $database_path exists
if [ ! -f "$database_path" ]; then
    echo "Error: $database_path does not exist"
    exit 1
fi

fasta=$3

# Check if $fasta exists
if [ ! -f "$fasta" ]; then
    echo "Error: $fasta does not exist"
    exit 1
fi
echo "Using fasta:" $fasta

lineage=$(basename "$fasta" .fasta)
known=$4

# Check if $known exists
if [ ! -f "$known" ]; then
    echo "Error: $known does not exist"
    exit 1
fi
echo "Using known repeats in:" $known

unknown=$5

# Check if $unknown exists
if [ ! -f "$unknown" ]; then
    echo "Error: $unknown does not exist"
    exit 1
fi
echo "Using unknown repeats in:" $unknown

output=$6
echo "Requested output file:" $output
repeat_dir=$(dirname "$output")
echo "Repeat directory:" $repeat_dir
echo "Getting absolute path of ${repeat_dir}"
repeat_dir_abs=$(realpath ${repeat_dir})
echo "Absolute path of ${repeat_dir}: ${repeat_dir_abs}"
database="database.fasta"
ln -s -r -f ${database_path} ${repeat_dir_abs}/${database}

# #### RepeatMasker ####
echo "Starting RepeatMasker"
echo "Creating output directories"
echo "mkdir -p ${repeat_dir_abs}/01_simple ${repeat_dir_abs}/02_complex ${repeat_dir_abs}/03_known ${repeat_dir_abs}/04_unknown"
mkdir -p ${repeat_dir_abs}/01_simple 
mkdir -p ${repeat_dir_abs}/02_complex 
mkdir -p ${repeat_dir_abs}/03_known 
mkdir -p ${repeat_dir_abs}/04_unknown

# round 1. Simple Repeats -a -noint -xsmall
echo "Identifying Simple Repeats"
echo "RepeatMasker -pa ${threads} -lib ${repeat_dir_abs}/${database} -a -e ncbi -dir ${repeat_dir_abs}/01_simple -noint -xsmall ${fasta}"
RepeatMasker -pa ${threads} -lib ${repeat_dir_abs}/${database} -a -e ncbi -dir ${repeat_dir_abs}/01_simple -noint -xsmall ${fasta}

# round 2: identify  elements sourced from Repbase using output from 1st round of RepeatMasker
echo "Identifying Complex Repeats"
echo "RepeatMasker -pa ${threads} -lib ${repeat_dir_abs}/${database} -a -e ncbi -dir ${repeat_dir_abs}/02_complex -nolow ${fasta}"
RepeatMasker -pa ${threads} -lib ${repeat_dir_abs}/${database} -a -e ncbi -dir ${repeat_dir_abs}/02_complex -nolow ${fasta}

# round 3: identify known elements
echo "Identifying Known Repeat Families"
echo "RepeatMasker -pa ${threads} -lib ${known} -a -e ncbi -dir ${repeat_dir_abs}/03_known -nolow  ${fasta}"
RepeatMasker -pa ${threads} -lib ${known} -a -e ncbi -dir ${repeat_dir_abs}/03_known -nolow  ${fasta}

# round 4: identify unknown elements
echo "Identifying Unknown Repeat Families"
echo "RepeatMasker -pa ${threads} -lib ${unknown} -a -e ncbi -dir ${repeat_dir_abs}/04_unknown -nolow  ${fasta}"
RepeatMasker -pa ${threads} -lib ${unknown} -a -e ncbi -dir ${repeat_dir_abs}/04_unknown -nolow  ${fasta}

#### Combine results of all rounds of masking ####
echo "Combining Results From Previous Rounds"

tail -n +4 ${repeat_dir_abs}/01_simple/${lineage}.fasta.out | awk '{print $5"\t"$6"\t"$7"\t"$11}' > ${repeat_dir_abs}/01_simple/${lineage}.bed
tail -n +4 ${repeat_dir_abs}/02_complex/${lineage}.fasta.out | awk '{print $5"\t"$6"\t"$7"\t"$11}' > ${repeat_dir_abs}/02_complex/${lineage}.bed
tail -n +4 ${repeat_dir_abs}/03_known/${lineage}.fasta.out | awk '{print $5"\t"$6"\t"$7"\t"$11}' > ${repeat_dir_abs}/03_known/${lineage}.bed
tail -n +4 ${repeat_dir_abs}/04_unknown/${lineage}.fasta.out | awk '{print $5"\t"$6"\t"$7"\t"$11}' > ${repeat_dir_abs}/04_unknown/${lineage}.bed

cat ${repeat_dir_abs}/01_simple/${lineage}.bed ${repeat_dir_abs}/02_complex/${lineage}.bed ${repeat_dir_abs}/03_known/${lineage}.bed ${repeat_dir_abs}/04_unknown/${lineage}.bed > ${output}.unsorted

bedtools sort -i ${output}.unsorted > ${output}

echo "Cleaning up"
find . -type d -empty -delete

# Check if $output file is empty, remove it and raise an error
if [ ! -s "$output" ]; then
    echo "Error: $output file is empty. Removing it."
    rm "$output"
    exit 1
fi
