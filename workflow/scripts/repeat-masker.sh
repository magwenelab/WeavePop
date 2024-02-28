#!/bin/bash

threads=$1 
echo $threads
lineage=$2
echo $lineage
repeat_dir=$3
echo "Using databse in:" $4
database_path=$4
database="database.fasta"
ln -s -r -f ${database_path} ${repeat_dir}/${database}
cd $repeat_dir
ln -s -r -f ../${lineage}.fasta .
echo "Working directory:" 
pwd

# database_path=config/RepBase.fasta
# lineage=VNI
# # mkdir -p results/references/${lineage}/repeats
# ln -s -r -f ${database_path} results/references/${lineage}/repeats/database.fasta
# database="database.fasta"
# cd results/references/${lineage}/repeats
# ln -s -r -f ../${lineage}.fasta .

# Download RepBase database ####
# wget https://www.girinst.org/server/RepBase/protected/RepBase29.01.fasta.tar.gz
# tar -xvzf RepBase29.01.fasta.tar.gz
# cat RepBase29.01.fasta/*.ref > RepBase.fasta
# cat RepBase29.01.fasta/appendix/*.ref >> RepBase.fasta
# rm -rf RepBase29.01.fasta/ RepBase29.01.fasta.tar.gz 

# #### RepeatMasker ####
echo "Starting RepeatMasker"
mkdir -p 01_simple 02_complex 03_known 04_unknown 05_full

# round 1. Simple Repeats -a -noint -xsmall
echo "Identifying Simple Repeats"
RepeatMasker -pa ${threads} -lib ${database} -a -e ncbi -dir 01_simple -noint -xsmall ${lineage}.fasta
# mv 01_simple/${lineage}.fasta.masked 01_simple/${lineage}
# mv 01_simple/${lineage}.fasta.out 01_simple/${lineage}.out
# for i in align tbl ori.out cat.gz; do mv 01_simple/*.$i 01_simple/${lineage}.$i; done

# round 2: annotate/mask fungi elements sourced from Repbase using output from 1st round of RepeatMasker
echo "Identifying Complex Repeats"
RepeatMasker -pa ${threads} -lib ${database} -a -e ncbi -dir 02_complex -nolow ${lineage}.fasta
#mv 02_complex/${lineage}.masked 02_complex/${lineage}

# round 3: annotate/mask known elements sourced from species-specific de novo repeat library using output from 2nd round of RepeatMasker
echo "Identifying Known Repeat Families"
RepeatMasker -pa ${threads} -lib ${lineage}_known.fa -a -e ncbi -dir 03_known -nolow  ${lineage}.fasta
# mv 03_known/${lineage}.masked 03_known/${lineage}

# round 4: annotate/mask unknown elements sourced from species-specific de novo repeat library using output from 3nd round of RepeatMasker
echo "Identifying Unknown Repeat Families"
RepeatMasker -pa ${threads} -lib ${lineage}_unknown.fa -a -e ncbi -dir 04_unknown -nolow  ${lineage}.fasta
# mv 04_unknown/${lineage}.masked 04_unknown/${lineage}

#### Combine results of all rounds of masking ####
echo "Combining Results From Previous Rounds"
# cat */${lineage}.cat.gz > 05_full/${lineage}.full_mask.cat.gz
# cat 01_simple/${lineage}.out <(cat 02_complex/${lineage}.out | tail -n +4) <(cat 03_known/${lineage}.out | tail -n +4) <(cat 04_unknown/${lineage}.out | tail -n +4) > 05_full/${lineage}.full_mask.out
# cat */${lineage}.align > 05_full/${lineage}.full_mask.align
# resummarize repeat compositions from combined analysis of all RepeatMasker rounds
# ProcessRepeats -a -species fungi 05_full/${lineage}.full_mask.cat.gz

# Create bed file for simple repeat regions
tail -n +4 01_simple/${lineage}.fasta.out | awk '{print $5"\t"$6"\t"$7"\t"$11}' > 01_simple/${lineage}.bed
tail -n +4 02_complex/${lineage}.fasta.out | awk '{print $5"\t"$6"\t"$7"\t"$11}' > 02_complex/${lineage}.bed
tail -n +4 03_known/${lineage}.fasta.out | awk '{print $5"\t"$6"\t"$7"\t"$11}' > 03_known/${lineage}.bed
tail -n +4 04_unknown/${lineage}.fasta.out | awk '{print $5"\t"$6"\t"$7"\t"$11}' > 04_unknown/${lineage}.bed

cat 01_simple/${lineage}.bed 02_complex/${lineage}.bed 03_known/${lineage}.bed 04_unknown/${lineage}.bed > 05_full/${lineage}_unsorted.bed

bedtools sort -i 05_full/${lineage}_unsorted.bed > 05_full/${lineage}.bed

# Create bed file for regions from round 2, 3 and 4
# cat 02_complex/${lineage}.out <(cat 03_known/${lineage}.out | tail -n +4) <(cat 04_unknown/${lineage}.out | tail -n +4) | tail -n +4 | awk '{print $5"\t"$6"\t"$7"\t"$11}' > 05_full/${lineage}.complex_mask.bed
# tail -n +4 05_full/${lineage}.full_mask.out | awk '{print $5"\t"$6"\t"$7"\t"$11}' > 05_full/${lineage}.full_mask.bed

# # Mask fasta first softmask with simple repeats, then hard mask with the other repeats
# echo "Creating Masked Reference Fasta"
# bedtools maskfasta -fullHeader -soft -fi ${lineage}.fasta -bed 01_simple/${lineage}.bed -fo 05_full/${lineage}.simple_mask.fasta
# bedtools maskfasta -fullHeader -fi 05_full/${lineage}.simple_mask.fasta -bed 05_full/${lineage}.complex_mask.bed -fo 05_full/${lineage}.full_mask.fasta


