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
mkdir -p 01_simple_out 02_complex_out 03_known_out 04_unknown_out 05_full_out

# round 1. Simple Repeats -a -noint -xsmall
echo "Identifying Simple Repeats"
RepeatMasker -pa ${threads} -lib ${database} -a -e ncbi -dir 01_simple_out -noint -xsmall ${lineage}.fasta
mv 01_simple_out/${lineage}.fasta.masked 01_simple_out/${lineage}
mv 01_simple_out/${lineage}.fasta.out 01_simple_out/${lineage}.out
for i in align tbl ori.out cat.gz; do mv 01_simple_out/*.$i 01_simple_out/${lineage}.$i; done

# round 2: annotate/mask fungi elements sourced from Repbase using output from 1st round of RepeatMasker
echo "Identifying Complex Repeats"
RepeatMasker -pa ${threads} -lib ${database} -a -e ncbi -dir 02_complex_out -nolow 01_simple_out/${lineage}
mv 02_complex_out/${lineage}.masked 02_complex_out/${lineage}

# round 3: annotate/mask known elements sourced from species-specific de novo repeat library using output froom 2nd round of RepeatMasker
echo "Identifying Known Repeat Families"
RepeatMasker -pa ${threads} -lib ${lineage}_known.fa -a -e ncbi -dir 03_known_out -nolow  02_complex_out/${lineage}
mv 03_known_out/${lineage}.masked 03_known_out/${lineage}

# round 4: annotate/mask unknown elements sourced from species-specific de novo repeat library using output from 3nd round of RepeatMasker
echo "Identifying Unknown Repeat Families"
RepeatMasker -pa ${threads} -lib ${lineage}_unknown.fa -a -e ncbi -dir 04_unknown_out -nolow  03_known_out/${lineage}
mv 04_unknown_out/${lineage}.masked 04_unknown_out/${lineage}

#### Combine results of all rounds of masking ####
echo "Combining Results From Previous Rounds"
cat *_out/${lineage}.cat.gz > 05_full_out/${lineage}.full_mask.cat.gz
cat 01_simple_out/${lineage}.out <(cat 02_complex_out/${lineage}.out | tail -n +4) <(cat 03_known_out/${lineage}.out | tail -n +4) <(cat 04_unknown_out/${lineage}.out | tail -n +4) > 05_full_out/${lineage}.full_mask.out
cat */${lineage}.align > 05_full_out/${lineage}.full_mask.align
# resummarize repeat compositions from combined analysis of all RepeatMasker rounds
ProcessRepeats -a -species fungi 05_full_out/${lineage}.full_mask.cat.gz

# Create bed file for simple repeat regions
tail -n +4 01_simple_out/${lineage}.out | awk '{print $5"\t"$6"\t"$7}' > 01_simple_out/${lineage}.bed
# Create bed file for regions from round 2, 3 and 4
cat 02_complex_out/${lineage}.out <(cat 03_known_out/${lineage}.out | tail -n +4) <(cat 04_unknown_out/${lineage}.out | tail -n +4) | tail -n +4 | awk '{print $5"\t"$6"\t"$7"\t"$11}' > 05_full_out/${lineage}.complex_mask.bed
tail -n +4 05_full_out/${lineage}.full_mask.out | awk '{print $5"\t"$6"\t"$7"\t"$11}' > 05_full_out/${lineage}.full_mask.bed

# Mask fasta first softmask with simple repeats, then hard mask with the other repeats
echo "Creating Masked Reference Fasta"
bedtools maskfasta -fullHeader -soft -fi ${lineage}.fasta -bed 01_simple_out/${lineage}.bed -fo 05_full_out/${lineage}.simple_mask.fasta
bedtools maskfasta -fullHeader -fi 05_full_out/${lineage}.simple_mask.fasta -bed 05_full_out/${lineage}.complex_mask.bed -fo 05_full_out/${lineage}.full_mask.fasta


