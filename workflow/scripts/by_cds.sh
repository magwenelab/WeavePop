#!/usr/bin/env bash

mkdir results/cds
cat ${snakemake_input[0]} | while read protein
do
   parallel "seqkit faidx analysis/{}/predicted_cds.fa $protein | seqkit replace -p '($)' -r ' sample={}'"  :::: ${snakemake_input[1]} > results/cds/$protein.fa
done 2> ${snakemake_log[0]}
touch ${snakemake_output[0]}

# mkdir cds/
# cat files/protein_list.txt | while read protein
# do
#     parallel -j 4 "seqkit faidx analysis/{}/predicted_cds.fa $protein | seqkit replace -p '($)' -r ' sample={}'"  :::: samples.txt > cds/$protein.fa
# done 2> logs/cds/cds.log
# touch cds/done.txt