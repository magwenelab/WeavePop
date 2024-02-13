lineage=$1

# Download RepBase database ####
wget https://www.girinst.org/server/RepBase/protected/RepBase29.01.fasta.tar.gz
tar -xvzf RepBase29.01.fasta.tar.gz
cat RepBase29.01.fasta/*.ref > RepBase.fasta
cat RepBase29.01.fasta/appendix/*.ref >> RepBase.fasta
rm -rf RepBase29.01.fasta/ RepBase29.01.fasta.tar.gz 


# RepeatModeler ####
mkdir logs ${lineage}_db 
BuildDatabase -name Crypto-${lineage} -engine ncbi ${lineage}.fasta &> logs/00_repeatmodeler.log
RepeatModeler -pa 36 -engine ncbi -database Crypto-${lineage} &> logs/01_repeatmodeler.log
mv Crypto-${lineage}.* ${lineage}_db 
cat Crypto-${lineage}-families.fa | seqkit fx2tab | awk '{ print "cryNeo1_"$0 }' | seqkit tab2fx > Crypto-${lineage}-families.prefix.fa # Add species prefix to headers
cat Crypto-${lineage}-families.prefix.fa | seqkit fx2tab | grep -v "Unknown" | seqkit tab2fx > Crypto-${lineage}_known.fa # Separate known repeat sequences
cat Crypto-${lineage}-families.prefix.fa | seqkit fx2tab | grep "Unknown" | seqkit tab2fx > Crypto-${lineage}_unknown.fa # Separate unknown repeat sequences

# RepeatMasker ####
mkdir -p 01_simple_out 02_complex_out 03_known_out 04_unknown_out

# round 1. Simple Repeats -a -noint -xsmall
RepeatMasker -pa 36 -lib RepBase.fasta -a -e ncbi -dir 01_simple_out -noint -xsmall ${lineage}.fasta &> logs/01_simplemask.log
mv 01_simple_out/${lineage}.fasta.masked 01_simple_out/${lineage}
mv 01_simple_out/${lineage}.fasta.out 01_simple_out/${lineage}.out
for i in align tbl ori.out cat.gz; do mv 01_simple_out/*.$i 01_simple_out/${lineage}.$i; done


# round 2: annotate/mask fungi elements sourced from Repbase using output from 1st round of RepeatMasker
RepeatMasker -pa 36 -lib RepBase.fasta -a -e ncbi -dir 02_complex_out -nolow 01_simple_out/${lineage} &> logs/02_complex.log
mv 02_complex_out/${lineage}.masked 02_complex_out/${lineage}

# round 3: annotate/mask known elements sourced from species-specific de novo repeat library using output froom 2nd round of RepeatMasker
RepeatMasker -pa 36 -lib Crypto-${lineage}_known.fa -a -e ncbi -dir 03_known_out -nolow  02_complex_out/${lineage} &> logs/03_known.log
mv 03_known_out/${lineage}.masked 03_known_out/${lineage}

# round 4: annotate/mask unknown elements sourced from species-specific de novo repeat library using output from 3nd round of RepeatMasker
RepeatMasker -pa 36 -lib Crypto-${lineage}_unknown.fa -a -e ncbi -dir 04_unknown_out -nolow  03_known_out/${lineage} &> logs/04_unknown.log
mv 04_unknown_out/${lineage}.masked 04_unknown_out/${lineage}

# Combine results of all rounds of masking ####
mkdir 05_full_out
cat *_out/${lineage}.cat.gz > 05_full_out/${lineage}.full_mask.cat.gz
cat 01_simple_out/${lineage}.out <(cat 02_complex_out/${lineage}.out | tail -n +4) <(cat 03_known_out/${lineage}.out | tail -n +4) <(cat 04_unknown_out/${lineage}.out | tail -n +4) > 05_full_out/${lineage}.full_mask.out
cat */${lineage}.align > 05_full_out/${lineage}.full_mask.align
# resummarize repeat compositions from combined analysis of all RepeatMasker rounds
ProcessRepeats -a -species fungi 05_full_out/${lineage}.full_mask.cat.gz &> logs/05_fullmask.log

# Create bed file for simple repeat regions
tail -n +4 01_simple_out/${lineage}.out | awk '{print $5"\t"$6"\t"$7}' > 01_simple_out/${lineage}.bed
# Create bed file for regions from round 2, 3 and 4
cat 02_complex_out/${lineage}.out <(cat 03_known_out/${lineage}.out | tail -n +4) <(cat 04_unknown_out/${lineage}.out | tail -n +4) | tail -n +4 | awk '{print $5"\t"$6"\t"$7}' > 05_full_out/${lineage}.complex_mask.bed
#tail -n +4 05_full_out/${lineage}.full_mask.out | awk '{print $5"\t"$6"\t"$7}' > 05_full_out/${lineage}.full_mask.bed

# Mask fasta first softmask with simple repeats, then hard mask with the other repeats
bedtools maskfasta -fullHeader -soft -fi ${lineage}.fasta -bed 01_simple_out/${lineage}.bed -fo 05_full_out/${lineage}.simple_mask.fasta
bedtools maskfasta -fullHeader -fi 05_full_out/${lineage}.simple_mask.fasta -bed 05_full_out/${lineage}.complex_mask.bed -fo 05_full_out/${lineage}.full_mask.fasta