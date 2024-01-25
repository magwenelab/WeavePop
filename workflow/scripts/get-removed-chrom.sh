#/bin/bash
#
# Remove chromosome/contig from reference genome (fasta and gff)
# Keeping the original as name.fasta.original and name.gff.original

scriptargs="path-to-fasta path-to-gff seq_id"
nexpargs=3
nargs=$#

if [[ $nargs -ne $nexpargs ]]
then
    echo
    echo "Usage: $(basename $0) $scriptargs"
    echo
    exit
fi

fasta=$1
gff=$2
seqid=$3

echo "Removing $seqid from $fasta"
seqkit grep -n -r -v -p ${seqid} ${fasta} > ${fasta}.modif
mv ${fasta} ${fasta}.original
mv ${fasta}.modif ${fasta}

echo "Removing $seqid from $gff"
grep -v $seqid $gff >${gff}.modif
mv ${gff} ${gff}.original
mv ${gff}.modif ${gff}
