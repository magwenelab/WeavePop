#/usr/bin/env bash
  
cut -d',' -f2 files/lineage_references.csv | tail -n +2 | while read file
    do 
        lineage=$(echo $file | cut -d'.' -f1)
        grep ">" references/$file | cut -d',' -f1 | cut -d'>' -f2 |  while read line
        do 
            accession=$(echo $line | cut -d' ' -f1)
            name=$(echo $line | rev | cut -d' ' -f1 | rev)
            echo $lineage,$accession,$name
        done
    done > files/chromosome_names.csv