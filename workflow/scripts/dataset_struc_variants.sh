# /usr/bin/bash
# This script is used to merge the structural variants from different samples into a single file

echo "Input files:" ${snakemake_input[*]} >> ${snakemake_log}
counter=0
for file in ${snakemake_input[*]}
do
    name=$(basename "$(dirname "$file")")
    echo "Processing" $name
    if [ $counter -eq 0 ]
    then
        awk -F'\t' 'NR==1 {OFS="\t"; print "sample", $1, $2, $3, $4, $5, $6, $7, $(NF-2), $(NF-1), $NF; exit}' "$file" > ${snakemake_output[0]} # print header
    fi
    awk -F'\t' -v name="$name" 'NR>1 {OFS="\t"; print name,  $1, $2, $3, $4, $5, $6, $7, $(NF-2), $(NF-1), $NF}' "$file" > "$file.tmp" # print data
    cat $file.tmp >> ${snakemake_output[0]} # merge data
    rm $file.tmp # remove temporary file
    counter=$((counter+1))
done >> ${snakemake_log}