samtools mpileup --output-extra MAPQ ${snakemake_input[0]} 2> ${snakemake_log[0]} | awk -v OFS='\t' '{ if ($4 != 0) {
        sum = 0
        split($7, values, ",")
        for (i = 1; i <= $4; i++) sum += values[i]
        avg = sum / $4
        end = $2 + 1
        print $1,$2,end,avg 
    } else {
        end = $2 + 1
        print $1,$2,end,"NA"
    }
    }' > ${snakemake_output[0]} 2> ${snakemake_log[0]}

gunzip ${snakemake_input[1]}  --stdout 2>> ${snakemake_log[0]} | cut -f1,2,3 2>> ${snakemake_log[0]} | bedtools map -a - -b ${snakemake_output[0]} -c 4 -o mean -null 0 > ${snakemake_output[1]} 2>> ${snakemake_log[0]}