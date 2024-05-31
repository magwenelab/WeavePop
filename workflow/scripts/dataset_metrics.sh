

echo "Concatenate copy number variants table"
cat ${snakemake_input[c]} | head -n 1 1> ${snakemake_output[allc]} 2> ${snakemake_log[0]}
tail -q -n +2 ${snakemake_input[c]}  1>> ${snakemake_output[allc]} 2>> ${snakemake_log[0]}
echo "Concatenate mapping stats"
cat ${snakemake_input[m]} | head -n 1 1> ${snakemake_output[allm]} 2> ${snakemake_log[0]}
tail -q -n +2 ${snakemake_input[m]}  1>> ${snakemake_output[allm]} 2>> ${snakemake_log[0]}
echo "Concatenate depth_by_chrom_good_normalized"
cat ${snakemake_input[n]} | head -n 1 1> ${snakemake_output[alln]} 2> ${snakemake_log[0]}
tail -q -n +2 ${snakemake_input[n]}  1>> ${snakemake_output[alln]} 2>> ${snakemake_log[0]}
echo "Concatenate mapq_depth"
cat ${snakemake_input[md]} | head -n 1 1> ${snakemake_output[allmd]} 2> ${snakemake_log[0]}
tail -q -n +2 ${snakemake_input[md]}  1>> ${snakemake_output[allmd]} 2>> ${snakemake_log[0]}
