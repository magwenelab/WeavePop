
# Concatenate good coverage tables of all samples into one
cat ${snakemake_input[g]} | head -n 1 1> ${snakemake_output[allg]} 2> ${snakemake_log[0]}
tail -q -n +2 ${snakemake_input[g]}  1>> ${snakemake_output[allg]} 2>> ${snakemake_log[0]}
# Concatenate raw coverage tables of all samples into one
cat ${snakemake_input[r]}  | head -n 1 1> ${snakemake_output[allr]} 2>> ${snakemake_log[0]}
tail -q -n +2 ${snakemake_input[r]} 1>> ${snakemake_output[allr]} 2>> ${snakemake_log[0]}
# Concatenate structural variants tables of all samples into one
cat ${snakemake_input[sv]} | head -n 1 1> ${snakemake_output[allsv]} 2> ${snakemake_log[0]}
tail -q -n +2 ${snakemake_input[sv]}  1>> ${snakemake_output[allsv]} 2>> ${snakemake_log[0]}
# Concatenate mapped reads tables of all samples into one
cat ${snakemake_input[m]} | head -n 1 1> ${snakemake_output[allm]} 2> ${snakemake_log[0]}
tail -q -n +2 ${snakemake_input[m]}  1>> ${snakemake_output[allm]} 2>> ${snakemake_log[0]}
# Concatenate mapqcov tables of all samples into one
cat ${snakemake_input[mc]} | head -n 1 1> ${snakemake_output[allmc]} 2> ${snakemake_log[0]}
tail -q -n +2 ${snakemake_input[mc]}  1>> ${snakemake_output[allmc]} 2>> ${snakemake_log[0]}