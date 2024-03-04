
cat ${snakemake_input[g]} | head -n 1 1> ${snakemake_output[allg]} 2> ${snakemake_log[0]}  && \
cat ${snakemake_input[r]}  | head -n 1 1> ${snakemake_output[allr]} 2>> ${snakemake_log[0]}   && \
tail -q -n +2 ${snakemake_input[g]}  1>> ${snakemake_output[allg]} 2>> ${snakemake_log[0]} && \
tail -q -n +2 ${snakemake_input[r]} 1>> ${snakemake_output[allr]} 2>> ${snakemake_log[0]} 