
## Inputs 

Input to each pipeline is tabular (CSV), typically `sample_table.csv`. First column always `sample`, other columns will vary depending on pipeline.

* For Snippy, columns might be: `sample, refgenome, fq1, fq2`

* For liftoff, columns might be:  `sample, samplegenome, refgenome, refgff`


Pipelines operate at a sample level.

Pipelines will write output to `OUTDIR/SAMPLE/TOOL/`

* For example, Snippy analysis of PMY1234 will write output to `./analysis/PMY1234/snippy`

Other than OUTDIR, avoid hard coding or requiring path info in config file, rather build scripts that generate appropriate paths into `sample_table.csv`. 

* snippy line:  

    ```
    PMY1234, refgenomes/H99_genome.fasta, fastqs/PMY1234_r1.fq.gz, fastqs/PMY1234_r1.fq.gz
    ```

    In the above `refgenomes/` and `fastqs/` are not part of the config. Instead, build scripts to generate appropriate table entries




