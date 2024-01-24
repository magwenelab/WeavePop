
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
    PMY1234, refgenomes/H99_genome.fasta, fastqs/PMY1234_r1.fq.gz, fastqs/PMY1234_r2.fq.gz
    ```

    In the above `refgenomes/` and `fastqs/` are not part of the config. Instead, build scripts to generate appropriate table entries






## 

Metadata File:

* sample, fq1, fq2, refgenome, refgff


Snippy:

* input table: sample, fq1, fq2, refgenome
* output: OUTDIR/snippy/SAMPLE/...{samplegenome = snps.consensus.fa}

Liftoff:

* input table: sample, samplegenome (snps.consensus.fa), refgenome, refgff
* output: OUTDIR/liftoff/SAMPLE/...{samplegff = lifted.gff_polished}

Agat:

* input table: sample, samplegenome (OUTDIR/snippy/SAMPLE/snps.consensus.fa), samplegff (OUTDIR/liftoff/SAMPLE/lifted.gff_polished)
* output: OUTDIR/agat/SAMPLE/...{samplecds = predicted_cds.fa, sampleprots = predicted_proteins.fa}
