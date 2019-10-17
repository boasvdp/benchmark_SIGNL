# Benchmark SIGNL
This repo contains the Snakemake pipeline used for the benchmarking exercise 2019 of the special interest group for bioinformaticians in medical microbiology in the Netherlands (SIGBMMNL? usually just SIGNL). 

## Installation and dependencies
The pipeline is meant to run on Linux and has three main dependencies:

1. Conda
2. Snakemake
3. A downloaded Kraken2 database

### Conda
(Mini)conda can be installed through https://docs.conda.io/en/latest/miniconda.html#linux-installers. 

### Snakemake
If Miniconda is installed, install snakemake using:

```
conda install -c bioconda -c conda-forge snakemake
```

### Kraken2 database
These can be downloaded from https://ccb.jhu.edu/software/kraken2/downloads.shtml. We used the MiniKraken2_v1_8GB database (no human genome included). You can also compile your own database.

## Pipeline
We used snakemake version 5.7.1 to run the pipeline.

The pipeline entails several steps, depicted in `rulegraph.svg`. In short: 
- Genome assembly using SKESA
- SNP alignment using SKA
- Snp-dists is used to get SNP counts
- QC using Kraken, fastp, Quast, and multiqc 
- AMR genes identification using ABRicate 
- Clustering/phylogenetics using iqtree (from ska alignment) and poppunk
- MLST mainly for backwards compatibility and as extra check

The results from these analyses will be summarised in a final `report.html`. This report contains a csv containing basic stats per strain, a phylogeny, multiqc and fastp reports, and the PopPUNK output. I create this report because it is relatively easy to share results like this (one file), but I will probably edit this in the future if this pipeline will be used more often. It is suboptimal but will do for now.
