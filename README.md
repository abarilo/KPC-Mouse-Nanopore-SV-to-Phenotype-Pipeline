# KPC Mouse Nanopore SV-to-Phenotype Pipeline

A streamlined long-read pipeline for the KPC mouse PDAC model: somatic and mosaic SVs are called with Sniffles2, high- and moderate- impact events are annotated via Ensembl VEP, and reveal enrichments in nucleotide biosynthesis, mitochondrial respiration, PI3K-AKT-mTOR signaling and oxidative-stress pathways, mapping - a coherent SV-driven cancer phenotype.

This Nextflow, Singularity/Docker pipeline performs:

1.  **QC** of long-read FASTQ (NanoPlot)\
2.  **Alignment** to mouse GRCm39 (minimap2 → sorted/indexed BAM)\
3.  **Structural variant calling** (Sniffles2 mosaic + cohort modes)\
4.  **Annotation** of SVs (Ensembl VEP “tab” output)\
5.  **Linking SVs to tumor phenotypes** via gene lists and public resources

------------------------------------------------------------------------

## 📁 Repository layout

```         
.
├── main.nf                 # Nextflow DSL2 pipeline
├── nextflow.config         # Config profiles + parameters
├── results/                # Nanoplot results, Sniffles2 graphs; Annotation                                   #outputs; top enrichment results
├── Case_study_overview.pdf # Pdf overview of the study
├── Enrichment.RMD          # Pathway enrichment analyis in R
├── Enrichment.html         # Knitted pathway enrichment analyis in R
├── qc-map.Dockerfile       # Dockerfile for QC and allignment
├── variant-call.Dockerfile # Dockerfile for variant calling and annotation
└── README.md               # This file
```

------------------------------------------------------------------------

## 🔧 Before running

1.  **Download data**

    ``` bash
    mkdir case
    cd case
    cd data

    \# Reference genome 
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz

    \# Nanopore reads 
    #Install SRA toolkit if needed 
    prefetch ERR4351540 --max-size 100G 
    prefetch ERR4351539 --max-size 100G 
    fasterq-dump ERR4351540.sra 
    fasterq-dump ERR4351539.sra 

    cd -
    ```

2.  **Download the Nextflow files**

    ``` bash
    wget https://raw.githubusercontent.com/abarilo/KPC-Mouse-Nanopore-SV-to-Phenotype-Pipeline/refs/heads/master/nextflow.config
    wget https://raw.githubusercontent.com/abarilo/KPC-Mouse-Nanopore-SV-to-Phenotype-Pipeline/refs/heads/master/main.nf
    ```

3.  **Pull Singularity/Docker images**

  Singularity:

    ``` bash
    #Singularity
    singularity pull qc-map.sif   docker://abarilo/qc-map:latest
    singularity pull variant-call.sif docker://abarilo/variant-call:latest
    singularity pull ensembl-vep_latest.sif docker://ensemblorg/ensembl-vep:latest
    
    ```
  Docker:
  
    ``` bash
    docker pull abarilo/qc-map:latest
    docker pull abarilo/variant-call:latest
    docker pull ensemblorg/ensembl-vep:latest

    ```
4.  **Install VEP cache**

  Singularity:

    ``` bash
    export VEP_CACHE_DIR=vep_cache
    mkdir -p $VEP_CACHE_DIR
    singularity exec \
    --bind $VEP\_CACHE\_DIR:/root/.vep \\
    ensembl-vep\_latest.sif \\
    perl /opt/vep/src/ensembl-vep/INSTALL.pl \\
     --AUTO cf \\
     --SPECIES mus\_musculus \\
     --ASSEMBLY GRCm39 \\
     --NO\_HTSLIB \\
     --NO\_TEST

    ```
  Docker:

    ``` bash
 
    export VEP_CACHE_DIR=$PWD/vep_cache
    mkdir -p $VEP_CACHE_DIR
    docker run --rm \
    -u root \
    -v $VEP_CACHE_DIR:/root/.vep \
    --entrypoint perl \
    ensemblorg/ensembl-vep:latest \
    /opt/vep/src/ensembl-vep/INSTALL.pl \
    --AUTO cf \
    --SPECIES mus_musculus \
    --ASSEMBLY GRCm39 \
    --NO_HTSLIB \
    --NO_TEST

    ```
------------------------------------------------------------------------

## ⚙️ Running the pipeline

**Local (standard)**

``` bash
nextflow run main.nf \
  --profile standard \
  --fastqs       'pth to *.fastq files' \
  --ref          'path to reference' \
  --outdir       ./results \
  --qc_image     abarilo/qc-map:latest \
  --sv_image     abarilo/variant-call:latest \
  --vepImage     ensemblorg/ensembl-vep:latest \
  --genomeBuild  GRCm39 \
  --vepSpecies   mus_musculus \
  --vepCache     'path to .vep' \
```

**HPC (SLURM)**

``` bash
nextflow run main.nf \
  --profile hpc \
  --fastqs       'pth to *.fastq files' \
  --ref          'path to reference' \
  --outdir       ./results \
  --qc_image     path to qc-map.sif   \
  --sv_image     path to variant-call.sif \
  --vepImage     path to ensembl-vep_latest.sif \
  --genomeBuild  GRCm39 \
  --vepSpecies   mus_musculus \
  --vepCache     'path to .vep' \
```

------------------------------------------------------------------------

## 🔍 Outputs

-   `results/nanoplot/` - QC plots & stats\
-   `results/bams/` - sorted & indexed BAMs\
-   `results/sv/mosaic/` - per-sample mosaic VCFs\
-   `results/sv/population/` - merged cohort VCF\
-   `results/sv/vep/mosaic/` - VEP annotation (VCF + TSV) for mosaic SVs\
-   `results/sv/vep/population/` - VEP annotation of merged cohort\
-   `results/svplots/` - SV plots

------------------------------------------------------------------------

## 📖 Post-processing

Use Enrichment.RMD script for pathway enrichment analysis

------------------------------------------------------------------------

## 🔍 Results folder overview

Case_study_overview.pdf - study overview and brief results

results/nanoplot - QC results

results/sniffles2_plots - Sniffles2 results plotted

results/Vep_results - annotation tables

results/Enrichment_results/top10_go.csv - Top GO enriched terms

results/Enrichment_results/top10_kegg.csv - Top KEGG enriched terms

results/Enrichment_results/top10_react.csv - Top REACTOME enriched terms

results/Enrichment_results/go_df.csv - all GO enriched terms

results/Enrichment_results/kegg_df.csv - all KEGG enriched terms

results/Enrichment_results/react_df.csv - all REACTOME enriched terms
