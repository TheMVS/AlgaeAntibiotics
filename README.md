# RNA-seq Analysis Pipeline

This document explains how to use the RNA-seq analysis pipeline implemented in R.

## Prerequisites

### R Packages

Make sure the following packages are installed. You can use BiocManager if needed:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("ShortRead", "Rsubread", "DESeq2", "pheatmap", "rtracklayer",
                       "clusterProfiler", "AnnotationDbi", "limma", "edgeR", "topGO",
                       "biomaRt", "org.At.tair.db"))
```

Also, make sure these packages are loaded in your R session:

```r
library(ShortRead)
library(Rsubread)
library(DESeq2)
library(pheatmap)
library(R.utils)
library(rtracklayer)
library(clusterProfiler)
library(AnnotationDbi)
library(dplyr)
library(edgeR)
library(limma)
library(topGO)
library(biomaRt)
library(org.At.tair.db)
```

## Project Setup

1. Define the paths for your data and outputs:

```r
# Path to raw data
base_path <- "/path/to/your/data"

# Path to store results
output_path <- "."
```

2. The script will automatically create the necessary folders:

```
output_base/ref_index
output_base/ref
output_base/output
```

## Input Files

Place the following CSV files in your data folder:

* `muestras.csv` : Sample information, including columns `sample_id`, `read1`, `read2`, `species`, `group`.
* `comparaciones.csv` : Pairwise comparisons between groups.

## Running the Pipeline

1. Normalize paths and load metadata.
2. Prepare references for each species (FASTA genome and GTF annotation). The script will download them if not present.
3. For each sample:

   * Perform quality control (QC) using `ShortRead::qa()`
   * Align reads with `Rsubread::align()` to generate BAM files
4. Quantify reads using `featureCounts`
5. Perform differential expression analysis using:

   * `DESeq2`
   * `limma-voom` (alternative method)
6. Generate plots:

   * MDS plots for sample clustering
   * Volcano plots for DE genes
   * Heatmaps for top variable genes
7. Annotate DE genes and optionally perform GO enrichment analysis using `clusterProfiler`

## Functions Overview

* `normalize_path_slash(path)`: Standardizes file paths.
* `preparar_referencia(especie)`: Downloads or detects reference genome and annotation, builds index.
* `hacer_qc_fastq(fq_path, outdir)`: Performs QC on FASTQ files.
* `alinear_muestra(fq1, fq2, index_path, output_bam)`: Aligns paired-end reads.
* `cuantificar_todas(bam_files, gtf_path, especie)`: Counts reads per gene.
* `comparaciones_deseq(conteos, diseno, comparaciones_esp, output_dir, gtf_file)`: DESeq2 analysis.
* `realizar_limma_voom(conteos, diseno, comparaciones_esp, gtf_file, output_dir, especie)`: Limma-voom analysis.
* `generar_heatmap_global(conteos, diseno, output_file)`: Heatmap of top variable genes.
* `generar_mds_plot(dds, output_file)`: MDS plot.
* `generar_tabla_anotada_y_go(res_df, gtf_file, output_prefix)`: Annotates DE genes and optional GO enrichment.

## Outputs

All outputs will be stored in `output_base/output/<species>`:

* BAM files for each sample
* Count matrices per species
* DESeq2 and Limma-voom results
* MDS and volcano plots
* Heatmaps of variable genes
* Annotated DE gene tables and GO enrichment results

## Notes

* Make sure you have a stable internet connection for downloading references.
* For large genomes, increase the download timeout with `options(timeout = 600)`.
* The script supports both `.fa` and `.fna` genome files.
* QC reports will be saved as HTML in each sample folder.

