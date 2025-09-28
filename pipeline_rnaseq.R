#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install(c("ShortRead", "Rsubread", "DESeq2", "pheatmap", "rtracklayer",
#                       "clusterProfiler", "AnnotationDbi", "limma", "edgeR", "topGO", 
#                       "biomaRt", "org.At.tair.db"))

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

# =========================================================
# 🔧 Función para normalizar rutas
# =========================================================
normalize_path_slash <- function(path) {
  gsub("\\\\", "/", path)
}

# =========================================================
# 📝 Variables de usuario
# =========================================================
base_path <- "/media/puente/sharge/BMK_DATA_20250611162313_1/Data"
output_path <- "/media/puente/Expansion/Resultados_Primers"
data_path <- normalize_path_slash(base_path)
output_base <- normalize_path_slash(output_path)

dir.create(file.path(output_base, "ref_index"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_base, "ref"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_base, "output"), showWarnings = FALSE, recursive = TRUE)

# =========================================================
# 📂 Lectura de metadata
# =========================================================
samples <- read.csv(file.path(".", "muestras.csv"), stringsAsFactors = FALSE)
comparaciones <- read.csv(file.path(".", "comparaciones.csv"), stringsAsFactors = FALSE)

# Normalizar rutas
if ("read1" %in% names(samples)) samples$read1 <- normalize_path_slash(samples$read1)
if ("read2" %in% names(samples)) samples$read2 <- normalize_path_slash(samples$read2)
if ("sample_id" %in% names(samples)) samples$sample_id <- normalize_path_slash(samples$sample_id)

# Validaciones
if (!"species" %in% names(samples)) stop("❌ Falta la columna 'species' en samples")
if (nrow(samples) < 2) stop("❌ Se requieren al menos 2 muestras")
species_list <- unique(samples$species)

# =========================================================
# 📥 URLs de referencia
# =========================================================
ref_urls <- list(
  chlamy = list(
    fasta = "ftp://ftp.ensemblgenomes.org/pub/plants/release-61/fasta/chlamydomonas_reinhardtii/dna/Chlamydomonas_reinhardtii.Chlamydomonas_reinhardtii_v5.5.dna.toplevel.fa.gz",
    gtf   = "ftp://ftp.ensemblgenomes.org/pub/plants/release-61/gtf/chlamydomonas_reinhardtii/Chlamydomonas_reinhardtii.Chlamydomonas_reinhardtii_v5.5.61.gtf.gz"
  ),
  lolium = list(
    fasta = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/359/855/GCF_019359855.2_Kyuss_2.0/GCF_019359855.2_Kyuss_2.0_genomic.fna.gz",
    gtf   = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/359/855/GCF_019359855.2_Kyuss_2.0/GCF_019359855.2_Kyuss_2.0_genomic.gtf.gz"
  )
)

# =========================================================
# 🧬 Función para preparar referencia
# =========================================================
preparar_referencia <- function(especie) {
  urls <- ref_urls[[especie]]
  ref_dir <- file.path(output_base, "ref")
  index_dir <- file.path(output_base, "ref_index")
  dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(index_dir, recursive = TRUE, showWarnings = FALSE)
  
  # FASTA
  fasta_file <- list.files(ref_dir, pattern = paste0("^", especie, "\\.(fa|fna)$"), full.names = TRUE)
  if (length(fasta_file) > 0) {
    message("⏩ FASTA ya existe para ", especie)
    fasta_file <- fasta_file[1]
  } else {
    tmp_file <- file.path(ref_dir, "tmp_genome.gz")
    options(timeout = 600)
    download.file(urls$fasta, destfile = tmp_file, mode = "wb")
    ext <- ifelse(grepl("\\.fna\\.gz$", urls$fasta, ignore.case = TRUE), "fna", "fa")
    fasta_file <- file.path(ref_dir, paste0(especie, ".", ext))
    gunzip(tmp_file, destname = fasta_file, remove = TRUE, overwrite = TRUE)
  }
  
  # GTF/GFF
  gtf_file <- list.files(ref_dir, pattern = paste0("^", especie, "\\.(gtf|gff|gff3)$"), full.names = TRUE)
  if (length(gtf_file) > 0) {
    message("⏩ GTF/GFF ya existe para ", especie)
    gtf_file <- gtf_file[1]
  } else {
    tmp_gtf <- file.path(ref_dir, "tmp_annotation.gz")
    download.file(urls$gtf, destfile = tmp_gtf, mode = "wb")
    descomprimido <- sub("\\.gz$", "", tmp_gtf)
    gunzip(tmp_gtf, destname = descomprimido, overwrite = TRUE)
    ext <- sub(".*\\.(gtf|gff3?|GTF|GFF3?)\\.gz$", "\\1", urls$gtf, ignore.case = TRUE)
    gtf_file <- file.path(ref_dir, paste0(especie, ".", ext))
    file.rename(descomprimido, gtf_file)
  }
  
  # Índice
  index_prefix <- file.path(index_dir, paste0(especie, "_index"))
  if (file.exists(paste0(index_prefix, ".reads"))) {
    message("⏩ Índice ya existe para ", especie)
  } else {
    buildindex(basename = index_prefix, reference = fasta_file)
  }
  
  list(fasta = fasta_file, gtf = gtf_file, index = index_prefix)
}

# =========================================================
# 🧪 Función QC FASTQ
# =========================================================
hacer_qc_fastq <- function(fq_path, outdir) {
  dir.create(outdir, showWarnings = FALSE)
  qa_res <- qa(fq_path, type = "fastq")
  report(qa_res, dest = outdir)
}

# =========================================================
# 🎯 Función alinear
# =========================================================
alinear_muestra <- function(fq1, fq2, index_path, output_bam) {
  fq1 <- normalize_path_slash(fq1)
  fq2 <- normalize_path_slash(fq2)
  output_bam <- normalize_path_slash(output_bam)
  if (!file.exists(fq1)) stop("❌ FASTQ1 not found: ", fq1)
  if (!file.exists(fq2)) stop("❌ FASTQ2 not found: ", fq2)
  align(index = index_path, readfile1 = fq1, readfile2 = fq2,
        output_file = output_bam, input_format = "gzFASTQ",
        output_format = "BAM", nthreads = 4)
}

# =========================================================
# 🔢 Función cuantificación
# =========================================================
cuantificar_todas <- function(bam_files, annotation_file, especie, output_base) {
  outdir_esp <- file.path(output_base, "output", especie)
  dir.create(outdir_esp, recursive = TRUE, showWarnings = FALSE)
  count_flag <- file.path(outdir_esp, "counts_done.flag")
  
  # Si ya existen los conteos, los leemos
  if (file.exists(count_flag)) {
    message("⏩ Conteos ya calculados para ", especie)
    # Leer todos los CSV de conteo y reconstruir la matriz
    count_files <- list.files(outdir_esp, pattern = "_counts\\.csv$", full.names = TRUE)
    counts_list <- lapply(count_files, function(f) {
      df <- read.csv(f, stringsAsFactors = FALSE)
      setNames(df$Count, df$Gene)
    })
    count_matrix <- do.call(cbind, counts_list)
    colnames(count_matrix) <- gsub("_counts\\.csv$", "", basename(count_files))
    return(count_matrix)
  }
  
  # Si no existen los conteos, calcular
  ext <- tolower(tools::file_ext(annotation_file))
  if (ext %in% c("gtf", "gff", "gff3")) {
    gtf_obj <- rtracklayer::import(annotation_file)
    attr_type <- ifelse("gene_id" %in% colnames(mcols(gtf_obj)), "gene_id", "ID")
    fc <- featureCounts(files = bam_files, annot.ext = annotation_file,
                        isGTFAnnotationFile = TRUE, GTF.featureType = "exon",
                        GTF.attrType = attr_type, isPairedEnd = TRUE, nthreads = 4)
  } else {
    fc <- featureCounts(files = bam_files, annot.ext = annotation_file,
                        isGTFAnnotationFile = FALSE, isPairedEnd = TRUE, nthreads = 4)
  }
  
  count_matrix <- fc$counts
  colnames(count_matrix) <- sub(".bam$", "", basename(bam_files))
  for (i in seq_len(ncol(count_matrix))) {
    sample_name <- colnames(count_matrix)[i]
    df <- data.frame(Gene = rownames(count_matrix), Count = count_matrix[, i])
    write.csv(df, file.path(outdir_esp, paste0(sample_name, "_counts.csv")), row.names = FALSE)
  }
  
  # Crear flag para indicar que los conteos ya se hicieron
  file.create(count_flag)
  return(count_matrix)
}

# =========================================================
# 📊 Resumen DEGs
# =========================================================
resumen_genes_diferenciales <- function(res_df, logfc_threshold = 1, padj_threshold = 0.05) {
  res_df <- na.omit(res_df)
  up <- sum(res_df$log2FoldChange > logfc_threshold & res_df$padj < padj_threshold)
  down <- sum(res_df$log2FoldChange < -logfc_threshold & res_df$padj < padj_threshold)
  total <- nrow(res_df)
  data.frame(Upregulated = up, Downregulated = down, Total = total)
}

# =========================================================
# 🧭 MDS Plot
# =========================================================
generar_mds_plot <- function(dds, output_file) {
  if (file.exists(output_file)) {
    message("⏩ MDS plot ya existe, se omite.")
    return()
  }
  vsd <- vst(dds, blind = TRUE)
  d <- dist(t(assay(vsd)))
  mds <- cmdscale(d)
  conds <- colData(vsd)$condition
  png(output_file, width = 2000, height = 1600, res = 300)
  plot(mds, col = as.numeric(factor(conds)), pch = 19,
       main = "MDS: Sample Clustering", xlab="Dim 1", ylab="Dim 2")
  legend("topright", legend = levels(factor(conds)), col = 1:length(levels(factor(conds))), pch = 19)
  dev.off()
}

# =========================================================
# 🔬 DESeq2 con Heatmap
# =========================================================
comparaciones_deseq <- function(conteos, diseno, comparaciones_esp, output_dir, gtf_file, muestras_esp) {
  dir.create(output_dir, showWarnings = FALSE)
  flag_file <- file.path(output_dir, "deseq_done.flag")
  if (file.exists(flag_file)) {
    message("⏩ DESeq2 ya ejecutado previamente en ", output_dir)
    return(invisible(NULL))
  }
  
  colnames(conteos) <- muestras_esp$sample_id
  colData <- data.frame(condition = factor(diseno), row.names = names(diseno))
  dds <- DESeqDataSetFromMatrix(countData = conteos, colData = colData, design = ~condition)
  dds <- DESeq(dds)
  
  # MDS
  generar_mds_plot(dds, file.path(output_dir, "MDS_plot.png"))
  
  # Heatmap
  heatmap_file <- file.path(output_dir, "heatmap_genes_mas_variables.png")
  if (!file.exists(heatmap_file)) {
    vsd <- vst(dds, blind = TRUE)
    topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
    mat <- assay(vsd)[topVarGenes, ]
    mat <- mat - rowMeans(mat)
    
    # ✅ FIX: asegurar que annotation_col tenga rownames que coincidan con columnas de mat
    annotation_col <- data.frame(condition = colData$condition)
    rownames(annotation_col) <- rownames(colData)
    
    png(heatmap_file, width = 2000, height = 2000, res = 300)
    pheatmap(mat,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             show_rownames = TRUE,
             annotation_col = annotation_col)
    dev.off()
  }
  
  # Loop comparaciones
  for (i in 1:nrow(comparaciones_esp)) {
    condA <- comparaciones_esp$control[i]
    condB <- comparaciones_esp$treat[i]
    
    if (!(condA %in% colData$condition) || !(condB %in% colData$condition)) {
      message("⚠️ Saltando comparación DESeq2 inválida: ", condA, " vs ", condB)
      next
    }
    
    csv_file <- file.path(output_dir, paste0("DESeq2_", condA, "_vs_", condB, ".csv"))
    if (file.exists(csv_file)) {
      message("⏩ Resultados DESeq2 ", condA, " vs ", condB, " ya existen.")
      next
    }
    
    res <- results(dds, contrast = c("condition", condA, condB))
    res_df <- as.data.frame(res)
    write.csv(res_df, csv_file)
    
    resumen <- resumen_genes_diferenciales(res_df)
    write.csv(resumen, file.path(output_dir, paste0("summary_", condA, "_vs_", condB, ".csv")), row.names = FALSE)
  }
  
  file.create(flag_file)
}


# =========================================================
# 💻 Limma-voom
# =========================================================
comparaciones_limma <- function(counts, diseno, comparaciones_esp, output_dir, muestras_esp) {
  flag_file <- file.path(output_dir, "limma_done.flag")
  if (file.exists(flag_file)) {
    message("⏩ Limma ya ejecutado previamente en ", output_dir)
    return(invisible(NULL))
  }
  
  colnames(counts) <- muestras_esp$sample_id
  design <- model.matrix(~0 + factor(diseno))
  colnames(design) <- levels(factor(diseno))
  
  y <- DGEList(counts = counts)
  y <- calcNormFactors(y)
  v <- voom(y, design, plot = FALSE)
  
  for (i in 1:nrow(comparaciones_esp)) {
    condA <- comparaciones_esp$control[i]
    condB <- comparaciones_esp$treat[i]
    
    if (!(condA %in% colnames(design)) || !(condB %in% colnames(design))) {
      message("⚠️ Saltando comparación Limma inválida: ", condA, " vs ", condB)
      next
    }
    
    contrast_matrix <- makeContrasts(contrasts = paste0(condB, "-", condA), levels = design)
    fit <- lmFit(v, design)
    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2)
    
    csv_file <- file.path(output_dir, paste0("limma_", condA, "_vs_", condB, ".csv"))
    if (!file.exists(csv_file)) {
      res_df <- topTable(fit2, number = Inf, adjust.method = "BH")
      write.csv(res_df, csv_file)
    }
  }
  
  file.create(flag_file)
}

# =========================================================
# 🚀 Loop principal
# =========================================================
for (esp in species_list) {
  message("🔍 Procesando especie: ", esp)
  
  muestras_esp <- samples[samples$species==esp,]
  diseno <- setNames(muestras_esp$group, muestras_esp$sample_id)
  
  # Preparar referencias
  ref <- preparar_referencia(esp)
  
  outdir_esp <- file.path(output_base, "output", esp)
  dir.create(outdir_esp, recursive = TRUE, showWarnings = FALSE)
  
  bam_files <- c()
  for (i in 1:nrow(muestras_esp)) {
    sid <- muestras_esp$sample_id[i]
    fq1 <- file.path(data_path, muestras_esp$read1[i])
    fq2 <- file.path(data_path, muestras_esp$read2[i])
    bam_out <- file.path(outdir_esp, paste0(sid, ".bam"))
    
    if (!file.exists(bam_out)) {
      alinear_muestra(fq1, fq2, ref$index, bam_out)
    }
    bam_files <- c(bam_files, bam_out)
  }
  
  # Cuantificación
  count_matrix <- cuantificar_todas(bam_files, ref$gtf, esp, output_base)
  
  # Subset de comparaciones válidas para esta especie
  compar_esp <- comparaciones[(comparaciones$control %in% diseno) & (comparaciones$treat %in% diseno),]
  
  # Ejecutar DESeq2
  comparaciones_deseq(count_matrix, diseno, compar_esp, outdir_esp, ref$gtf, muestras_esp)
  
  # Ejecutar Limma-voom
  comparaciones_limma(count_matrix, diseno, compar_esp, outdir_esp, muestras_esp)
}

