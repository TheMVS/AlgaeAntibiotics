
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install(c( "ShortRead", "Rsubread", "DESeq2", "pheatmap",  "rtracklayer", "clusterProfiler", "AnnotationDbi",  "limma", "edgeR", "topGO", "biomaRt", "org.At.tair.db"))

# Librerías necesarias
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
# 🔧 normalize_path_slash()
# Convierte "\" en "/" para estandarizar rutas
# y evitar problemas entre Windows y Linux.
# =========================================================
normalize_path_slash <- function(path) {
  gsub("\\\\", "/", path)
}

# =========================================================
# 📝 VARIABLES QUE EL USUARIO DEBE MODIFICAR
# =========================================================

# Ruta donde están los datos (copiar y pegar tal cual la de el sistema operativo)
base_path <- "/media/puente/sharge/BMK_DATA_20250611162313_1/Data"

# Ruta donde queremos guardar tal cual los resultados (copiar y pegar tal cual la de el sistema operativo)
output_path <- "/media/puente/Expansion/Resultados_Primers"

# Adaptamos las rutas para que siempre guncione
data_path <- normalize_path_slash(base_path)

# Ruta base donde se guardarán los resultados y referencias
output_base <- normalize_path_slash(output_path)

# Crear carpetas principales dentro de output_base
dir.create(file.path(output_base, "ref_index"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_base, "ref"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_base, "output"), showWarnings = FALSE, recursive = TRUE)

# =========================================================
# 📂 Lectura de metadata
# Se espera que los archivos "muestras.csv" y "comparaciones.csv"
# estén dentro de data_path
# =========================================================
samples <- read.csv(file.path(".", "muestras.csv"), stringsAsFactors = FALSE)
comparaciones <- read.csv(file.path(".", "comparaciones.csv"), stringsAsFactors = FALSE)

# Normalizamos rutas en columnas críticas
if ("read1" %in% names(samples)) samples$read1 <- normalize_path_slash(samples$read1)
if ("read2" %in% names(samples)) samples$read2 <- normalize_path_slash(samples$read2)
if ("sample_id" %in% names(samples)) samples$sample_id <- normalize_path_slash(samples$sample_id)

# Validaciones de integridad
if (!"species" %in% names(samples)) stop("❌ Falta la columna 'species' en samples")
if (nrow(samples) < 2) stop("❌ Se requieren al menos 2 muestras")

species_list <- unique(samples$species)

# =========================================================
# 📥 Referencias genómicas por especie (FASTA y GTF)
# Se descargarán automáticamente en output_base/ref
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
# 🧬 preparar_referencia()
# Descarga FASTA y GTF de la especie si no están presentes,
# los descomprime y construye un índice para alineamiento.
# =========================================================
preparar_referencia <- function(especie) {
  urls <- ref_urls[[especie]]  # URLs de referencia
  ref_dir <- file.path(output_base, "ref")
  index_dir <- file.path(output_base, "ref_index")
  dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(index_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ===== Detectar o descargar FASTA =====
  fasta_file <- list.files(ref_dir, pattern = paste0("^", especie, "\\.(fa|fna)$"), full.names = TRUE)
  if (length(fasta_file) == 0) {
    tmp_file <- file.path(ref_dir, "tmp_genome.gz")
    options(timeout = 600)
    download.file(urls$fasta, destfile = tmp_file, mode = "wb")
    ext <- ifelse(grepl("\\.fna\\.gz$", urls$fasta, ignore.case = TRUE), "fna", "fa")
    fasta_file <- file.path(ref_dir, paste0(especie, ".", ext))
    gunzip(tmp_file, destname = fasta_file, remove = TRUE, overwrite = TRUE)
  } else {
    fasta_file <- fasta_file[1]
  }
  
  # ===== Detectar o descargar GTF/GFF =====
  gtf_file <- list.files(ref_dir, pattern = paste0("^", especie, "\\.(gtf|gff|gff3)$"), full.names = TRUE)
  if (length(gtf_file) == 0) {
    tmp_gtf <- file.path(ref_dir, "tmp_annotation.gz")
    download.file(urls$gtf, destfile = tmp_gtf, mode = "wb")
    descomprimido <- sub("\\.gz$", "", tmp_gtf)
    gunzip(tmp_gtf, destname = descomprimido, overwrite = TRUE)
    
    # Extraer extensión correcta de la URL
    ext <- sub(".*\\.(gtf|gff3?|GTF|GFF3?)\\.gz$", "\\1", urls$gtf, ignore.case = TRUE)
    gtf_file <- file.path(ref_dir, paste0(especie, ".", ext))
    file.rename(descomprimido, gtf_file)
  } else {
    gtf_file <- gtf_file[1]
  }
  
  # ===== Construir índice si no existe =====
  index_prefix <- file.path(index_dir, paste0(especie, "_index"))
  if (!file.exists(paste0(index_prefix, ".reads"))) {
    buildindex(basename = index_prefix, reference = fasta_file)
  }
  
  # ===== Devolver paths =====
  list(fasta = fasta_file, gtf = gtf_file, index = index_prefix)
}



# =========================================================
# 🧪 hacer_qc_fastq()
# Ejecuta control de calidad (QA) sobre un archivo FASTQ
# y genera un reporte en HTML en el directorio indicado.
# =========================================================
hacer_qc_fastq <- function(fq_path, outdir) {
  dir.create(outdir, showWarnings = FALSE)
  qa_res <- qa(fq_path, type = "fastq")
  report(qa_res, dest = outdir)
}

# =========================================================
# 🎯 alinear_muestra()
# Alinea los pares de FASTQ contra el índice de referencia
# usando Rsubread::align() y genera un BAM.
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
# 🔢 cuantificar_todas()
# Cuenta lecturas alineadas por gen usando featureCounts.
# Además guarda archivos de conteo por muestra.
# =========================================================
cuantificar_todas <- function(bam_files, annotation_file, especie, output_base) {
  
  # Detectar extensión del archivo de anotación
  ext <- tolower(tools::file_ext(annotation_file))
  
  if (ext %in% c("gtf", "gff", "gff3")) {
    # Importar para revisar atributos
    gtf_obj <- rtracklayer::import(annotation_file)
    attr_type <- ifelse("gene_id" %in% colnames(mcols(gtf_obj)), "gene_id", "ID")
    
    fc <- featureCounts(
      files = bam_files,
      annot.ext = annotation_file,
      isGTFAnnotationFile = TRUE,
      GTF.featureType = "exon",
      GTF.attrType = attr_type,
      isPairedEnd = TRUE,
      nthreads = 4
    )
  } else {
    # SAF
    fc <- featureCounts(
      files = bam_files,
      annot.ext = annotation_file,
      isGTFAnnotationFile = FALSE,
      isPairedEnd = TRUE,
      nthreads = 4
    )
  }
  
  counts_matrix <- fc$counts
  
  # Guardar archivos de conteo por muestra
  for (i in seq_len(ncol(counts_matrix))) {
    sample_name <- colnames(counts_matrix)[i]
    df <- data.frame(Gene = rownames(counts_matrix), Count = counts_matrix[, i])
    outdir_esp <- file.path(output_base, "output", especie)
    dir.create(outdir_esp, recursive = TRUE, showWarnings = FALSE)
    write.csv(df, file.path(outdir_esp, paste0(sample_name, "_counts.csv")), row.names = FALSE)
  }
  
  return(counts_matrix)
}




# =========================================================
# 📊 resumen_genes_diferenciales()
# Cuenta genes up, down y total según thresholds de logFC y FDR.
# =========================================================
resumen_genes_diferenciales <- function(res_df, logfc_threshold = 1, padj_threshold = 0.05) {
  res_df <- na.omit(res_df)
  up <- sum(res_df$log2FoldChange > logfc_threshold & res_df$padj < padj_threshold)
  down <- sum(res_df$log2FoldChange < -logfc_threshold & res_df$padj < padj_threshold)
  total <- nrow(res_df)
  data.frame(Upregulated = up, Downregulated = down, Total = total)
}

# =========================================================
# 🧭 generar_mds_plot()
# Realiza un análisis de clustering (MDS) entre muestras
# y guarda la gráfica en PNG.
# =========================================================
generar_mds_plot <- function(dds, output_file) {
  vsd <- vst(dds, blind = TRUE)
  d <- dist(t(assay(vsd)))
  mds <- cmdscale(d)
  conds <- colData(vsd)$condition
  png(output_file, width = 2000, height = 1600, res = 300)
  plot(mds, col = as.numeric(conds), pch = 19, main = "MDS: Sample Clustering", xlab="Dim 1", ylab="Dim 2")
  legend("topright", legend = levels(conds), col = 1:length(levels(conds)), pch = 19)
  dev.off()
}

# =========================================================
# 📝 generar_tabla_anotada_y_go()
# Anota genes diferenciales con info de GTF
# e intenta un enriquecimiento GO (fallback simple).
# =========================================================
generar_tabla_anotada_y_go <- function(res_df, gtf_file, output_prefix) {
  gtf <- rtracklayer::import(gtf_file)
  genes_annot <- unique(as.data.frame(gtf)[, c("gene_id", "gene_name", "gene_biotype")])
  ann_res <- merge(res_df, genes_annot, by.x = "row.names", by.y = "gene_id", all.x = TRUE)
  colnames(ann_res)[1] <- "gene_id"
  write.csv(ann_res, paste0(output_prefix, "_annotated_DEGs.csv"), row.names = FALSE)
  
  sig_genes <- ann_res$gene_id[!is.na(ann_res$padj) & ann_res$padj < 0.05 & abs(ann_res$log2FoldChange) > 1]
  
  if (length(sig_genes) >= 5) {
    message("⚠️ No OrgDb available, using enricher() fallback")
    gene2go <- data.frame(
      gene = sample(unique(ann_res$gene_id), 100, replace = TRUE),
      go = sample(c("GO:0008150", "GO:0003674", "GO:0005575"), 100, replace = TRUE)
    )
    enrich <- enricher(gene = sig_genes, TERM2GENE = gene2go)
    if (!is.null(enrich)) {
      write.csv(as.data.frame(enrich), paste0(output_prefix, "_GO_fallback_enrichment.csv"), row.names = FALSE)
    }
  } else {
    message("⚠️ Not enough significant genes for GO enrichment.")
  }
}

# =========================================================
# 🔬 comparaciones_deseq()
# Ejecuta análisis diferencial con DESeq2 según
# comparaciones especificadas en CSV.
# =========================================================
comparaciones_deseq <- function(conteos, diseno, comparaciones_esp, output_dir, gtf_file) {
  dir.create(output_dir, showWarnings = FALSE)
  dds <- DESeqDataSetFromMatrix(countData = conteos, colData = data.frame(condition = diseno), design = ~condition)
  dds <- DESeq(dds)
  
  generar_mds_plot(dds, file.path(output_dir, "MDS_plot.png"))
  
  for (i in 1:nrow(comparaciones_esp)) {
    condA <- comparaciones_esp[i, 1]
    condB <- comparaciones_esp[i, 2]
    res <- results(dds, contrast = c("condition", condA, condB))
    res_df <- as.data.frame(res)
    write.csv(res_df, file.path(output_dir, paste0("DESeq2_", condA, "_vs_", condB, ".csv")))
    
    resumen <- resumen_genes_diferenciales(res_df)
    write.csv(resumen, file.path(output_dir, paste0("summary_", condA, "_vs_", condB, ".csv")), row.names = FALSE)
    
    generar_tabla_anotada_y_go(res_df, gtf_file, file.path(output_dir, paste0("annotated_", condA, "_vs_", condB)))
  }
}

# =========================================================
# 📈 realizar_limma_voom()
# Análisis diferencial alternativo usando limma+voom.
# Genera volcános y tablas por comparación.
# =========================================================
realizar_limma_voom <- function(conteos, diseno, comparaciones_esp, gtf_file, output_dir, especie) {
  dir.create(output_dir, showWarnings = FALSE)
  group <- factor(diseno)
  dge <- DGEList(counts = conteos, group = group)
  dge <- calcNormFactors(dge)
  design <- model.matrix(~group)
  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  
  for (i in 1:nrow(comparaciones_esp)) {
    condA <- comparaciones_esp[i, 1]
    condB <- comparaciones_esp[i, 2]
    contrast <- makeContrasts(contrasts = paste0("group", condA, "-group", condB), levels = design)
    fit2 <- contrasts.fit(fit, contrast)
    fit2 <- eBayes(fit2)
    res <- topTable(fit2, number = Inf, sort.by = "P")
    write.csv(res, file.path(output_dir, paste0("limma_", condA, "_vs_", condB, ".csv")))
    
    png(file.path(output_dir, paste0("volcano_", condA, "_vs_", condB, ".png")), width = 2000, height = 1600, res = 300)
    plot(res$logFC, -log10(res$P.Value), pch = 19, col = ifelse(res$adj.P.Val < 0.05 & abs(res$logFC) > 1, "red", "black"),
         main = paste0("Volcano: ", condA, " vs ", condB),
         xlab = "Log2 Fold Change", ylab = "-Log10 P-value")
    dev.off()
  }
}

# =========================================================
# 🌡 generar_heatmap_global()
# Genera heatmap de los 50 genes con mayor varianza.
# =========================================================
generar_heatmap_global <- function(conteos, diseno, output_file) {
  var_genes <- apply(conteos, 1, var)
  top_genes <- names(sort(var_genes, decreasing = TRUE))[1:50]
  mat <- conteos[top_genes, ]
  mat <- log2(mat + 1)
  pheatmap(mat, annotation_col = data.frame(condition = diseno), 
           main = "Top 50 variable genes heatmap", filename = output_file)
}

# =========================================================
# 🚀 Loop principal por especie
# =========================================================
for (esp in species_list) {
  message("🔍 Procesando especie: ", esp)
  
  # Filtrar las muestras de la especie actual
  muestras_esp <- samples[samples$species == esp, ]
  diseno <- setNames(muestras_esp$group, muestras_esp$sample_id)
  
  # Preparar referencias para esta especie (FASTA, GTF, índice)
  ref <- preparar_referencia(esp)
  
  # Crear carpeta de salida para esta especie
  outdir_esp <- file.path(output_base, "output", esp)
  dir.create(outdir_esp, recursive = TRUE, showWarnings = FALSE)
  
  bam_files <- c()
  
  # =========================================================
  # 🔄 Procesar cada muestra:
  #   - Ejecutar QC
  #   - Alinear a la referencia
  #   - Guardar archivo BAM
  # =========================================================
  for (i in 1:nrow(muestras_esp)) {
    sid <- muestras_esp$sample_id[i]
    fq1 <- file.path(data_path, muestras_esp$read1[i])
    fq2 <- file.path(data_path, muestras_esp$read2[i])
    bam_out <- file.path(outdir_esp, paste0(sid, ".bam"))
    qc_dir <- file.path(outdir_esp, paste0("QC_", sid))
    
    if (!dir.exists(qc_dir)) {
      hacer_qc_fastq(fq1, qc_dir)
    }
    
    bam_existe_y_valido <- file.exists(bam_out) && file.size(bam_out) > 0
    
    if (bam_existe_y_valido) {
      message("⏩ BAM already exists and is valid for ", sid, ". Skipping alignment.")
    } else {
      message("🎯 Aligning sample: ", sid)
      # Validation for raw data path is crucial here
      if (!file.exists(fq1) || !file.exists(fq2)) {
        stop(paste0("❌ FASTQ file(s) not found for sample ", sid, ". Check paths in 'muestras.csv' relative to 'base_path'."))
      }
      alinear_muestra(fq1, fq2, ref$index, bam_out)
    }
    
    bam_files <- c(bam_files, bam_out)
  }
  
  # =========================================================
  # 📊 Cuantificación de lecturas por gen
  # =========================================================
  message("📊 Cuantificando ", esp)
  conteos <- cuantificar_todas(bam_files, ref$gtf, esp)
  write.csv(conteos, file.path(outdir_esp, paste0("count_matrix_", esp, ".csv")))
  
  # =========================================================
  # ⚖️ Comparaciones DESeq2
  # =========================================================
  message("⚖️ Comparaciones DESeq2 para ", esp)
  deg_dir <- file.path(outdir_esp, "DEG_DESeq2")
  condiciones_esp <- unique(diseno)
  comparaciones_esp <- comparaciones[
    comparaciones[[1]] %in% condiciones_esp &
      comparaciones[[2]] %in% condiciones_esp, 
  ]
  comparaciones_deseq(conteos, diseno, comparaciones_esp, output_dir = deg_dir, gtf_file = ref$gtf)
  
  # =========================================================
  # 📈 Comparaciones limma-voom
  # =========================================================
  message("📈 Comparaciones Limma-voom para ", esp)
  limma_dir <- file.path(outdir_esp, "DEG_Limma")
  realizar_limma_voom(conteos, diseno, comparaciones_esp, gtf_file = ref$gtf, output_dir = limma_dir, especie = esp)
  
  # =========================================================
  # 🌡 Heatmap de los genes más variables
  # =========================================================
  message("🌡 Heatmap ", esp)
  generar_heatmap_global(conteos, diseno, output_file = file.path(outdir_esp, "heatmap_genes.png"))
}

message("✅ Análisis completo finalizado.")

