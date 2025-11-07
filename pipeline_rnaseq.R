# =========================================================
# 🧩 Instalación de librerías necesarias (solo 1ª vez)
# =========================================================
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install(c("ShortRead", "Rsubread", "DESeq2", "pheatmap", "rtracklayer", "clusterProfiler", "AnnotationDbi", "edgeR", "biomaRt"), force = TRUE)
#install.packages(c("openxlsx", "dplyr", "R.utils", "ggplot2", "stringr", "matrixStats"))

# =========================================================
# 📚 Librerías
# =========================================================
library(ShortRead)
library(Rsubread)
library(DESeq2)
library(pheatmap)
library(rtracklayer)
library(clusterProfiler)
library(AnnotationDbi)
library(dplyr)
library(edgeR)
library(biomaRt)
library(openxlsx)
library(R.utils)
library(ggplot2)
library(stringr)
library(matrixStats) # Asegurando que rowVars esté disponible

# =========================================================
# 🔬 Parámetros Experimentales (Ajustables por Biólogo)
# =========================================================

# 💻 Parámetros de Rendimiento (CPU/Tiempo)
NUM_THREADS <- 4                  # 🔹 Número de núcleos de CPU a usar para alineamiento y conteo.
TIMEOUT_DOWNLOAD_SEC <- 600       # 🔹 Tiempo máximo (segundos) para la descarga de referencias FASTA/GTF.

# 🧹 Parámetros de Filtrado Inicial (edgeR)
UMBRAL_EXPRESION_MINIMA <- 10     # 🔹 Mínimo de conteos (CPM o valor) que un gen debe alcanzar para ser considerado expresado y retenido.

# 🧪 Parámetros de Significación (DESeq2/Resultados)
ALFA_P_ADJ <- 0.05                # 🔹 Umbral de p-valor ajustado (FDR) para considerar un gen como diferencialmente expresado (DE).
UMBRAL_LOG2FC <- 1                # 🔹 Umbral de log2(Fold Change). Se usa en el Volcano Plot para marcar los genes DE (si es > 1).
UMBRAL_GENES_ENRIQ <- 1           # 🔹 Mínimo de genes DE significativos necesarios para ejecutar el análisis de enriquecimiento KEGG.

# 🧬 Parámetros de Enriquecimiento Funcional (clusterProfiler/KEGG)
KEGG_P_CUTOFF <- 0.05             # 🔹 P-valor de corte para los resultados del análisis KEGG.

# =========================================================
# 🔧 Función para normalizar rutas
# =========================================================
normalize_path_slash <- function(path) gsub("\\\\", "/", path)

# =========================================================
# 📝 Variables de usuario (Rutas)
# =========================================================
# ⚠️ AJUSTAR ESTAS RUTAS ANTES DE EJECUTAR
base_path <- "/media/puente/sharge/BMK_DATA_20250611162313_1/Data"
output_path <- "/media/puente/Expansion/Resultados_Primers_Filtrados"

data_path <- normalize_path_slash(base_path)
output_base <- normalize_path_slash(output_path)

dir.create(file.path(output_base, "ref_index"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_base, "ref"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_base, "output"), showWarnings = FALSE, recursive = TRUE)

# =========================================================
# 📂 Lectura de metadata
# =========================================================
# Se asume que los archivos 'muestras.csv' y 'comparaciones.csv'
# están en el directorio de trabajo donde se ejecuta el script.
samples <- read.csv("muestras.csv", stringsAsFactors = FALSE)
comparaciones <- read.csv("comparaciones.csv", stringsAsFactors = FALSE)

if (!"species" %in% names(samples)) stop("❌ Falta columna 'species' en samples")
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
# 🧬 Preparar referencias
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
    ext <- tolower(tools::file_ext(fasta_file))
  } else {
    tmp_file <- file.path(ref_dir, "tmp_genome.gz")
    options(timeout = TIMEOUT_DOWNLOAD_SEC) # Usando variable global
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
    ext <- tolower(tools::file_ext(gtf_file))
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
    # ⚠️ RECALCULA EL ÍNDICE si falta
    message("⏳ Generando índice con Rsubread::buildindex...")
    buildindex(basename = index_prefix, reference = fasta_file)
    message("✅ Índice generado.")
  }
  
  list(fasta = fasta_file, gtf = gtf_file, index = index_prefix)
}
# =========================================================
# 🧪 QC FASTQ
# =========================================================
hacer_qc_fastq <- function(fq_path, outdir) {
  dir.create(outdir, showWarnings = FALSE)
  qa_res <- qa(fq_path, type = "fastq")
  report(qa_res, dest = outdir)
}

# =========================================================
# 🎯 Alineación
# =========================================================
alinear_muestra <- function(fq1, fq2, index_path, output_bam) {
  if (file.exists(output_bam)) {
    message("⏩ BAM existente: ", output_bam)
    return()
  }
  align(index = index_path, readfile1 = fq1, readfile2 = fq2,
        output_file = output_bam, input_format = "gzFASTQ",
        output_format = "BAM", nthreads = NUM_THREADS) # Usando variable global
}

# =========================================================
# 🔢 Cuantificación
# =========================================================
cuantificar_todas <- function(bam_files, annotation_file, especie, output_base) {
  outdir_esp <- file.path(output_base, "output", especie)
  dir.create(outdir_esp, recursive = TRUE, showWarnings = FALSE)
  
  count_file <- file.path(outdir_esp, "count_matrix.csv")
  if (file.exists(count_file)) {
    message("⏩ Conteos existentes, se cargan.")
    return(read.csv2(count_file, row.names = 1))
  }
  
  # Detectar formato
  ext <- tolower(tools::file_ext(annotation_file))
  
  # Leer primeras líneas para decidir cómo proceder
  # n = 200 es un valor fijo de control para la detección de formato
  gtf_lines <- readLines(annotation_file, n = 200) 
  has_gene_id <- any(grepl("gene_id", gtf_lines))
  has_parent <- any(grepl("Parent=", gtf_lines))
  has_id <- any(grepl("ID=", gtf_lines))
  
  # Elegir configuración según tipo
  if (ext %in% c("gtf")) {
    # NCBI usa gene_id
    attr_type <- "gene_id"
    feature_type <- "exon"
    message("🧬 Archivo GTF detectado, usando attrType='", attr_type, "' y featureType='", feature_type, "'")
    
  } else if (ext %in% c("gff", "gff3")) {
    # Phytozome usa Parent e IDs con exones
    attr_type <- "Parent"
    feature_type <- "exon"
    message("🧬 Archivo GFF3 detectado, usando attrType='", attr_type, "' y featureType='", feature_type, "'")
    
  } else {
    stop("❌ No se reconoce el formato de anotación: ", annotation_file)
  }
  
  # Ejecutar featureCounts
  fc <- featureCounts(
    files = bam_files,
    annot.ext = annotation_file,
    isGTFAnnotationFile = TRUE,
    GTF.featureType = feature_type,
    GTF.attrType = attr_type,
    isPairedEnd = TRUE,
    nthreads = NUM_THREADS # Usando variable global
  )
  
  # Limpiar IDs si vienen con '.exon.#'
  rownames(fc$counts) <- sub("\\.exon\\..*$", "", rownames(fc$counts))
  
  # Guardar resultados
  count_matrix <- fc$counts
  colnames(count_matrix) <- sub(".bam$", "", basename(bam_files))
  count_df <- data.frame(GeneID = rownames(count_matrix), count_matrix, check.names = FALSE)
  write.table(count_df, count_file, sep=";", col.names=NA)
  openxlsx::write.xlsx(count_df, file.path(outdir_esp, "count_matrix.xlsx"))
  
  return(count_matrix)
}

# =========================================================
# 🧹 Filtrado de genes poco expresados (edgeR)
# =========================================================
filtrar_genes_bajos <- function(count_matrix, group, especie, output_base, umbral = UMBRAL_EXPRESION_MINIMA) { # Usando variable global como default
  message("🔎 Filtrando genes con umbral de expresión mínimo = ", umbral)
  
  # Crear objeto DGEList
  y <- DGEList(counts = count_matrix, group = group)
  
  total_genes <- nrow(y$counts)
  
  # Filtrado con umbral configurable
  keep <- filterByExpr(y, group = group, min.count = umbral)
  filtered <- y$counts[keep, ]
  
  genes_filtrados <- nrow(filtered)
  porcentaje <- round((genes_filtrados / total_genes) * 100, 2)
  
  message("✅ Genes retenidos: ", genes_filtrados, " / ", total_genes, " (", porcentaje, "%)")
  
  # Crear resumen
  resumen <- data.frame(
    Especie = especie,
    Genes_totales = total_genes,
    Genes_retenidos = genes_filtrados,
    Porcentaje_retenidos = porcentaje,
    Umbral_usado = umbral
  )
  
  # Guardar resumen
  outdir_esp <- file.path(output_base, "output", especie)
  dir.create(outdir_esp, recursive = TRUE, showWarnings = FALSE)
  
  resumen_csv <- file.path(outdir_esp, paste0("resumen_filtrado_", especie, ".csv"))
  resumen_xlsx <- file.path(outdir_esp, paste0("resumen_filtrado_", especie, ".xlsx"))
  
  write.csv(resumen, resumen_csv, row.names = FALSE)
  openxlsx::write.xlsx(resumen, resumen_xlsx, row.names = FALSE)
  
  return(filtered)
}


# =========================================================
# 🌿 Anotación y enriquecimiento funcional
# =========================================================
anotar_y_enriquecer <- function(res_df, especie, output_dir,
                                umbral_min_genes = UMBRAL_GENES_ENRIQ, # Usando variable global
                                alfa_padj = ALFA_P_ADJ, # Usando variable global
                                kegg_p_cutoff = KEGG_P_CUTOFF) { # Usando variable global
  
  sig <- res_df %>% filter(!is.na(padj), padj < alfa_padj)
  if (nrow(sig) < umbral_min_genes) {
    message("⚠️ Se han encontrado ", nrow(sig)," genes significativos, que es menos de ", umbral_min_genes, "genes significativos requeridos. Se omite enriquecimiento.")
    return(NULL)
  }
  
  # ---------------------------------------------------------
  # CHLAMYDOMONAS: normalización robusta de gene names
  # ---------------------------------------------------------
  if (tolower(especie) == "chlamy") {
    message("🔗 Conectando a Ensembl Plants para Chlamydomonas reinhardtii...")
    plant_mart <- useMart(biomart = "plants_mart", host = "https://plants.ensembl.org")
    dataset <- "creinhardtii_eg_gene"
    org <- "cre"
    mart <- useDataset(dataset = dataset, mart = plant_mart)
    
    # Normalización robusta: (mismo código)
    sig <- sig %>%
      mutate(
        gene_id = as.character(gene_id),
        gene_base = str_extract(gene_id, "Cre[A-Za-z0-9]+\\.g\\d+"),
        note = NA_character_
      )
    
    idx_na <- which(is.na(sig$gene_base))
    if (length(idx_na) > 0) {
      sig$gene_base[idx_na] <- sub("_\\d+(?:\\.\\d+)*\\.v[0-9\\.]+$", "",
                                   sig$gene_id[idx_na], perl = TRUE)
      sig$gene_base[idx_na] <- str_extract(sig$gene_base[idx_na], "Cre[A-Za-z0-9]+\\.g\\d+")
      still_na <- idx_na[is.na(sig$gene_base[idx_na])]
      if (length(still_na) > 0) {
        sig$note[still_na] <- "no_pattern_extracted"
      }
    }
    
    cp_idx <- grep("^CreCp\\.", sig$gene_id)
    if (length(cp_idx) > 0) {
      sig$note[cp_idx] <- ifelse(is.na(sig$note[cp_idx]), "discarded_contig_CreCp", paste(sig$note[cp_idx], "discarded_contig_CreCp", sep=";"))
      sig$gene_base[cp_idx] <- NA_character_
    }
    
    message("🧩 Genes que intentaremos mapear (formato limpio, NA=descartado):")
    for (i in seq_len(nrow(sig))) {
      message(sprintf("   • original: %s  -> base: %s  note: %s", sig$gene_id[i], ifelse(is.na(sig$gene_base[i]), "<NA>", sig$gene_base[i]), ifelse(is.na(sig$note[i]), "-", sig$note[i])))
    }
    
    # Mapear solo los gene_base no NA
    values_to_map <- na.omit(unique(sig$gene_base))
    genes_mapeados <- data.frame()
    if (length(values_to_map) > 0) {
      genes_mapeados <- getBM(
        attributes = c("ensembl_gene_id", "external_gene_name"),
        filters = "external_gene_name",
        values = values_to_map,
        mart = mart
      )
    }
    
    if (nrow(genes_mapeados) == 0) {
      message("⚠️ getBM devolvió 0 filas para los nombres proporcionados (posible mismatch de versiones).")
    } else {
      message("🔍 Resultados de mapeo Ensembl (primeras filas):")
      print(head(genes_mapeados))
    }
    
    # Relacionar y marcar mapped TRUE/FALSE
    sig$mapped_gene <- NA_character_
    if (nrow(genes_mapeados) > 0) {
      sig$mapped_gene <- genes_mapeados$ensembl_gene_id[match(sig$gene_base, genes_mapeados$external_gene_name)]
    }
    sig$mapped <- !is.na(sig$mapped_gene)
    
    # Guardar tabla de mapeo (incluye nota y estado)
    mapping_file <- file.path(output_dir, paste0("gene_mapping_", especie, ".csv"))
    mapping_table <- sig %>% dplyr::select(gene_id, gene_base, mapped_gene, mapped, note)
    write.table(mapping_table, mapping_file, sep = ";", row.names = FALSE, quote = FALSE)
    message("💾 Archivo de mapeo guardado en: ", mapping_file)
    
    # Preparar vector para enrichKEGG: usamos mapped_gene (IDs Ensembl) como antes
    mapped_genes <- unique(na.omit(sig$mapped_gene))
    if (length(mapped_genes) == 0) {
      message("⚠️ No se encontraron genes mapeados a Ensembl. Se omite enriquecimiento KEGG para Chlamy.")
      return(mapping_table)  # devolvemos el mapping para que lo revises
    }
    
    message("✅ Genes mapeados correctamente (count): ", length(mapped_genes))
    for (g in mapped_genes) message("   • ", g)
    
    message("🔬 Ejecutando enriquecimiento KEGG para Chlamy (", org, ")...")
    ekegg <- enrichKEGG(
      gene = mapped_genes,
      organism = org,
      pvalueCutoff = kegg_p_cutoff # Usando variable global
    )
    
    if (is.null(ekegg) || nrow(as.data.frame(ekegg)) == 0) {
      message("⚠️ No se encontraron rutas KEGG significativas para Chlamy.")
      return(list(mapping = mapping_table, kegg = NULL))
    }
    
    ekegg_df <- as.data.frame(ekegg)
    write.table(ekegg_df, file.path(output_dir, paste0("KEGG_", especie, ".csv")), sep = ";", row.names = FALSE)
    openxlsx::write.xlsx(ekegg_df, file.path(output_dir, paste0("KEGG_", especie, ".xlsx")))
    message("✅ Enriquecimiento completado para ", especie)
    return(list(mapping = mapping_table, kegg = ekegg_df))
    
    
    # =========================================================
    # 🌱 LOLIUM PERENNE (RefSeq v6.1)
    # =========================================================
  } else if (tolower(especie) == "lolium") {
    message("🌿 Procesando *Lolium perenne* (RefSeq v6.1)...")
    
    org <- "lper"
    
    # Extraer IDs numéricos de GeneID o LOC
    sig <- sig %>%
      mutate(entrez_id = sub("^LOC", "", gene_id)) %>%
      mutate(entrez_id = sub("^GeneID:", "", entrez_id))
    
    entrez_ids <- na.omit(unique(sig$entrez_id))
    
    message("🧩 GeneIDs detectados para enriquecimiento:")
    for (g in entrez_ids) message("   • ", g)
    
    # Guardar tabla de mapeo
    mapping_file <- file.path(output_dir, paste0("gene_mapping_", especie, ".csv"))
    write.table(
      data.frame(gene_id = sig$gene_id, entrez_id = sig$entrez_id),
      mapping_file, sep = ";", row.names = FALSE
    )
    message("💾 Archivo de mapeo guardado en: ", mapping_file)
    
    # Enriquecimiento KEGG
    message("🔬 Ejecutando enriquecimiento KEGG para Lolium (", org, ")...")
    ekegg <- enrichKEGG(
      gene = entrez_ids,
      organism = org,
      keyType = "ncbi-geneid",
      pvalueCutoff = kegg_p_cutoff # Usando variable global
    )
    
  } else {
    message("⚠️ Especie no reconocida: ", especie)
    return(NULL)
  }
  
  # =========================================================
  # ✅ Procesar resultados y exportar
  # =========================================================
  if (is.null(ekegg) || nrow(as.data.frame(ekegg)) == 0) {
    message("⚠️ No se encontraron rutas KEGG significativas para ", especie)
    return(NULL)
  }
  
  ekegg_df <- as.data.frame(ekegg)
  message("📈 Rutas KEGG significativas encontradas: ", nrow(ekegg_df))
  
  write.table(
    ekegg_df,
    file.path(output_dir, paste0("KEGG_", especie, ".csv")),
    sep = ";", row.names = FALSE
  )
  openxlsx::write.xlsx(
    ekegg_df,
    file.path(output_dir, paste0("KEGG_", especie, ".xlsx"))
  )
  
  message("✅ Enriquecimiento completado para ", especie)
  return(ekegg_df)
}

# =========================================================
# 📊 Volcán y Heatmap
# =========================================================
plot_volcano <- function(res_df, output_file) {
  res_df <- na.omit(res_df)
  # Usando variables globales para los umbrales de color:
  gg <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = padj < ALFA_P_ADJ & abs(log2FoldChange) > UMBRAL_LOG2FC)) +
    theme_minimal() +
    xlab("log2 Fold Change") + ylab("-log10(padj)") +
    scale_color_manual(values = c("grey", "red")) +
    theme(legend.position = "none")
  ggsave(output_file, gg, width = 6, height = 5, dpi = 300)
}

plot_heatmap <- function(dds, output_file) {
  vsd <- vst(dds, blind = TRUE)
  # Asegurando que matrixStats::rowVars esté disponible.
  topVarGenes <- head(order(matrixStats::rowVars(assay(vsd)), decreasing = TRUE), 50)
  mat <- assay(vsd)[topVarGenes, ]
  mat <- mat - rowMeans(mat)
  pheatmap(mat, cluster_rows = TRUE, cluster_cols = TRUE,
           show_rownames = TRUE, filename = output_file, width = 8, height = 8)
}

# =========================================================
# 🔬 DESeq2 + anotación + gráficos
# =========================================================
comparaciones_deseq <- function(conteos, diseno, comparaciones_esp, output_dir, especie) {
  
  # Flag de control: si el archivo existe, significa que ya se ejecutó DESeq2 previamente
  res_flag <- file.path(output_dir, "deseq_done.flag")
  if (file.exists(res_flag))  
    return(message("⏩ DESeq2 ya ejecutado."))
  
  # Crea la tabla de metadatos de las muestras (colData)
  colData <- data.frame(condition = factor(diseno), row.names = names(diseno))
  
  # Crea el objeto DESeqDataSet con la matriz de conteos y las condiciones
  dds <- DESeqDataSetFromMatrix(countData = conteos, colData = colData, design = ~condition)
  
  # Ejecuta todo el pipeline interno de DESeq2:
  dds <- DESeq(dds)
  
  # Genera un heatmap de los 50 genes con mayor varianza entre muestras
  plot_heatmap(dds, file.path(output_dir, "heatmap_top50.png"))
  
  # Recorre cada comparación indicada en la tabla 'comparaciones_esp'
  for (i in 1:nrow(comparaciones_esp)) {
    condA <- comparaciones_esp$control[i]  # condición control
    condB <- comparaciones_esp$treat[i]    # condición tratada
    
    # Extrae los resultados del contraste (condB vs condA)
    res <- results(dds, contrast = c("condition", condA, condB))
    
    # Convierte los resultados en un data.frame para manipularlos fácilmente
    res_df <- as.data.frame(res)
    # Añade una columna con los IDs de los genes
    res_df$gene_id <- rownames(res_df)
    
    # Define nombre de salida para los resultados (CSV y Excel)
    csv_file <- file.path(output_dir, paste0("DESeq2_", condA, "_vs_", condB, ".csv"))
    
    # Guarda los resultados solo si no existen previamente
    if (!file.exists(csv_file)) {
      write.table(res_df, csv_file, sep=";", row.names=FALSE)
      write.xlsx(res_df, sub(".csv", ".xlsx", csv_file), row.names=FALSE)
    }
    
    # Crea el gráfico tipo Volcano plot, usando los umbrales globales
    plot_volcano(res_df, file.path(output_dir, paste0("Volcano_", condA, "_vs_", condB, ".png")))
    
    # Ejecuta la anotación funcional y el análisis de enriquecimiento KEGG
    # Los parámetros ahora se pasan de forma explícita/usando los valores globales por defecto
    anotar_y_enriquecer(res_df, especie, output_dir)
  }
  
  # Marca la ejecución completada creando un archivo “flag”
  file.create(res_flag)
}


# =========================================================
# 🚀 Loop principal
# =========================================================
for (esp in species_list) {
  message("🔍 Procesando especie: ", esp)
  muestras_esp <- samples[samples$species==esp,]
  diseno <- setNames(muestras_esp$group, muestras_esp$sample_id)
  
  ref <- preparar_referencia(esp)
  outdir_esp <- file.path(output_base, "output", esp)
  dir.create(outdir_esp, recursive = TRUE, showWarnings = FALSE)
  
  bam_files <- c()
  for (i in 1:nrow(muestras_esp)) {
    sid <- muestras_esp$sample_id[i]
    fq1 <- file.path(data_path, muestras_esp$read1[i])
    fq2 <- file.path(data_path, muestras_esp$read2[i])
    bam_out <- file.path(outdir_esp, paste0(sid, ".bam"))
    alinear_muestra(fq1, fq2, ref$index, bam_out)
    bam_files <- c(bam_files, bam_out)
  }
  
  counts <- cuantificar_todas(bam_files, ref$gtf, esp, output_base)
  
  # El filtrado usa el UMBRAL_EXPRESION_MINIMA definido al inicio
  counts <- filtrar_genes_bajos(counts, diseno, esp, output_base) 
  
  compar_esp <- comparaciones[(comparaciones$control %in% diseno) & (comparaciones$treat %in% diseno),]
  
  comparaciones_deseq(counts, diseno, compar_esp, outdir_esp, esp)
  
  # Ejecutar Limma-voom
  #comparaciones_limma(count_matrix, diseno, compar_esp, outdir_esp, muestras_esp)
}