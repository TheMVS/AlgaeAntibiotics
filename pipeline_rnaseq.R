############################################################
### RNA-seq pipeline FINAL y ROBUSTO
### Chlamydomonas reinhardtii + Lolium perenne
############################################################

options(stringsAsFactors = FALSE)
options(scipen = 999)

suppressPackageStartupMessages({
  library(ShortRead)
  library(Rsubread)
  library(DESeq2)
  library(edgeR)
  library(clusterProfiler)
  library(biomaRt)
  library(dplyr)
  library(openxlsx)
  library(R.utils)
})

############################################################
## PARÁMETROS GLOBALES
############################################################
NUM_THREADS <- 4
MIN_COUNTS  <- 1
ALFA_P_ADJ  <- 0.05

############################################################
## RUTAS
############################################################
base_path   <- "/media/puente/sharge/BMK_DATA_20250611162313_1/Data"
output_path <- "/media/puente/Expansion/Resultados_Raw"

dir.create(output_path, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(output_path,"ref"),       showWarnings=FALSE)
dir.create(file.path(output_path,"ref_index"), showWarnings=FALSE)
dir.create(file.path(output_path,"output"),    showWarnings=FALSE)

############################################################
## METADATA
############################################################
samples       <- read.csv("muestras_raw.csv")
comparaciones <- read.csv("comparaciones.csv")
species_list  <- unique(samples$species)

############################################################
## REFERENCIAS
############################################################
ref_urls <- list(
  chlamy = list(
    fasta = "ftp://ftp.ensemblgenomes.org/pub/plants/release-61/fasta/chlamydomonas_reinhardtii/dna/Chlamydomonas_reinhardtii.Chlamydomonas_reinhardtii_v5.5.dna.toplevel.fa.gz",
    gtf   = "ftp://ftp.ensemblgenomes.org/pub/plants/release-61/gtf/chlamydomonas_reinhardtii/Chlamydomonas_reinhardtii.Chlamydomonas_reinhardtii_v5.5.61.gtf.gz",
    kegg  = "cre"
  ),
  lolium = list(
    fasta = "ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/lolium_perenne/dna/Lolium_perenne.MPB_Lper_Kyuss_1697.dna.toplevel.fa.gz",
    gtf   = "ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/gtf/lolium_perenne/Lolium_perenne.MPB_Lper_Kyuss_1697.62.gtf.gz",
    tsv   = "ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/tsv/lolium_perenne/Lolium_perenne.MPB_Lper_Kyuss_1697.62.uniprot.tsv.gz",
    kegg  = "lper"
  )
)

############################################################
## FUNCIONES
############################################################

preparar_referencia <- function(esp){
  message("🔧 [1/8] Preparando referencia (", esp, ")")
  fa  <- file.path(output_path,"ref",paste0(esp,".fa"))
  gtf <- file.path(output_path,"ref",paste0(esp,".gtf"))
  idx <- file.path(output_path,"ref_index",paste0(esp,"_index"))
  
  if(!file.exists(fa)){
    message("  ⬇️ Descargando FASTA...")
    download.file(ref_urls[[esp]]$fasta, paste0(fa,".gz"), mode="wb")
    gunzip(paste0(fa,".gz"))
  }
  if(!file.exists(gtf)){
    message("  ⬇️ Descargando GTF...")
    download.file(ref_urls[[esp]]$gtf, paste0(gtf,".gz"), mode="wb")
    gunzip(paste0(gtf,".gz"))
  }
  if(!file.exists(paste0(idx,".reads"))){
    message("  🔨 Creando índice Rsubread...")
    buildindex(basename = idx, reference = fa)
  } else {
    message("  ⏩ Índice ya existente")
  }
  list(fasta=fa, gtf=gtf, index=idx)
}

hacer_qc <- function(fq,outdir){
  sumfile <- file.path(outdir,"qaSummary.txt")
  if(file.exists(sumfile)) return()
  dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
  x <- qa(fq, type="fastq")
  report(x, dest=outdir)
  ## pequeño marcador para no repetir
  cat("QC done\n", file=sumfile)
}

alinear <- function(fq1,fq2,idx,bam){
  if(file.exists(bam)){
    message("  ⏩ BAM existente: ", basename(bam))
    return()
  }
  message("  ▶️ Alineando: ", basename(bam))
  align(
    index       = idx,
    readfile1   = fq1,
    readfile2   = fq2,
    output_file = bam,
    type        = "rna",
    nthreads    = NUM_THREADS
  )
}

contar <- function(bams,gtf,outdir){
  fcsv <- file.path(outdir,"counts.csv")
  fxls <- file.path(outdir,"counts.xlsx")
  
  if(file.exists(fcsv)){
    message("  ⏩ Conteos existentes")
    return(as.matrix(read.csv(fcsv,row.names=1,check.names=FALSE)))
  }
  
  message("🔢 [4/8] Conteo de lecturas (featureCounts)")
  fc <- featureCounts(
    files                = bams,
    annot.ext            = gtf,
    isGTFAnnotationFile  = TRUE,
    GTF.featureType      = "exon",
    GTF.attrType         = "gene_id",
    isPairedEnd          = TRUE,
    countMultiMappingReads = TRUE,
    nthreads             = NUM_THREADS
  )
  
  counts <- fc$counts
  counts <- counts[rowSums(counts) > 0,, drop=FALSE]
  
  write.csv(counts, fcsv, quote = FALSE)
  write.xlsx(as.data.frame(counts), fxls)
  message("  ✔ Guardado counts.csv y counts.xlsx")
  counts
}

filtrar <- function(counts,group){
  message("🧹 [5/8] Filtrado edgeR")
  y <- DGEList(counts=counts, group=group)
  keep <- filterByExpr(y, group=group, min.count=MIN_COUNTS)
  if(sum(keep) == 0)
    stop("❌ Ningún gen pasa el filtrado (filterByExpr)")
  y$counts[keep,,drop=FALSE]
}

run_go <- function(sig_genes, esp, outdir, comp){
  message("🧬 [7/8] GO: ", comp)
  
  if(length(sig_genes) == 0){
    message("  ⚠️ Sin genes para GO (lista vacía)")
    write.xlsx(data.frame(), file.path(outdir,paste0("GO_",comp,".xlsx")))
    return()
  }
  
  go_rds <- file.path(output_path,"ref",paste0(esp,"_go_map.rds"))
  
  if(file.exists(go_rds)){
    go_map <- readRDS(go_rds)
  } else {
    mart <- useDataset(
      ifelse(esp=="chlamy","creinhardtii_eg_gene","lperenne_eg_gene"),
      useMart("plants_mart",host="https://plants.ensembl.org")
    )
    go_map <- getBM(
      c("ensembl_gene_id","go_id","name_1006"),
      mart=mart
    )
    saveRDS(go_map, go_rds)
  }
  
  go_map <- go_map[go_map$go_id != "" & !is.na(go_map$go_id), , drop=FALSE]
  if(nrow(go_map) == 0){
    message("  ⚠️ Mapa GO vacío")
    write.xlsx(data.frame(), file.path(outdir,paste0("GO_",comp,".xlsx")))
    return()
  }
  
  term2gene <- go_map |>
    dplyr::select(ID = go_id, gene = ensembl_gene_id)
  
  term2name <- go_map |>
    dplyr::filter(!is.na(name_1006), name_1006 != "") |>
    dplyr::distinct(go_id, .keep_all = TRUE) |>
    dplyr::select(ID = go_id, Description = name_1006)
  
  ego <- tryCatch(
    enricher(
      sig_genes,
      TERM2GENE   = term2gene,
      TERM2NAME   = term2name,
      pvalueCutoff = 1,
      qvalueCutoff = 1
    ),
    error = function(e){
      message("  ⚠️ Error en enricher GO: ", e$message)
      NULL
    }
  )
  
  if(is.null(ego)){
    write.xlsx(data.frame(), file.path(outdir,paste0("GO_",comp,".xlsx")))
    return()
  }
  
  df <- as.data.frame(ego)
  if(nrow(df) > 0){
    df$Significant_pvalue  <- df$pvalue  < ALFA_P_ADJ
    df$Significant_padjust <- df$p.adjust < ALFA_P_ADJ
  }
  
  write.xlsx(df, file.path(outdir,paste0("GO_",comp,".xlsx")))
}

run_kegg <- function(sig_genes, esp, outdir, comp){
  message("🧬 [8/8] KEGG: ", comp)
  
  if(length(sig_genes) == 0){
    message("  ⚠️ Sin genes para KEGG (lista vacía)")
    write.xlsx(data.frame(), file.path(outdir,paste0("KEGG_",comp,".xlsx")))
    return()
  }
  
  ## Para Lolium descargamos el TSV (aunque KEGG puede seguir sin mapear)
  if(esp == "lolium"){
    tsv <- file.path(output_path,"ref","lolium_uniprot.tsv")
    if(!file.exists(tsv)){
      message("  ⬇️ Descargando mapping UniProt Lolium...")
      download.file(ref_urls$lolium$tsv, paste0(tsv,".gz"), mode="wb")
      gunzip(paste0(tsv,".gz"))
    }
    ## Aquí podrías en el futuro hacer un mapeo extra Ensembl -> UniProt -> KEGG
  }
  
  ek <- tryCatch(
    enrichKEGG(
      gene          = sig_genes,
      organism      = ref_urls[[esp]]$kegg,
      pvalueCutoff  = 1,
      qvalueCutoff  = 1
    ),
    error = function(e){
      message("  ⚠️ Error en enrichKEGG: ", e$message)
      NULL
    }
  )
  
  if(is.null(ek)){
    write.xlsx(data.frame(), file.path(outdir,paste0("KEGG_",comp,".xlsx")))
    return()
  }
  
  df <- as.data.frame(ek)
  if(nrow(df) > 0){
    df$Significant_pvalue  <- df$pvalue  < ALFA_P_ADJ
    df$Significant_padjust <- df$p.adjust < ALFA_P_ADJ
  }
  
  write.xlsx(df, file.path(outdir,paste0("KEGG_",comp,".xlsx")))
}

run_deseq <- function(counts, group, comps, outdir, esp){
  message("📊 [6/8] DESeq2")
  
  colData <- data.frame(
    condition = factor(group),
    row.names = colnames(counts)
  )
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = colData,
    design    = ~ condition
  )
  dds <- DESeq(dds)
  
  for(i in 1:nrow(comps)){
    cA   <- comps$control[i]
    cB   <- comps$treat[i]
    comp <- paste0(cB,"_vs_",cA)
    message("  🔎 Comparación: ", comp)
    
    res <- results(dds, contrast = c("condition", cB, cA))
    df  <- as.data.frame(res)
    df$gene_id <- rownames(df)
    
    ## Guardar DESeq2 completo (todos los genes, sin filtrar)
    out_de <- file.path(outdir, paste0("DESeq2_",comp,".xlsx"))
    write.xlsx(df, out_de)
    
    ## Para GO y KEGG: TODOS los genes con padj no NA (como pediste)
    sig_ids <- df$gene_id[!is.na(df$padj)]
    
    run_go(sig_ids,  esp, outdir, comp)
    run_kegg(sig_ids, esp, outdir, comp)
  }
}

############################################################
## LOOP PRINCIPAL
############################################################
for(esp in species_list){
  
  message("\n==================== ", toupper(esp)," ====================")
  
  outdir <- file.path(output_path,"output",esp)
  dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
  
  ref <- preparar_referencia(esp)
  smp <- samples[samples$species==esp,]
  group <- setNames(smp$group, smp$sample_id)
  
  ## QC
  qc_flag <- file.path(outdir,"qc_done.flag")
  if(!file.exists(qc_flag)){
    message("🧪 [2/8] QC FASTQ")
    for(i in 1:nrow(smp)){
      hacer_qc(file.path(base_path,smp$read1[i]),
               file.path(outdir,"QC",paste0(smp$sample_id[i],"_R1")))
      hacer_qc(file.path(base_path,smp$read2[i]),
               file.path(outdir,"QC",paste0(smp$sample_id[i],"_R2")))
    }
    file.create(qc_flag)
  } else {
    message("⏩ QC ya realizado")
  }
  
  ## Alineación
  message("🧬 [3/8] Alineación")
  bams <- character()
  for(i in 1:nrow(smp)){
    bam <- file.path(outdir, paste0(smp$sample_id[i],".bam"))
    alinear(
      fq1 = file.path(base_path, smp$read1[i]),
      fq2 = file.path(base_path, smp$read2[i]),
      idx = ref$index,
      bam = bam
    )
    bams <- c(bams, bam)
  }
  
  ## Conteo
  counts   <- contar(bams, ref$gtf, outdir)
  counts_f <- filtrar(counts, group)
  
  ## Comparaciones válidas (usa los niveles presentes en 'group')
  comps <- comparaciones[
    comparaciones$control %in% group &
      comparaciones$treat   %in% group,
  ]
  
  if(nrow(comps) == 0){
    message("⚠️ No hay comparaciones válidas para ", esp)
  } else {
    run_deseq(counts_f, group, comps, outdir, esp)
  }
}
