#' create annotation df
#' 
#' Create transcript annotation file to map transcripts to genes. Either with ensembldb or read from supplied gtf path. Used gtf file in referenceData dir.
#' @param version string ensembldb version (just the number). Default=75
#' @param g_scource character string genome annotation source, ensembldb or ucsc (->gtf file). Remember there were some modifications needed to be done for the gtf file.
#' @param ntx_filter integer threshold for maximum number of transcripts per gene. Genes with more transcripts will be removed from the annotation 
#' @return annotation dataframe as required by DRIMSeq in DTU::prep_dtu() to build DRIMSeqData object
#' @importFrom rlang .data 
#' @export
get_annot <- function(version = "75", g_source = "ucsc", gtf_path = NULL, ntx_filter = NULL) {
 if (g_source == "ensembldb") {
  db <- paste0("EnsDb.Hsapiens.v", version)
  if (!(db %in% rownames(utils::installed.packages()))) {
   BiocManager::install(db)
  }
  if (!("ensembldb" %in% rownames(utils::installed.packages()))) {
   BiocManager::install("ensembldb")
  }
  library(db, character.only = TRUE)
  # Ensembl mapping between genes and transcripts
  edb <- get(db)
  txdb <- as.data.frame(ensembldb::transcriptsBy(edb, by = "gene",
   columns = c("gene_id", "tx_id", "gene_name", "tx_biotype", "entrezid")))
  txdb <- subset(txdb, ensembldb::seqnames %in% c(1:22, "X", "Y", "MT") & startsWith(gene_id, "ENSG"))
  txdb <- txdb[, c("tx_id", "gene_id", "gene_name", "tx_biotype", "entrezid")]
  names(txdb) <- c("TXNAME", "GENEID", "GENENAME", "BIOTYPE", "ENTREZ")
  annot <- txdb[, c(2, 1)]
  tab <- table(annot$GENEID)
  annot$ntx <- tab[match(annot$GENEID, names(tab))]
  if (!(is.null(ntx_filter))) {
   annot <- subset(annot, ntx <= ntx_filter)
  }
  cat(crayon::blue("ensembl db obj info"))
  print(edb)
  return(annot)
 }
 else if (g_source == "ucsc") {
  if (is.null(gtf_path)) {
   cat(crayon::blue("you chose to load annotation from ucsc gtf file, you didnt provide a path to your gtf file (e.g. ucsc.hg19.gtf)"))
   return(NULL)
  } else {
   txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_path)
   annot <- ensembldb::select(txdb, ensembldb::keys(txdb, "GENEID"), "TXNAME", "GENEID")
   tab <- table(annot$GENEID)
   annot$ntx <- tab[match(annot$GENEID, names(tab))]
   if (!(is.null(ntx_filter))) {
    annot <- subset(annot, ntx <= ntx_filter)
   }
   return(annot)
  }
 } else {
  cat(crayon::blue("The database you want to use does not exist or not offered to be used within this function"))
  cat(crayon::blue("Will return NULL"))
  return(NULL)
 }
}
