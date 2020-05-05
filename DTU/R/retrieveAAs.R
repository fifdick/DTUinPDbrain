#' Retrieve protein sequence
#'
#' Get protein sequence of transcript isoforms from ensembldb 
#' @param tx_ids string vector with transcript ids identifying transcripts for which protein sequence should be extracted. Must be compatible with ensembldb release version specified in eV.
#' @param eV string specifying ensembldb release version to be used
#' @param out string specifying file to which output should be written.
#' @param append bool specifying whether protein sequences should be appended to content in file or if file should be overwritten. Default is FALSE.
#' @return dataframe with protein sequences
#' @importFrom rlang .data
#' @importFrom dplyr %>%
#' @export 
retrieveAAs <- function(tx_ids, eV="75", out = "./test.fasta", append=F) {
 if (!(paste0("EnsDb.Hsapiens.v", eV) %in% rownames(utils::installed.packages()))) { 
  BiocManager::install(paste0("EnsDb.Hsapiens.v", eV)) 
 }
 library(paste0("EnsDb.Hsapiens.v", eV), character.only = TRUE)
 edb <- get(paste0("EnsDb.Hsapiens.v", eV))
 stopifnot(ensembldb::hasProteinData(edb))
 txs <- ensembldb::transcripts(edb, filter = list(AnnotationFilter::TxIdFilter(tx_ids),
  ensembldb::UniprotMappingTypeFilter("DIRECT")), columns = c("gene_id", "uniprot_id", "tx_biotype", "tx_id", "uniprot_db", "protein_sequence", "gene_name"))
 txs <- as.data.frame(txs) %>% dplyr::filter(!(is.na(.data$uniprot_id)), .data$uniprot_db == "SWISSPROT") %>% dplyr::mutate(nAA = nchar(.data$protein_sequence))
 if (append) {
  seqinr::write.fasta(sequences = as.list(txs$protein_sequence), names = paste0(txs$tx_id, "_", txs$gene_id), as.string = FALSE, file.out = out, open = "a")
 } else {
  seqinr::write.fasta(sequences = as.list(txs$protein_sequence), names = paste0(txs$tx_id, "_", txs$gene_id), as.string = FALSE, file.out = out, open = "w")
 }
 return(txs)
}

