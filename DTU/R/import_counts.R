#' tximport wrapper
#'
#' Import transcript-level counts from Salmon output directories using sampe IDs from the info dataframe from DTU::create_sample_info
#' @param salmon_data base directory of folders with salmon output (one folder per sample), e.g. \"myHome/salmonOutput/\"
#' @param info dataframe or list of dataframes holding sample meta data created from DTU::create_sample_info(). One dataframe per cohort.
#' @param quant_tool character string indicating which tarnscript quantification tool was used. This is redundant. Always used salmon. Default: salmon
#' @param scale_method character string indicating which scaling method should be used by tximport . Default scaledTPM as suggested for DTU analysis.
#' @param annot annotation df holding annotations for the transcripts. Created by DTU::get_annot() using a reference file (used the .gtf file in the referenceData dir)
#' @return list of tximport objects, one for each cohort.
#' @export
import_counts <- function(salmon_data, info, quant_tool="salmon", scale_method="scaledTPM", annot) {
 if (is.list(info)) {
  files_per_cohort <- lapply(info, function(inf) {
   files <- file.path(salmon_data, inf$sample_id, "quant.sf")
   names(files) <- inf$sample_id
   return(files)
  })
  names(files_per_cohort) <- names(info)
  txi_lst <- lapply(files_per_cohort, function(fs) {
   tx2gene <- annot[, c(2, 1)]
   colnames(tx2gene) <- c("tx", "gene")
   txi <- tximport::tximport(fs, type = quant_tool, txOut = TRUE, tx2gene = tx2gene, countsFromAbundance = scale_method)
   cts <- txi$counts
   #filter out non existant transcripts ( 0 for all samples)
   cts <- cts[rowSums(cts) > 0, ]
   txi$counts <- cts
   return(txi)
  })
  names(txi_lst) <- names(info)
  print("DEBUG")
  a <- apply(txi_lst[[1]]$cts,1, function(tx) {cor(tx, info[[1]]$Oligo_Genes)})
  print(summary(a))

  return(txi_lst)
 } else {
  files_per_cohort <- file.path(salmon_data, info$sample_id, "quant.sf")
  names(files_per_cohort) <- info$sample_id
  #import
  tx2gene=annot[, c(2, 1)]
  colnames(tx2gene) <- c("tx", "gene")
  txi <- tximport::tximport(files_per_cohort, type = quant_tool, tx2gene = tx2gene, txOut = TRUE, countsFromAbundance = scale_method)
  cts <- txi$counts
  #filter out non existant transcripts ( 0 for all samples)
  cts <- cts[rowSums(cts) > 0, ]
  txi$counts <- cts
  return(txi)
 }
}


