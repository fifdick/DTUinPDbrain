#' Extract counts from tool object
#' 
#' Used for plotting. Extract counts from respective tools
#' @param obj either DRIMSeq or DEXSeq object
#' @param tool <"drim"|"dex"> 
#' @param gene_info df with gene annotation from DTU::get_gene_info()
#' @return dataframe with raw counts
#' @export
#' @importFrom rlang .data
#' @importFrom dplyr %>% 
extract_raw_counts <- function(obj, tool, gene_info) {
 # check which tool it is to know how to get counts
 # then extract the counts on transcript level
 if (grepl("dex", tool, fixed = TRUE)) {
  print("extracting counts from DEXSeq")
  cs <- DEXSeq::DEXSeqResults(obj, independentFiltering = F)$countData
  cs <- as.data.frame(cs)
  cs$transcript_id <- as.character(rownames(cs))
  cs <- tidyr::separate(data = cs, col = transcript_id, into = c("gene_id", "transcript_id"), sep = ":")
  cs <- subset(cs, gene_id %in% gene_info$gene_id)
  count_mats <- cs
 } else if (grepl("drim", tool, fixed = TRUE)) {
  cs <- DRIMSeq::counts(obj)
  colnames(cs)[2] <- "transcript_id"
  cs <- subset(cs, gene_id %in% gene_info$gene_id)
  count_mats <- cs
 }
 # melt to prepare for plotting with ggplot
 count_mats <- reshape2::melt(count_mats, id.vars = c("gene_id", "transcript_id"), variable.name = "sample_id")
 # DEXSeq new version changed colnames...
 if (tool == "dex") {
  colnames(count_mats)[3] <- "sample"
  count_mats$sample_id <- dplyr::left_join(count_mats, as.data.frame(DEXSeq::sampleAnnotation(obj)[, 
  c(1, 2)]), by = "sample")$sample_id
 }
 # replace observedCounts by the proportions (txcount/totalTxOutput of gene)
 count_mats %<>% dplyr::group_by(.data$sample_id, .data$gene_id) %>%
  dplyr::mutate(frac = .data$value / (sum(.data$value))) %>%
  dplyr::select("gene_id", "transcript_id", "sample_id", "frac", "value")
 return(count_mats)
}





