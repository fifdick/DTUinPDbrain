#' Prepare df to plot fitted proportions
#' 
#' decide which tool to extract fitted values from
#' extract them
#' if there is a DESeq2 obj (dge), extract fitted values from there aswell
#' bind fitted, (dge) and observed 
#' @param ds_obj result object from either DRIMSeq or DEXSeq
#' @param dge_obj DESeq2 result object
#' @param observed_counts dataframe of counts unfitted  (from DTU::extract_raw_counts())
#' @param tool <"drim"|"dex"> specifying the tool from which to extract the fitted values
#' @export
#' @importFrom dplyr %>%
#' @importFrom dplyr case_when
#' @importFrom magrittr %<>%
#' @importFrom rlang .data
extract_fitted_bind_observed <- function(ds_obj, dge_obj = NULL, observed_counts, tool) {
 if (grepl("dex", tool, fixed = TRUE)) {
  print("extracting fitted couns from DEXSeq")
  # DEXSeq
  # why dexseq 98 cols? and rownames of drim proportions?"!
  n_samples <- ncol(DEXSeq::featureCounts(ds_obj))
  fitted.common.scale <- as.data.frame(t(t(SummarizedExperiment::assays(ds_obj)[["mu"]][, 1:n_samples]) / ds_obj$sizeFactor[1:n_samples]))
  colnames(fitted.common.scale) <- DEXSeq::sampleAnnotation(ds_obj)$sample_id
  fitted.common.scale$transcript_id <- sub(".*?:", "", rownames(fitted.common.scale))
  fitted.common.scale$gene_id <- sub("?:.*", "", rownames(fitted.common.scale))
  fitted_props <- reshape2::melt(fitted.common.scale, id.vars = c("transcript_id", "gene_id"), variable.name = "sample_id") %>%
   dplyr::group_by(.data$sample_id, .data$gene_id) %>% dplyr::mutate(frac = .data$value / sum(.data$value))
 } else if (grepl("drim", tool, fixed = TRUE)) { 
  #DRIMSeq
  fitted_props <- DRIMSeq::proportions(ds_obj)
  fitted_props <- reshape2::melt(fitted_props, id.vars = c("gene_id", "feature_id"), variable.name = "sample_id") %>% dplyr::group_by(.data$sample_id, .data$gene_id)
  colnames(fittedProps) <- c("gene_id", "transcript_id", "sample_id", "frac")
 }
 if (!is.null(dge_obj)) {
  dds <- dge_obj
  fitted.common.scale_dge <- t(t(SummarizedExperiment::assays(dds)[["mu"]]) / DESeq2::sizeFactors(dds)) %>%
   reshape2::melt(.)
  colnames(fitted.common.scale_dge) <- c("gene_id", "sample_id", "frac")
  # for the dge dataframe add just the gene_id as transcript_id to be able to bind the rows with the fittedProps
  fitted.common.scale_dge %<>% dplyr::mutate(transcript_id = .data$gene_id) %>% dplyr::select(.data$gene_id, .data$transcript_id, .data$sample_id, .data$frac)
  # resiudals_dge <- counts(dds, normalized=TRUE) - fitted.common.scale
  observed_dge <- t(DESeq2::counts(dds, normalized = TRUE)) %>%
   reshape2::melt(.) 
  colnames(observed_dge) <- c("sample_id", "gene_id", "frac")
  observed_dge %<>% dplyr::mutate(transcript_id = .data$gene_id) %>% dplyr::select(.data$gene_id, .data$transcript_id, .data$sample_id, .data$frac)
  df <- dplyr::bind_rows("fitted_tx" = fitted_props, "fitted_dge" = fitted.common.scale_dge, "observed_tx" = observed_counts, "observed_dge" = observed_dge, .id ="countType")
 } else {
  df <- dplyr::bind_rows("fitted_tx" = fitted_props, "observed_tx" = observed_counts, .id = "countType")
 }
 return(df)
}



