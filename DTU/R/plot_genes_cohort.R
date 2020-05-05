#' Extract counts and plot DTU 
#' 
#' extracts fitted and observed counts and plots relative transcript expressions for specific tool and specific cohort
#' @param out string specifying output directory for pdf of plot 
#' @param obj result data object (see makeFigures.Rmd)
#' @param want_jitter bool whether to plot fitted and observed as points jittered (old, not used anymore)
#' @param selected_samples string vector of sample ids, of which to extract the counts for plotting
#' @param genes_toplot df ; if is.null(gene_name) provide the df from DTU::get_gene_info with all transcripts that should be plotted.
#' (Df can have rows with transcript from more than one gene, then a list of plots, one for each gene will be returned)
#' @param gene_name string specifying the gene which should be plotted
#' @param tool <"drim"|"dex">
#' @param plot bool, if TRUE, print plot
#' @param only_nom_sig bool whether to mark sepcific samples in the plot with a txt label
#' @export
#' @importFrom rlang .data
#' @importFrom dplyr %>%
#' @importFrom magrittr %<>%
plot_genes_cohort <- function(out = "./", obj, want_jitter = TRUE, conditions = c("Control", "Case"),
 cohort = NULL, selected_samples = NULL, genes_toplot = NULL, gene_name = NULL, tool = "dex", plot = T, only_nom_sig = T) {   

 counts <- obj$Ds  # accessed by counts[[Ds_tool]][[cohort]]             
 counts <- counts[[paste0("Ds_", tool)]][[cohort]] #for counts 
 # DGE obj from DESeq2
 dge_obj <- obj$dge[[cohort]]
 dge_results <- obj$dge_results[[cohort]]
 main_df <- obj$main_df[[cohort]][[tool]]
 # add dge significance
 main_df <- dplyr::left_join(main_df, dge_results[, c("gene_id", "padj", "pvalue")], by = "gene_id") %>%
  dplyr::rename(padj_dge = .data$padj, pval_dge = .data$pvalue)
 info <- obj$info[[cohort]]
 if (is.null(genes_toplot) && !is.null(gene_name)) {
  genes_toplot <- get_gene_info(subset(main_df, gene_name == gene_name)$tx_id, tx = TRUE)       
 } else {
  if (!is.data.frame(genes_toplot) || !(c("tx_id", "gene_name") %in% colnames(genes_toplot))) {
   stop("Please provide correct gene info dataframe with tx_id and gene_name column without NA")
  }	
 }
 if (tool == "dex" || tool == "drim") {
  print(paste0("extracting counts from the provided count object assuming it is a ", tool, " object"))
  observed_counts <- extract_raw_counts(obj = counts, tool = tool, gene_info = genes_toplot)
  # extract the fitted values from the required tool and bind with the observed counts 
  all_vals <- extract_fitted_bind_observed(ds_obj = counts, dge_obj = dge_obj, observed_counts = observed_counts, tool)
  all_vals %<>% dplyr::filter(.data$gene_id %in% genes_toplot$gene_id)
 } else {
  print("Tool can either be dex or drim, you either didnt specify the function argument or provided a string that is not recognized, will assume you meant dex")
  tool <- "dex"
 }    
 # add significant and condition column using info df and result dataframe of desired tool    
 all_vals <- add_condition_and_sig_tx(info = info, countsmelted = all_vals, main_df = main_df) 
 # plot and generate nice df
 plot_data <- lapply(unique(genes_toplot$gene_id), function(p_gene_id) {
  result <- plot_gene(want_jitter = want_jitter, selected_samples = selected_samples,
		      all_vals = all_vals, p_gene_id = p_gene_id, genes_toplot = genes_toplot, plot = plot,
		      out = out, tool = tool, only_nom_sig = only_nom_sig, conditions = conditions)
  return(result)
 })
 names(plot_data) <- unique(genes_toplot$gene_id)
 return(plot_data) 
}
