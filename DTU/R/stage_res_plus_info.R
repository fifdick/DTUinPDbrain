#' Generate main result dataframe
#' 
#' Makes a df to more easily compare tool and cohort results. Df have same colnames so can be melted etc for plotting
#'  requires stageR result lists, the drimseqdexseq objects and which cohort to use
#' @param list of stageR result dataframes
#' @param list of DRIMSeq and DEXSeq objects
#' @param cohort string cohort name to be used to extract from the lists
#' @return list of dataframes (one for each tool) holding all relevant result information 
#' (nominal p-values, stageR screened p-values, (one for each tool) condition coefficient (adjusted for DRIMSeq (log2(exp(coefficient)))))
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export 
stage_res_plus_info <- function(genelists, ds, cohort) {
# starting with DRIMSeq:
nom_pval <- DRIMSeq::results(ds$Ds_drim[[cohort]], level = "feature") %>%
 dplyr::select(.data$feature_id, .data$pvalue, .data$adj_pvalue) %>%
 dplyr::rename(nom_pval = .data$pvalue, tx_id = .data$feature_id, toolAdjPval = .data$adj_pvalue) 
gene_info <- get_gene_info(id_list = nom_pval$tx_id, tx = TRUE) %>%
 dplyr::select(.data$gene_id, .data$tx_id, .data$gene_name, .data$tx_biotype)
nom_pval <- dplyr::left_join(nom_pval, gene_info, by = "tx_id")
sigdrim <- genelists[[cohort]]$drim_gene_lsts %>%
 dplyr::select(.data$txID, .data$gene, .data$transcript)
colnames(sigdrim) <- c("tx_id", "gene_pvalueStageR", "tx_pvalueStageR")
sigdrim <- dplyr::right_join(sigdrim, nom_pval, by = "tx_id") %>%
 dplyr::arrange(.data$nom_pval)
coeff_drim <- DRIMSeq::coefficients(ds$Ds_drim[[cohort]], level = "feature") %>%
 dplyr::select(.data$feature_id, .data$condition1) %>%
 dplyr::rename(tx_id = .data$feature_id, lnfc = .data$condition1) 
sigdrim <- dplyr::left_join(sigdrim, coeff_drim, by = "tx_id") %>%
 dplyr::mutate(l2fc = log2(exp(.data$lnfc)))
sigdrim <- sigdrim %>% dplyr::select(-.data$lnfc)

# DEXSeq:
nom_pval <- as.data.frame(DEXSeq::DEXSeqResults(ds$Ds_dex[[cohort]], independentFiltering = FALSE))  %>%
 dplyr::select(.data$featureID, .data$pvalue, .data$padj) %>%
 dplyr::rename(nom_pval = .data$pvalue, tx_id = .data$featureID, toolAdjPval = .data$padj) 
gene_info <- get_gene_info(nom_pval$tx_id, tx = TRUE) %>%
 dplyr::select(.data$gene_id, .data$tx_id, .data$gene_name, .data$tx_biotype)
nom_pval <- dplyr::left_join(nom_pval, gene_info, by = "tx_id")
sigdex <- genelists[[cohort]]$dex_gene_lsts %>% dplyr::select(.data$txID, .data$gene, .data$transcript)
colnames(sigdex) <- c("tx_id", "gene_pvalueStageR", "tx_pvalueStageR")
sigdex <- dplyr::right_join(sigdex, nom_pval, by = "tx_id") %>%
 dplyr::arrange(.data$nom_pval)
dex_coeffs <- as.data.frame(SummarizedExperiment::mcols(ds$Ds_dex[[cohort]])) %>%
 dplyr::select(.data$featureID, .data$exonthis.condition1) %>%
 dplyr::rename(tx_id = .data$featureID, l2fc = .data$exonthis.condition1)
sigdex <- dplyr::left_join(sigdex, dex_coeffs, by = "tx_id")

# column shows which rank this event has in the alternative tool result
sigdex <- sigdex %>% dplyr::mutate(rankDrim = match(.data$tx_id, sigdrim$tx_id))
sigdrim <- sigdrim %>% dplyr::mutate(rankDex = match(.data$tx_id, sigdex$tx_id))

return(list("drim" = sigdrim, "dex" = sigdex))
}


