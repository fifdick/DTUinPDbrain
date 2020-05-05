#' Identify isoform switch
#' 
#' Defined a gene to have an isoform switch, if the top expressed transcript differs between groups. Requires a minimum of 2 DTU events per gene.
#' @param df dataframe resulting from DTU::plot_genes_cohort(). (p[[1]]$data)
#' @param gene string gene identifier (Gene symbol)
#' @return list with bool specifying whether or not there is an isoform switch, and if yes, the transcript id and biotype of the top (relative) expressed transcript of the control and case group separately
#' @export
#' @importFrom dplyr %>%
#' @importFrom magrittr %<>%
#' @importFrom rlang .data
# TODO improve this documentation, it sucks.
isISwitch <- function(df, gene) {
 control_biotype <- ""
 case_biotype <- ""	
 df_m <- df %>% dplyr::group_by(.data$gene_id, .data$transcript_id, .data$condition) %>%
  dplyr::distinct(.data$transcript_id, .keep_all=T)
 df_m %<>% dplyr::group_by(.data$condition) %>%
  dplyr::summarize(sumNomSig = sum(.data$nomSig))
 # if there are at least two nominally significant transcript in this gene, and both conditions have the same number of nom sig transcripts (is this necessary?!)
 if ((df_m[1, 2] == df_m[2, 2]) && (df_m[1, 2] >= 2)) {      
  isoform_switch  <-  FALSE 
  # -->commented following line to see what changes 
  #	df <- dplyr::filter(df,nomSig==1)
  df %<>% dplyr::group_by(.data$condition, .data$transcript_id) %>%
   dplyr::mutate(median = stats::median.default(.data$frac))
  # interested in ties aswell
  # distinct to select on random sample row representing the transcript. The median in this transcript and condition group will be the same in all sample rows
  # group by condition (again) and run ranks on it. Return all entries with max rank (can be more than one)
  df %<>% dplyr::distinct(.data$transcript_id, .keep_all = T) %>%
   dplyr::group_by(.data$condition) %>%
   dplyr::slice(which(rank(.data$median, ties.method = "min") == max(rank(median, ties.method = "min"))))
  # add requirement that at least one transcript in each condition that was assigned max rank has to be nomSig==1 could be higher than 1 if there are ties
  main_expr_sig <- df %>% dplyr::group_by(.data$condition) %>% dplyr::summarize(sumNomSig = sum(.data$nomSig))
  if (subset(main_expr_sig, condition == 0)$sumNomSig >= 1 && subset(main_expr_sig, condition == 1)$sumNomSig >= 1) {
   if (nrow(df) > 2 && ((nrow(df) %% 2) != 0)) {		
    print(paste0("In Gene ", gene, " in one of the conditions, there is more than one dominant transcript"))
    print(df)
   } else if (nrow(df) > 2 && ((nrow(df)%%2) == 0)) {
    print(paste0("In Gene ", gene, ", both conditions show more than one dominant transcript"))
    print(df)
   } else if (nrow(df) == 1) {
    print(paste0("In Gene ", gene, ": for one of the conditions the dominant transcript cant be determined"))
    print(df)
   } else if (nrow(df) == 0) {
    print(paste0("In Gene ", gene, " in neither condition a dominant transcript was determined"))
    print(df)
   } else if (nrow(df) == 2) {
    if (df[1, "transcript_id"] == df[2, "transcript_id"]) {
     isoform_switch <- FALSE
     print(paste0("Both conditions share the same dominant transcript - NO isoform switch detected!"))
    } else {
     print(df)			
     isoform_switch <- TRUE
     print(paste0("Conditions show different dominant transcripts - isoform switch detected!"))
     control_biotype <- gsub(".*\n", "", dplyr::filter(df, .data$condition == 0)$tx_id_label)
     case_biotype <- gsub(".*\n", "", dplyr::filter(df, .data$condition == 1)$tx_id_label) 
     print(paste0("Switch from CT: ", control_biotype, " to CS: ", case_biotype, " transcript"))
     # TODO: fix this
     print(df[,c(1, 2, 6, 7, 8, 12)])	
    }
   }
  } else {
   print(paste0("Gene ", gene, " doesnt show the main expressed transcript in the DTU event"))
  }
 } else {
  print(paste0("Gene ", gene, " doesnt have at least 2 nominally significant transcripts"))
  isoform_switch <- FALSE
 }
 return(list(isIS = isoform_switch, CT.tx_id = dplyr::filter(df, .data$condition == 0)$transcript_id, CT.tx_biotype = control_biotype, CS.tx_id = dplyr::filter(df, .data$condition == 1)$transcript_id, CS.tx_biotype = case_biotype))
}

