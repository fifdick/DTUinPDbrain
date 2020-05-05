#' extract significant hits
#'
#' Extract significant events
#' 
#' Go through list of stageR result dataframes to extract significant events. Not sure where I use this, since stageR has a function for this and results from stageR are only valid with the FDR that was used to run stageR
#' @param obj list of stageR result dataframes, where transcript is the transcript p-value adjusted by stageR and gene is the column for the gene-level adjusted pvalue 
#' @return list of dataframes holding only the significant events
#' @export
extract_sig <- function(obj, alpha=0.05) {
 obj_sig <- lapply(seq(1, length(obj)), function(o) {
  s <- subset(obj[[o]], transcript < alpha & gene < alpha)
  return(s)
 })
 names(obj_sig) <- names(obj)
 return(obj_sig)
}
