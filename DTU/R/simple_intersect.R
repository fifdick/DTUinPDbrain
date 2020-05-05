#' Intersect
#'
#' Determine intersection of elements in vectors in a list
#' @param res_obj list of dataframes resulting from stageR
#' @param tx bool, if TRUE, intersects transcript ids, else gene ids
#' @param gene_col character string of gene id column name
#' @param tx_col character string of transcrip id column name
#' @return character string holding strings which overlap
#' @export
simple_intersect <- function(res_obj, tx = TRUE, gene_col = "geneID", tx_col = "txID") {
 if (tx == TRUE) {
  txlists <- lapply(res_obj, function(o) {
   o[, tx_col]
  })
 } else {
  txlists <- lapply(res_obj, function(o) {
   o[, gene_col]
  })
 }
 names(txlists) <- names(res_obj)
 intersection <- Reduce(intersect, txlists)
 return(intersection)
}
