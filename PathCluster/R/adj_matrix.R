#' Calculate a pair-wise adjacency matrix for a list of pathways
#'
#' From a list of pathways (in which each element is a named vector of genes
#'     with binary values if the genes were differentially expressed or not),
#'     calculates a binary pairwise adjacency matrix based on the significance
#'     of a hypergeometric test of the intersections
#'
#' @param alist \code{list()} with one element for each pathway. Each element
#'     is a binary named vector, where the names are the names/IDs of the genes
#'     that belong to the pathway and the values represent nominal significance
#'     (or lack thereof)
#'
#' @param alpha hard threshold to transform the p-values of the test to a
#'     binary adjacency matrix
#'
#' @return Returns a (square) \code{data.frame} with the p-values of the
#'     pair-wise intersections using a hypergeometric test
#'
#' @export

adj_matrix <- function(alist, alpha=0.05){
    sig_mat <- intersect_signif_matrix(alist) # Get p-values of pair-wise intersections
    adj_matrix <- apply(sig_mat, c(1,2), function(x) ifelse(x < alpha, 1, 0))
    rownames(adj_matrix) <- rownames(sig_mat)
    colnames(adj_matrix) <- colnames(sig_mat)
    return(adj_matrix)
}

