#' Network from enriched pathways
#'
#' Creates a netowork from a list of enriched pathways, evaluates the
#' significance of all non-empty intersections between two pathways
#'
#' @param alist \code{list()} with one element for each pathway. Each element
#'     is a binary named vector, where the names are the gene names belonging
#'     to the pathway and the values represent nominal significance (or lack
#'     thereof)
#'
#' @return Returns a (square) \code{data.frame} with the p-values of the
#'     pair-wise intersections using a hypergeometric test
#'
#' @export

intersect_signif_matrix <- function(alist){
    df <- sapply(1:length(alist), function(i){
        one <- alist[[i]] # take one pathway
        sig_values <- sapply(1:length(alist), function(j) {
            int <- intersect(names(one), names(alist[[j]])) # compare to all other pathways and take intersection of genes
            N_i <- length(int) # number of intersecting genes 
            S_i <- sum(one[int]) # total number of genes from the first pathway that intersect and are significant
            S_u <- sum(one) + sum(alist[[j]]) # number of genes that are significant in the first and second pathways
            N_u <- length(one) + length(alist[[j]])
            pvalue <- stats::phyper(S_i, S_u, N_u - S_u, N_i)
            #pvalue <- fisher.test(matrix(c(S_i,S_u-S_i,N_i-S_i,(N_u-S_i)-S_u),2,2),alternative="less")$value
            return(pvalue)
        })
        sig_values
    })
    rownames(df) <- names(alist)
    colnames(df) <- names(alist)
    return(df)
}

