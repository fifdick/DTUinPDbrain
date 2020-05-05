#' Update similarity matrix
#'
#' Updates similarity matrix keeping a member of a pathway pair and discarding
#' the other
#'
#' @param failed name of the pathway that will be discarded (character)
#'
#' @param passed name of the pathway that will prevail (character)
#'
#' @param simmat similarity matrix between pairs of pathways
#'
#' @param alist named list of pathways with the hard-thresholded significance at the gene-level
#'
#' @return Returns an updated similarity matrix

update_simmat <- function(failed, passed, simmat, alist){
    diag(simmat) <- 1
    gene_pathway_mat <- gene_pathway_matrix(alist)
    genes_in_cluster <- rep(0, length(rownames(gene_pathway_mat)))
    names(genes_in_cluster) <- rownames(gene_pathway_mat)
    merged <- unique(c(names(alist[[failed]]), names(alist[[passed]])))
    #genes_in_cluster <- sapply(names(genes_in_cluster), function(z) ifelse(z %in% merged, 1, 0))
    genes_in_cluster <- ifelse(names(genes_in_cluster) %in% merged, 1, 0)
    names(genes_in_cluster) <- rownames(gene_pathway_mat)
    new_kappa_col <- sapply(1:nrow(simmat), function(x){
            irr::kappa2(data.frame(gene_pathway_mat[,x],genes_in_cluster))$value
    })
    simmat[,passed] <- new_kappa_col
    simmat <- simmat[-which(rownames(simmat)==failed), -which(rownames(simmat)==failed)]
    diag(simmat) <- NA
    return(simmat)
}
