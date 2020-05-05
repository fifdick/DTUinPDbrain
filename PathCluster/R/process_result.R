#' Extract cluster titles
#'
#' Gathers titles for the pathway clusters
#'
#' @param alist \code{list()} with one element for each pathway, named with the
#'     corresponding pathway name. Each element on the list is a vector of
#'     genes. Only the names of the elements in the vectors will be used (the
#'     gene ID/names)
#'
#' @param amatrix cluster membership matrix (actually, a \code{data.frame}).
#'     1's mean that the pathways belong to the same cluster, 0's if not. The
#'     diagonal of the matrix is 1 if the pathway will be the title of the
#'     cluster, -1 if it was discarded as the title
#'
#' @param pvalues named numeric vector with the pathway-level p-values. The
#'     names correspond to pathway names
#'
#' @return \code{data.frame()} with the following columns: (1) "items" -
#'     pathway names; (2) "title" - name of the cluster (representative
#'     pathway); (3) "abslog10" - -log10(pvalue) for the pathway
#'
#' @export

process_result <- function(alist, amatrix, pvalues){
    message("Processing results...")
    cluster_titles <- names(which(diag(as.matrix(amatrix)) != -1, arr.ind=TRUE))
    clusters <- lapply(cluster_titles, function(title){
                           rownames(amatrix)[which(amatrix[title]==1)]
                })
    names(clusters) <- cluster_titles
    clusters <- qdapTools::list2df(clusters, col1="items", col2="title")
    clusters$items <- as.character(clusters$items)
    pvalues <- data.frame(items=as.character(names(pvalues)), pvalue=as.numeric(pvalues))
    #clusters <- dplyr::left_join(clusters, pvalues, by="items")
    clusters <- merge(clusters, pvalues, by="items", all.x=TRUE, all.y=FALSE)
    clusters$abslog10 <- abs(log10(clusters$pvalue))
    return(clusters)
}
