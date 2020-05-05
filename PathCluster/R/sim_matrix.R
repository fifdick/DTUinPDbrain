#' Calculate pair-wise similarity between pathways based on their members
#'
#' From a list of pathways (with gene ID/names as character vectors) outputs a
#' similarity matrix between each pair of pathways
#'
#' @param alist \code{list()} with one element for each pathway, named with the
#'     corresponding pathway name. Each element on the list is a vector of gene
#'     names/IDs ('\code{character}) 
#'
#' @return Returns a (square) \code{data.frame} with the numeric similarity
#'     values (Kappa scores)
#'
#' @export

sim_matrix <- function (alist){
    #require(parallel)
    amatrix <- gene_pathway_matrix(alist)
    #cl <- makeCluster(8)
    df <- sapply(1:ncol(amatrix), function(i){
                     one <- as.data.frame(amatrix[,i])
                     sapply(1:ncol(amatrix), function(j){
                                irr::kappa2(data.frame(one, amatrix[,j]))$value
                     })
          })
    rownames(df) <- colnames(amatrix)
    colnames(df) <- colnames(amatrix)
    return(df)
}
