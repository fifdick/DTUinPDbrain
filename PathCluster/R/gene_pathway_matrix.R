#' Calculate membership matrix from pathway list
#'
#' From a \code{list()} of N pathways (with names corresponding to pathway
#' names/IDs) in which each element is a vector of gene names/IDs, creates an
#' NxM matrix (\code{data.frame}, actually) where N is the number of pathways
#' and M is the "uniqued" number of genes. The values are 1 or 0 designating
#' the presence of the gene in the patwhay
#'
#' @param alist \code{list()} with one element for each pathway, named with the
#'     corresponding pathway name. Each element on the list is a vector of
#'     genes. Only the names of the elements in the vectors will be used (the
#'     gene ID/names)
#'
#' @return Returns a numeric \code{data.frame} with 1's indicating membership
#'     of the gene to the pathway and 0 the lack thereof. Columns are pathways,
#'     rows are genes
#'
#' @export

gene_pathway_matrix <- function(alist){
    gene_names <- unique(unlist(lapply(alist, names)))
    pathway_names <- names(alist)
    
    df <- sapply(1:length(pathway_names), function(i){
        vec <- alist[[i]]
        sapply(1:length(gene_names), function(j){
            ifelse(gene_names[j] %in% names(vec), 1, 0)
        })
    })
    
    colnames(df) <- pathway_names
    rownames(df) <- gene_names
    return(df)
}
