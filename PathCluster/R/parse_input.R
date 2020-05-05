#' Parse the main input to the run_pathway-vis function
#'
#' Parses the input object to reduce the total number of pathway to analyse and
#' to discretise significance levels to binary (either significant or not) to
#' simplify downstream analysis
#'
#' @param obj list with two elements. The first element is a named list of
#'     data.frames, with the names corresponding to pathway names. Each
#'     of the data.frames contains the genes of the pathway (rows) and
#'     p-values, ids, and log-fold changes for each of them (columns). The
#'     second element of the object is a named list of p-values, where the
#'     names correspond to the pathway names, and the p-values to the
#'     significance of the enrichment of the pathway
#'
#' @param subsetsize number of top pathways that will be processed. Normally,
#'     only the top are interesting anyway, and saves a lot of time (default
#'     all of them)
#'
#' @param alpha level of significance to threshold gene expression differences
#'     (default 0.05)
#'
#' @param colPval column name in the data.frames from the second slot of the
#'     \code{obj} that contains the p-values (column name or position)
#'
#' @param colGeneID column name in the data.frames from the second slot of the
#'     \code{obj} that contains the gene ID/names (column name of position)
#'
#' @return Returns a \code{list} with each entry being a named vector. The name
#'     of the element in the list corresponds to the pathway name and the
#'     values of the vector are the gene ID/names
#'
#' @export

parse_input <- function(obj, subsetsize=length(obj[[1]]), alpha=0.05, colPval="pvalue", colGeneID="gene_id"){

    stopifnot(is.list(obj))
    stopifnot(length(obj)==2)
    stopifnot(is.list(obj[[1]]))
    stopifnot(is.list(obj[[2]]))
    stopifnot(!(is.null(names(obj[[1]]))))
    stopifnot(!(is.null(names(obj[[2]]))))
    stopifnot(all.equal(names(obj[[1]]), names(obj[[2]])))
    stopifnot(sum(sapply(obj[[2]], is.data.frame)) == length(obj[[2]]))
    stopifnot(all(sapply(obj[[2]], function(x){all(c(colPval, colGeneID) %in% colnames(x))})))

    subsetsize <- min(subsetsize, length(obj[[1]]))

    pathways <- obj[[1]]
    geneDFs <- obj[[2]]
    
    # sort by increasing pvalue
    pathways <- pathways[order(unlist(pathways))]
    geneDFs <- geneDFs[names(pathways)]
    
    input <- lapply(1:subsetsize, function(i){
        vec <- ifelse(geneDFs[[i]][[colPval]] >= alpha, 0, 1)
        names(vec) <- geneDFs[[i]][[colGeneID]]
        return(vec)
    })
    names(input) <- names(geneDFs)[1:subsetsize]

    pvalues <- pathways[1:subsetsize]
    names(pvalues) <- names(pathways)[1:subsetsize]
    return(list("geneset"=input,"pvals"=pvalues))
}


#parse_input_go <- function(df,subsetsize)
#{
#  library(sets)
#  library(stringr)
#  input <- lapply(seq(1,subsetsize),function(y){
#    names <- df$genesets[[y]]$gene_name
#    vec <- ifelse(df$genesets[[y]]$pvalue>0.05,0,1)
#    names(vec) <- names
#    return(vec)
#  })
#  names(input) <- names(df$genesets)[1:subsetsize]
#  names(input) <- sapply(seq(1:length(names(input))),function(x){
#    gsub("^GO:[0-9]*","",names(input)[x])
#  })
#  pvalues <- df$p.values[1:subsetsize]
#  names(pvalues) <- names(input)
#  return(list("geneset"=input,"pvals"=pvalues))
#}

#parse_input_kegg <- function(df)
#{
#  library(sets)
#  library(stringr)
#  input <- lapply(seq(1,length(df$genesets)),function(y){
#    names <- df$genesets[[y]]$gene_name
#    vec <- ifelse(df$genesets[[y]]$pvalue>0.05,0,1)
#    names(vec) <- names
#    return(vec)
#  })
#  names(input) <- names(df$genesets)
#  names(input) <- sapply(seq(1:length(names(input))),function(x){
#    gsub("^hsa[0-9]*","",names(input)[x])
#  })
#  pvalues <- df$p.values
#  names(pvalues) <- names(input)
#  return(list("geneset"=input,"pvals"=pvalues))
#}

