#' Wrapper function to remove pathway redundancy
#'
#' Given an input structure with pathway- and gene-level information, remove
#' the redundant pathways (and plot if desired)
#'
#' @param input list with two elements. The first element is a named list of
#'     data.frames, with the names corresponding to pathway names. Each of the
#'     data.frames contains the genes of the pathway (rows) and p-values, ids,
#'     and log-fold changes for each of them (columns). The second element of
#'     the object is a named list of p-values, where the names correspond to
#'     the pathway names, and the p-values to the significance of the
#'     enrichment of the pathway (required)
#'
#' @param subsetsize number of top pathways that will be processed. Normally,
#'     only the top are interesting anyway, and saves a lot of time (default
#'     all of them)
#'
#' @param threshold similarity threshold (kappa score). Pairs of pathways (or
#'     clusters of pathways) above this threshold will be iteratively merged
#'     until the clusters cannot be further merged. The higher the threshold,
#'     the less the pathways will be clustered (default=0.4)
#'
#' @param lowThreshold minimum pathway size (default=50)
#'
#' @param highThreshold maximum pathway size (default=1000)
#'
#' @param plot whether to plot a figure to a PDF file (default=TRUE)
#'
#' @param outPDF directory for output PDF figure (default="./treemap.pdf")
#'
#' @param plotTitle title for the plot (default="Treemap")
#'
#' @return matrix with the clusters of "aggregated" pathways
#'
#' @export

run_pathway_vis <- function(input, subsetsize=length(input$geneset), threshold=0.4,
                            lowThreshold=50, highThreshold=1000, 
                            plot=TRUE, outPDF="./treemap.pdf", plotTitle="Treemap"){
    # max number of pathways allowed are 1000
    if (subsetsize > 1000){
        warning("Maximum subset size allowed is 1000, setting subsetsize <- 1000")
    	subsetsize <- 1000
    }
    
    # parse input to desired format and subset to contain only subsetsize many pathways
    inputParsed <- parse_input(obj=input, subsetsize=subsetsize)
    
    start_time <- Sys.time()
    
    filtered <- redundant_filter(alist=inputParsed$geneset, threshold=threshold,
                                 pvalues=inputParsed$pvals, lowThreshold=lowThreshold, highThreshold=highThreshold)
    
    end_time <- Sys.time()
    message(paste0("Time used: ", end_time - start_time))
    
    clusters <- process_result(alist=inputParsed$geneset, amatrix=filtered$clusterMat, pvalues=inputParsed$pvals)

    if (plot == TRUE){
        grDevices::pdf(outPDF)
    	print(print_map(clusters, title=plotTitle))
    	grDevices::dev.off()
    }
    
    return(clusters)
}

