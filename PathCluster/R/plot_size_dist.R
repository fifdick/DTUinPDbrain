#' Plot distribution of pathway sizes
#'
#' Plots (using ggplot2) the distribution of pathway sizes (or returns the
#' data.frame that makes it)
#'
#' @param alist \code{list()} with one element for each pathway, named with the
#'     corresponding pathway name. Each element on the list is a vector of gene
#'     names/IDs ('\code{character}) 
#'
#' @param lowThreshold minimum number of genes in a pathway to pass the filter
#'     (default=50)
#'
#' @param highThreshold maximum number of genes in a pathway to pass the filter
#'     (default=1000)
#'
#' @param returnData if TRUE, returns the data.frame instead of the plot
#'     (default=FALSE)
#'
#' @return ggplot with the  distribution of pathway sizes
#'
#' @export

plot_size_dist <- function(alist, lowThreshold=50, highThreshold=1000, returnData=FALSE){
    #require("ggplot2")

    bin_mat <- gene_pathway_matrix(alist) # membership matrix
    gene_count <- colSums(bin_mat) # number of genes per pathway
    filt <- ifelse((gene_count > lowThreshold) & (gene_count < highThreshold), "passed", "failed")

    gene_count <- data.frame(count=gene_count, filtered=filt)
    cols <- c(passed="palegreen2", failed="grey68")
    if (returnData) {
        return(gene_count)
    } else {
        ggplot2::ggplot(data=gene_count, ggplot2::aes(x="count", fill="filtered")) +
            #ggplot2::geom_histogram(ggplot2::aes(y = ..density..)) +
            ggplot2::geom_histogram() +
            ggplot2::labs(title="Pathway size distribution", x="Number of genes per pathway", y="density", fill="Hard filter") +
            ggplot2::geom_vline(xintercept=lowThreshold, col="red") +
            ggplot2::geom_vline(xintercept=highThreshold, col="red") +
            ggplot2::scale_fill_manual(labels=c("failed", paste0(lowThreshold, " < x < ", highThreshold)), values=cols) +
            ggplot2::stat_function(fun=stats::dnorm, args=list(mean=mean(gene_count$count), sd=stats::sd(gene_count$count)), lwd=1, col='red')
    }
}

