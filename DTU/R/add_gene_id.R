#' add gene id
#'
#' not used?
#' @importFrom rlang .data
add_gene_id <- function(info, obj, counts) {
    annot <- obj$isoformFeatures[, c("isoform_id", "gene_id")]
    counts$gene_id <- sapply(as.character(counts$isoform_id), function(iso_id) {
        gene_id <- subset(annot,isoform_id == iso_id)$gene_id
        return(gene_id)
    })
    
    return(counts)
}
