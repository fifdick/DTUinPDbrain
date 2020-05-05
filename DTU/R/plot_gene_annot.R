#' plot transcript body
#' 
#' Using biomaRt and GenomeGraphs to plot transcript body 
#' @param gene_id ensembl gene id of the gene to be plotted
#' @param tx_ids string vector of transcript ids (ensembl) specifying which of the transcripts of the gene should be plotted
#' @param gene_name string gene name for the plot tilte
#' @param biotype string specifying transcript biotypes as annotation in the plot. Default is NULL.
#' @return GenomeGraphs plot object
#' @export
plot_gene_annot <- function(gene_id, tx_ids, gene_name, biotype = NULL) {
 mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", host = "http://feb2014.archive.ensembl.org", dataset = "hsapiens_gene_ensembl")
 #Draw the gene and then all transcripts of that gene
 gene <- GenomeGraphs::makeGene(id = gene_id, type = "ensembl_gene_id", biomart = mart)
 transcripts <- lapply(tx_ids, function(id) { 
  GenomeGraphs::makeTranscript(id = id, type = "ensembl_transcript_id", biomart = mart,
			       dp = GenomeGraphs::DisplayPars(plotId = TRUE, cex = 1))
 })
 tracks <- list(GenomeGraphs::makeTitle(gene_name), "Gene" = gene)
 if(!(is.null(biotype))) {
  names(transcripts) <- biotype
 }
 tracks <- c(tracks, transcripts)
 tracks <- c(tracks)
 p <- GenomeGraphs::gdPlot(tracks)
 return(p)
}

