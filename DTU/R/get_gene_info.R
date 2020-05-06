#' Retrieve gene annotation from ensembldb
#' 
#' Especially used for making the paper figures. To retrieve transcript biotypes and gene names based on either transcript id or gene id.
#' @param id_list character vector of ids used to retrieve gene annotation, if tx==T then it must be ensembl transcript ids that are valid in the ensembl release that is specified with the parameter version. Else it must be ensembl gene ids.
#' @param version string specifying ensembl release version (default is 75 which was used within this analysis, the last of grch37)
#' @param g_source char string specifing the database used to retrieve annotation (this is now redundant, removed functionality to use ucsc)
#' @param tx bool specifying whether the vector provided with id_list are transcript (tx==TRUE) or gene ids
#' @param exons
#' @param uniprot
#' @param symbol
#' @param promoters
#' @param upBP
#' @param downBPi
#' @return dataframe of requested gene annotation.
#' @importFrom dplyr %>%
#' @export
get_gene_info <- function(id_list, version = "75", g_source = "ensembldb", tx = FALSE, exons=F, uniprot = FALSE, symbol = FALSE, promoters = FALSE, upBP=2000, downBP=400) {
 if (g_source == "ensembldb") {
  db <- paste0("EnsDb.Hsapiens.v", version)   
  if (!(db %in% rownames(utils::installed.packages()))) {
   BiocManager::install(db)
  }
  if (!("ensembldb" %in% rownames(utils::installed.packages()))) {
   BiocManager::install("ensembldb") 
  }
  #not sure how to solve the porblem of having a library call in here
  library(db, character.only = TRUE)
  # Ensembl mapping between genes and transcripts
  edb <- get(db)
  if (symbol == TRUE) {
   print("Using GeneNameFilter")
   genes <- as.data.frame(ensembldb::genes(edb, filter = AnnotationFilter::GeneNameFilter(id_list), columns = c("gene_id", "gene_name", "gene_biotype", "entrezid")))
   return(genes)
  }
  if (tx == FALSE) {
   genes_gr <- ensembldb::genes(edb, filter = AnnotationFilter::GeneIdFilter(id_list), columns = c("gene_id", "gene_name", "gene_biotype", "entrezid"))
   genes <- as.data.frame(genes_gr)
   genes <- subset(genes, .data$seqnames %in% c(1:22, "X", "Y", "MT") & startsWith(gene_id, "ENSG"))
   if (promoters == TRUE) {
    return(promoters(genes_gr, upstream = upBP, downstream = downBP))	
   } else {
    return(genes)
   }
  } else {
   if (exons == FALSE) {
    if (uniprot == TRUE) {
     txdb <- as.data.frame(ensembldb::transcriptsBy(edb, by = "gene", filter = AnnotationFilter::TxIdFilter(id_list), columns = c("gene_id", "tx_id", "gene_name", "tx_biotype", "entrezid", "tx_name", "protein_id", "uniprot_id")))
     txdb <- txdb %>% dplyr::group_by(.data$tx_id) %>% dplyr::mutate(uniprotID = as.character(paste(.data$uniprot_id, collapse = ","))) %>% dplyr::select(-c(.data$uniprot_id, .data$entrezid)) %>% dplyr::distinct()
    } else {
     txdb <- as.data.frame(ensembldb::transcriptsBy(edb, by = "gene", filter = AnnotationFilter::TxIdFilter(id_list), columns = c("gene_id", "tx_id", "gene_name", "tx_biotype", "entrezid", "tx_name")))
    }
    txdb <- subset(txdb, seqnames %in% c(1:22, "X", "Y", "MT") & startsWith(gene_id, "ENSG"))
    tab <- table(txdb$gene_id)
    txdb$ntx <- tab[match(txdb$gene_id, names(tab))]
   } else {
    if (!(tx == TRUE)) {
     print("please enter tx_ids and put tx==T to get exons as they are retrieved from ensembldb by tx_id")
    } else {
     txdb <- as.data.frame(ensembldb::exonsBy(edb, by = "tx", filter = AnnotationFilter::TxIdFilter(id_list), columns = c("gene_id", "tx_id", "gene_name", "tx_biotype", "exon_id", "exon_idx", "exon_seq_start", "exon_seq_end", "seq_strand")))
     txdb <- subset(txdb, seqnames %in% c(1:22, "X", "Y", "MT") & startsWith(gene_id, "ENSG"))
     txdb <- txdb %>% dplyr::group_by(.data$tx_id)
    }    
   }
   return(txdb)
  }
 }
}
