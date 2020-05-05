#' prepare data for DTU analysis
#'
#' Used in DTU::run_dex and DTU::run_drim.
#' @param param_obj opt parse object parsed from the command line (runDTU.sh) 
#' @return tx_lst is a list of tximport results (one for each cohort)
#' annot is a list of annotation (one for each cohort)
#' Ds is an object for each cohort that contains the filtered DRIMSeq obj, the pre-filter DRIMSeq obj and the filter parameters used
#' @export

prep_dtu <- function(param_obj) {
 #load pheno
 info <- create_sample_info(pheno = param_obj$m, cohort_names = split_colon(param_obj$n), conditions = split_colon(param_obj$c), batch_vars = as.list(split_colon(param_obj$b)))
 #get annot
 annot_general <- get_annot(version = param_obj$e, gtf_path = param_obj$g, ntx_filter = param_obj$x)
  #import counts
 tx_lst <- import_counts(salmon_data = param_obj$s, info = info, scale_method = param_obj$z, annot = annot_general)
 #check libsizes (the colsums of the count data in each cohort are equal to the number od mapped paired end reads per experiment (in million)
 print("range of library sizes")
 lapply(tx_lst, function(tx) { print(range(colSums(tx$counts)) / 1e6) })
 #check if all transcripts appear in annot
 print("Number of transcripts whose ids are NOT represented in the annotation df")
 lapply(tx_lst, function(tx) { print(sum(!(rownames(tx$counts) %in% annot_general$TXNAME))) })
 #exclude all transcript ids from annotation that dont appear in data and check if both dfs are equal
 #check if cohorts should use different annot dfs. Length of annot after filtering out transcript ids that dont appear in respective count data
 annot <- lapply(tx_lst, function(tx) {	
  a <- annot_general[match(rownames(tx$counts), annot_general$TXNAME), ]
  print("Dim of annot:")
  print(dim(a))
  return(a)
 })
 names(annot) <- names(info)
 print("Are annotation and data transcript ids equal?")
 lapply(seq(1, length(annot)), function(i) {
   print(all(rownames(tx_lst[[i]]$counts) == annot[[i]]$TXNAME))
 })
 #generate count matrix and apply DRIMseq filter (for this we create a DRIMSeq obj (even when running DEXSeq))
 count_data <- lapply(seq(1, length(tx_lst)), function(i) {
  counts <- data.frame(gene_id = annot[[i]]$GENEID, feature_id = annot[[i]]$TXNAME, tx_lst[[i]]$counts, check.names = F)
    #quick fix for the problem of "data.frame()" replacing the hyphens in the sample_ids in the columnnames of the countmatrix with a dot, and subsequently resulting in a mismatch of sample_ids of meta data and count matrix
   # counts  <- counts %>% dplyr::rename_(.dots=setNames(names(.), gsub("\\.", "-",names(.))))
  return(counts)
 })
 Ds <- lapply(seq(1, length(info)), function(i) {
  d <- DRIMSeq::dmDSdata(counts = count_data[[i]], samples = info[[i]])
  print(d)
  return(d)
 })
 names(Ds) <- names(info)
 
 filter_params <- as.numeric(split_colon(param_obj$f))
 Ds <- lapply(Ds, function(d) {
  n.small <- min(table(DRIMSeq::samples(d)$condition))
  n <- nrow(DRIMSeq::samples(d))
  filt_info <- c("n.small" = n.small, "n" = n, "min_feature_expr" = filter_params[1], "min_feature_prop" = 0.1, "min_gene_expr" = filter_params[2])
  d_before_filt <- d
  d <- DRIMSeq::dmFilter(d,
   min_samps_feature_expr = n.small, min_feature_expr = filter_params[1],
   min_samps_feature_prop = n.small, min_feature_prop = 0.1,
   min_samps_gene_expr = n, min_gene_expr = filter_params[2])
  return(list("Ds" = d, "Ds_unfilt" = d_before_filt, "filt_info" = filt_info))
 })
 names(Ds) <- names(info)
 return(list(tx_lst = tx_lst, annot = annot, info = info, obj = Ds))
}
