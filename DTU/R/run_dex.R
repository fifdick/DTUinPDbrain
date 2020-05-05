#' Wrapper for DEXSeq
#'
#' Used in Rscript scripts/"dexStage.R". Builds model matrix from sample info and runs DEXSeq. Countdata is ecxtracted from DRIMSeq object from DTU::prep_dtu() after filtering.
#' @param info <df> or list of df with 3 columns (sample_id,condition,rin) and additional covariates 
#' @param obj list of DRIMSeqData objects created by DTU::prep_dtu(). One for each cohort.
#' @param covariates <list> names of variables to include into model design (actual values must be present in info)
#' @return list of DEXSeq objects after DTU analysis
#' TODO: if covariates is null
#' @export
#' @importFrom stats formula
run_dex <- function(info, obj, covariates=list("rin")) {
 stopifnot(is.list(covariates))
 sapply(info, function(i) {
  stopifnot(covariates %in% colnames(i))
 })
 #prepare covariate names for formula
 covariates <- lapply(covariates, function(c) {
  paste0("+ ", as.character(c), ":exon")
 })
 DXs <- lapply(obj, function(d) {
  design <- formula(paste("~sample + exon ", stringi::stri_flatten(covariates), " + condition:exon"))
  print("DESGIN:")
  print(design)
  reduced_design <- formula(paste("~sample + exon ", stringi::stri_flatten(covariates)))
  print("DESGIN reduced:")
  print(reduced_design)
  sample_data <- DRIMSeq::samples(d$Ds)
  count_data <- round(as.matrix(DEXSeq::counts(d$Ds)[, -c(1:2)]))
  dx <- DEXSeq::DEXSeqDataSet(countData = count_data,
   sampleData = sample_data,
   design = design,
   featureID = DEXSeq::counts(d$Ds)$feature_id, groupID = DEXSeq::counts(d$Ds)$gene_id)
  system.time({
   dx <- DEXSeq::estimateSizeFactors(dx)
   dx <- DEXSeq::estimateDispersions(dx, quiet = FALSE, maxit = 500)
   #dx <- nbinomLRT(dx,reduced=~sample +exon)
   #updated version of "swimming" suggests to use testForDEU insteadt of negative binomial liklihood ratio test
   #actually doesnt matter just internally calls the same nbinom func from DESeq2
   dx <- DEXSeq::testForDEU(dx, reduced = reduced_design)
   #I dont think that i use this fc measure in the analysis. As it is calculated without any covariates by refitting the model. Instead late in the postprocessing I use the coefficient of the condition variable as fold change
   dx <- DEXSeq::estimateExonFoldChanges(dx, fitExpToVar = "condition") 
  })
  return(dx)
 })
 names(DXs) <- names(info)
 return(DXs)
}
