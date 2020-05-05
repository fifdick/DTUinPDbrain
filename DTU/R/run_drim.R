#' Wrapper for DRIMSeq
#'
#' Used in Rscript scripts/"drimStage.R". Builds model matrix from sample info and runs DRIMSeq. Countdata is ecxtracted from DRIMSeq object from DTU::prep_dtu() after filtering.
#' @param info <df> or list of df with 3 columns (sample_id,condition,rin) and additional covariates 
#' @param obj list of DRIMSeqData objects created by DTU::prep_dtu(). One for each cohort.
#' @param covariates <list> names of variables to include into model design (actual values must be present in info)
#' @return list of DRIMSeq objects after DTU analysis
#' TODO: if covariates is null
#' @export
#' @importFrom stats formula

run_drim <- function(info, obj, covariates) {
 stopifnot(is.list(covariates))
 sapply(info, function(i) { stopifnot(covariates %in% colnames(i)) })
 #prepare covariate names for model formula
 covariates <- lapply(covariates, function(c) { paste0(" + ", as.character(c)) })
 Ds <- lapply(obj, function(d) {
  #model design
  design_formula <- formula(paste("~condition", stringi::stri_flatten(covariates)))
  design_full <- stats::model.matrix(design_formula, data = DRIMSeq::samples(d$Ds))
  print("DESIGN: Formula!!")
  print(design_formula)
  print("Design matrix:")
  print(design_full)
  #run analysis
  set.seed(1)
  system.time({
   #estimate precision
   d <- DRIMSeq::dmPrecision(d$Ds, design = design_full)
   #fit model
   d <- DRIMSeq::dmFit(d, design = design_full, verbose = 1)
   #coef indicates which columns should be removed from the full design to get the null model (the condition )
   d <- DRIMSeq::dmTest(d, coef = "condition1")
  })
  return(d)
 })
 names(Ds) <- names(info)
 return(Ds)
}
