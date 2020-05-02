#!/usr/bin/env Rscript

# DRIMSeq and stageR pipeline for DTU analysis
# make sure you have Rscript in your path and the R package optparse and the R packes DTU (from this git repo) installed
# arguments -f, -n, -c and -b need to be in string form separated by colons. E.g.: -n "cohort1:cohort2", e.g.: -c "Control:Case" e.g.: -b "rin:sex:age_years"
# Run > Rscript drimStage.R -o <outputdir> \
#        -f <min_transcript_expression:min_gene_expression> \
#        -g <reference_gtf_file> \
#        -m <metadatafile> \
#        -s <salmon_outputs_dir/> \
#        -z <tximport_scaling_method> \
#        -n <names_of_cohorts_as_in_metadata_file> \
#        -c <names_of_group_conditions> \
#        -d <bool_postHocFilter> \
#        -t <timestamp_or_outputfile_description> \
#        -b <covariates> \


suppressPackageStartupMessages(require(optparse))

option_list <-  list(
  make_option(c("-x", "--xNoTrFilt"), action = "store", default = NULL, type = "character", help = "Sets the number of transcripts a gene can maximally have to not be filtered out prior to analyses. To reduce complexity."),
  make_option(c("-t", "--today"), action = "store", default = NA, type = "character", help = "set timestamp to be used in the output filenames"),
  make_option(c("-s", "--salmonDir"), action = "store", default = NA, type = "character",
              help = "provide the path to your salmon data directory (containing each sample data in separate folders"),
  make_option(c("-m", "--metaData"), action = "store", type = "character", default = NA,
              help = "provide the path to your phenotype information file (which should include columns sample_id,rin, condition and procedence, if the levels of condition dont match (Control, Case) specify the condition names with the option -c."),
  make_option(c("-v", "--verbose"), action = "store_true", default = TRUE,
              help = "Should the program print extra stuff out? [default %default]"),
  make_option(c("-q", "--quiet"), action = "store_false", dest = "verbose",
              help = "Make the program not be verbose."),
  make_option(c("-e", "--ensemblVersion"), action = "store", default = "75",
              help = "which version of ensembl db should be used to retrieve gene names or in some cases the whole annotation if gtf file isnt used [default %default]"),
  make_option(c("-n", "--namesCohorts"), action = "store",
              help = "Define cohort names in string format separated by colons (e.g. \"NBB:PV\" [default %default]"),
  make_option(c("-c", "--conditionGroups"), action = "store", default = "Control:Case",
              help = "Define condition names in correct order (reference group will be the first, i.e. in default control). Specify as string separated by a colon (\"Control:Case\"  [default %default]"),
  make_option(c("-b", "--batchVars"), action = "store", default = "rin",
              help = "Define variable names of covariates (must be colnames in metadata file) to be included in model design (if more than one, separate by colon (\"rin:age\")) [default %default]"),
make_option(c("-g", "--gtfFile"), action = "store", default = "NA",
              help = "provide path to the gtf file used for annotating the transcripts, should contain all transcripts that are listed in the quant.sf files produced by salmon, preferably the same gtf file used for salmon"),
make_option(c("-f", "--filterParams"), action = "store", default = "10:10",
              help = "option to define filter threshold, enter minimum transcript expression followed by minimum gene expression separated by colon. This will be used as arguments for the DRIMSeq::DMFilter function (min_feature_expr and min_samps_gene_expr). Specify as string separated by a colon.  [default %default]"),
make_option(c("-o", "--outDir"), action = "store", default = "./",
              help = "provide directory in which all output objects will be stored [default %default]"),
make_option(c("-d", "--postHocFilter"), action = "store", default = "TRUE",
            help = "Boolean, TRUE: applies post hoc filter, setting pvalues to 1 for all transcripts in sample that have lower than 0.1 sd before 2stage testing with stageR is applied (according to \"Swimming downstream: statistical analysis of differential transcript usage following Salmon quantification\" Love et. al [default %default]")  
,
make_option(c("-z", "--zcalingMethod"), action = "store", default = "scaledTPM",  type = "character", help = "String specifying scaling method used to scale TPMs with tximport (possible: <scaledTPM,dtuScaledTPM>[default %default]")  
)

# parse commandline arguments
opt <- parse_args(OptionParser(option_list = option_list))
if (opt$v) {
 print(opt)
}

# load dependencies
if (!require(ddpcr)) { 
 install.packages(c("ddpcr", "v1.9"), repos = "http://cran.us.r-project.org")
}
load_dependencies<-function() {
 if (!require(dplyr)) { install.packages("dplyr") }
 if (!require(readr)) { install.packages("readr") }
 if (!require(rlang)) { install.packages("rlang") }
 if (!require(DTU)) { install.packages("../DTU/DTU_1.0.tar.gz", repos = NULL) }
 if (!requireNamespace("BiocManager", quietly = TRUE)) { install.packages("BiocManager") }
 if (!require(DRIMSeq)) { BiocManager::install("DRIMSeq") }
 if (!require(stageR)) { BiocManager::install("stageR") }
}
ddpcr::quiet(load_dependencies())

# create output directory structure
print("creating output dirs:")
# for the output of the DTU analysis tool
dir.create(file.path(opt$o, "Ds"))
print(file.path(opt$o, "Ds"))
# for the ouput of stageR (and its object)
dir.create(file.path(opt$o, "Ss"))
print(file.path(opt$o, "Ss"))
# for the result gene lists after stageR but with both sig. and non significant transcripts
dir.create(file.path(opt$o, "genes"))
print(file.path(opt$o, "genes"))

# tool specific helper functions
no.na <- function(x) ifelse(is.na(x), 1, x)

# extract DRIMSeq result to be able to use stageR
extract_res <- function(Ds) {
# extract results from drimseq obj, on gene and transcript level, set NA pvalues to 1
 Rs <- lapply(Ds,function(d) {
  res <- DRIMSeq::results(d)
  res.txp <- DRIMSeq::results(d, level = "feature")
  res$pvalue <- no.na(res$pvalue)
  res.txp$pvalue <- no.na(res.txp$pvalue)
  return(list("txp" = res.txp, "res" = res))
 })
 names(Rs) <- names(Ds)
 return(Rs)
}

# helper function for optional postHoc filter
# as described in "Swimming downstream: statistical analysis of differential transcript usage following Salmon quantification" Love et al. 
getSampleProportions <- function(d) {
 cts <- as.matrix(subset(DRIMSeq::counts(d), dplyr::select = -c(.data$gene_id, .data$feature_id)))
 gene.cts <- rowsum(cts, DRIMSeq::counts(d)$gene_id)
 total.cts <- gene.cts[match(counts(d)$gene_id, rownames(gene.cts)), ]
 cts / total.cts
}

after_drim_filter <- function(Rs,Ds) {
 Rs <- lapply(seq(1, length(Ds)), function(i) {
  prop.d <- getSampleProportions(Ds[[i]])
  res.txp.filt <- Rs[[i]]$txp
  res.txp.filt$prop.sd <- sqrt(rowVars(prop.d))
  res.txp.filt$pvalue[res.txp.filt$prop.sd < .1] <- 1
  res.txp.filt$adj_pvalue <- p.adjust(res.txp.filt$pvalue, method = "BH")
  return(list("res" = Rs[[i]]$res, "txp" = res.txp.filt))
 })
 return(Rs)
}

stage_test <- function(Ds, alpha = 0.05, post_filt = FALSE) {
 #extract result
 Rs <- extract_res(Ds)
 #apply optional postHoc filter
 if  (post_filt == TRUE) {
  Rs <- after_drim_filter(Rs,Ds)
  }
 Ss <- lapply(Rs, function(res) {
  #retrieve the gene-level pvalues from DRIMSeq they are already aggregated
  #but not corrected
  pScreen <- no.na(res$res$pvalue)
  names(pScreen) <- res$res$gene_id
  #get the transcript-level pvalues
  pConfirmation <- matrix(res$txp$pvalue, ncol = 1)
  rownames(pConfirmation) <- res$txp$feature_id
  tx2gene <- res$txp[, c("feature_id", "gene_id")]
  stageRObj <- stageR::stageRTx(pScreen = pScreen, pConfirmation = pConfirmation, pScreenAdjusted = FALSE, tx2gene = tx2gene)
  stageRObj <- stageR::stageWiseAdjustment(stageRObj, method = "dtu", alpha = alpha)
  return(stageRObj)
 })
 names(Ss) <- names(Rs)
 return(Ss)
}

# parameters for DRIMSeq::DMFilter
filter_params <- as.numeric(DTU::split_colon(opt$f))
# covariates which will be added to the model
covariates <- as.list(DTU::split_colon(opt$b))

# run all functions
if (!is.na(opt$s) & !is.na(opt$m) & !is.na(opt$g)) {
 #create sample info matrix, create annotation matrix, create and import counts
 data <- DTU::prep_dtu(opt)
 #separate unfiltered objects and filter parameter for downstream analysis
 Ds_unfilt <-  lapply(data$obj, function(cohort_obj) { return(cohort_obj$Ds_unfilt) })
 names(Ds_unfilt) <- names(data$info)
 filt_info <- lapply(data$obj, function(cohort_obj) { return(cohort_obj$filt_info) })
 names(filt_info) <- names(data$info)
 #filter model and fit with DRIMSeq
 Res <- DTU::run_drim(data$info, data$obj, covariates = covariates)
 #set timestamp for saving objects
 today <- as.character(opt$t)
 #run stageR
 Ss <- stage_test(Res, post_filt = as.logical(opt$postHocFilter))
 names(Ss) <- names(Res)
 gene_lsts <- lapply(Ss, function(s) {
  lst <- stageR::getAdjustedPValues(s, onlySignificantGenes = TRUE, order = TRUE)
  return(lst)
 })
 names(gene_lsts) <- names(Ss)
 #save objects
 saveRDS(Res, file = paste0(opt$o, "/Ds/Ds_", today, ".rds"))
 saveRDS(Ds_unfilt, file = paste0(opt$o, "/Ds/Ds_preFilt_", today, ".rds"))
 saveRDS(filt_info, file = paste0(opt$o, "/Ds/paramFilt_", today, ".rds"))
 saveRDS(Ss, file = paste0(opt$o, "/Ss/Ss_", today, ".rds"))
 saveRDS(gene_lsts, file = paste0(opt$o, "/genes/genelist_", today, ".rds"))
 }

