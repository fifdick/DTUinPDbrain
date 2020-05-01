#!/usr/bin/env Rscript 
library(DTU)

suppressPackageStartupMessages(require(optparse)) # don't say "Loading required package: optparse"
option_list = list(
  make_option(c("-g", "--geneListPaths"), action="store", default=NA, type='character',
              help="path to gene lists separated by colon :"),
 make_option(c("-n", "--cohortNames"), action="store", default="NBB:Norway", type='character',
              help="provide cohortNames separated by colon :"),
 make_option(c("-f", "--filterInfo"), action="store", default=NA, type='character',
              help="list of filter parameters used with DRIMSeq"),
  make_option(c("-u", "--unfiltDs"), action="store", default=NA, type='character',
              help="raw DRIMSeq object, just scaled counts after tximport, no filtering applied"), 
make_option(c("-c", "--conditionNames"), action="store", default="Control:Case", type='character',
             help="provide condition names separated by colon :"),
 make_option(c("-o", "--outDir"), action="store", default="./", type='character',
              help="porivde path where output should be written to"),
 make_option(c("-d", "--DsObjPaths"), action="store", default=NA, type='character',
              help="porivde path to pipeline objects"),
 make_option(c("-s", "--SsObjPaths"), action="store", default=NA, type='character',
              help="porivde path to pipeline objects"),
 make_option(c("-p", "--pheno"), action="store", default="./", type='character',
              help="porivde path to phenoData file"),
 make_option(c("-v", "--version"), action="store", default=75,
             help="which ensembl version should be used for annotating [default %default]"))

opt = parse_args(OptionParser(option_list=option_list))

####################################LOAD LIBRARIES############################################
if(!require(ddpcr)){install.packages(c("ddpcr", "v1.9"),repos='http://cran.us.r-project.org')}
load_dependencies<-function()
{
  if(!require(readr)){install.packages("readr")}
 if(!require(DRIMSeq)){
    source("https://bioconductor.org/biocLite.R")
    biocManager::install("DRIMSeq")}
  if(!require(stageR)){source("https://bioconductor.org/biocLite.R")
    biocManager::install("stageR")}
  if(!require(DEXSeq)){source("https://bioconductor.org/biocLite.R")
    biocManager::install("DEXSeq")}
  }
quiet(load_dependencies())

# read gene lists from tools (after stageR)
objs <- lapply(split_colon(opt$g),readRDS)
print(opt$g)
print(split_colon(opt$g))
objs <- lapply(objs,extract_sig)
names(objs) <- c("drim_gene_lsts","dex_gene_lsts")     

#split by cohorts
objs_byCohort <-split_by_cohorts(cohort_names=split_colon(opt$cohortNames),objs=objs)

#intersect and get gene annotation (additional info) from ensembldb
toolInt<-apply_tool_intersect(objs_by_cohort=objs_byCohort,tx=T,opt=opt)

# Ds unfiltered, both tools together 
Ds_unfilt <- lapply(split_colon(opt$u),readRDS)
names(Ds_unfilt) <- c("Ds_drim_unfilt","Ds_dex_unfilt")

#DS both tools together as DS$tool$cohort
Ds <- lapply(split_colon(opt$DsObjPaths),readRDS)
names(Ds) <- c("Ds_drim","Ds_dex")

#filtInfo
filtInfo <- lapply(split_colon(opt$f),readRDS)
names(filtInfo)<-c("filtInfo_drim","filtInfo_dex")

#Results dfs separated toolRes$cohort
dexRes<-lapply(Ds$Ds_dex,function(ds){
return( as.data.frame(DEXSeqResults(ds,independentFiltering=F)))})

drimRes<-lapply(Ds$Ds_drim,function(ds){
return( as.data.frame(DRIMSeq::results(ds,level="feature")))})

saveRDS(objs_byCohort,file=paste0(opt$o,"objs_byCohort.rds"))
saveRDS(toolInt,file=paste0(opt$o,"drim_dex_toolInt.rds"))
saveRDS(dexRes,file=paste0(opt$o,"dexRes.rds"))
saveRDS(drimRes,file=paste0(opt$o,"drimRes.rds"))
saveRDS(Ds,file=paste0(opt$o,"Ds.rds"))
saveRDS(Ds_unfilt,file=paste0(opt$o,"Ds_unfilt.rds"))
saveRDS(filtInfo,file=paste0(opt$o,"filtInfo.rds"))

