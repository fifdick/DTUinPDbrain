# Differential transcript usage in the Parkinson's disease brain 

This repository holds all code and data to reproduce results, figures and tables presented in "Differential transcript usage in the Parkinson's disease brain". The main analysis was performed with
[DEXSeq](https://bioconductor.org/packages/release/bioc/html/DEXSeq.html), [DRIMSeq](https://bioconductor.org/packages/release/bioc/html/DRIMSeq.html) and [stageR](https://www.bioconductor.org/packages/release/bioc/html/stageR.html).  

## Prerequisites

* R >= 3.6.0  
* Rscript in path 
* bash 
* All dependencies listed in DESCRIPTION of the DTU R package

## How to
```

```

## Content

* Transcript quantifications obtained with [Salmon](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5600148/)
for all samples in the study are in `./rawData/`. Where each subdirectory holds the estimates for the sample identified with the subdir name.  
* For annotation of the transcripts during DTU analysis with DEXSeq and DRIMSeq (when applying tximport::tximport for the tx2gene parameter), we used the reference file provided in "referenceData".
(Zipped to be able to fit on github)  
* Meta data including all covariates used in the model design are supplied in the csv in "metaData". The sample IDs match the subdirectory names in "rawData/". The cell type estimations listed as covariates ("Microglia\_genes" and "Oligo\_genes") have been estimated as described in the paper. The code for the calculation can be found on: https://git.app.uib.no/neuromics/cell-composition-rna-pd. The estimation has be performed on the replication and discovery cohort separately, and all samples of each cohort where included.  
* The log of final run of the pipeline can be found in "logs"
* Zipped .rds objects are in "results". These hold all results that are necessary to generate tables and figures presented in the paper (except the qPCR wetlab results). They were obtained by running the script
* "runDTU.sh". All arguments and parameters that were used to generate these results can be found in the log file "log/04-05-20.txt".  
* The directory "Rcode" contains a R packages "DTU" that contains all helper functions which are used in both the general pipeline (runDTU.sh) and in the .Rmd file that generates tables and figures. 
* The directory "scripts" holds 3 Rscripts which are called by "runDTU.sh" to perferm the general DTU analysis. Basically just wrapper functions around DRIMSeq and DEXSeq functions. The results of these are then assembled by "bindToolRes.R" to generate the results found in "rds". These hold all R objects needed to create figures and tables in "makeFigures.Rmd".  
* Code used to generate plots and tables is in makeFigures.Rmd  

## Built with 

### bash
5.0.3(1)-release

### R sessionInfo() output  


