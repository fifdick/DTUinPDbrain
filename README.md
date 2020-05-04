# Differential transcript usage in the Parkinson's disease brain 

This repository holds all code and data to reproduce results, figures and tables presented in "Differential transcript usage in the Parkinson's disease brain". The main analysis was performed with
[DEXSeq](https://bioconductor.org/packages/release/bioc/html/DEXSeq.html), [DRIMSeq](https://bioconductor.org/packages/release/bioc/html/DRIMSeq.html) and [stageR](https://www.bioconductor.org/packages/release/bioc/html/stageR.html). The application of these tools for DTU analysis was performed as suggested in [Swimming downstream: statistical analysis of differential transcript usage following Salmon quantification.](https://www.ncbi.nlm.nih.gov/pubmed/30356428).  

## Prerequisites

* [R >= 3.6.0](https://www.r-project.org/)
* Rscript in "$PATH" 
* [Bash](https://www.gnu.org/software/bash/) 
* All dependencies listed in DESCRIPTION of the DTU R package

## How to..  

### Use uploaded result objects (.rds) to generate figures and tables:    
In your terminal execute:
```
git clone https://github.com/fifdick/DTU_in_PD_brain.git
cd DTU_in_PD_brain
tar -xzf ./results/rds/rds.tar.gz
R CMD INSTALL DTU 
R CMD check DTU
R CMD build DTU
```
In your R editor open `makeFigures.Rmd` and follow the instructions there.

### Rerun the analysis from scratch:  

```
git clone https://github.com/fifdick/DTU_in_PD_brain.git
cd DTU_in_PD_brain
# unpack the reference .gtf file 
tar -xzf ./referenceData/extract_me_gtf.gz
# install the R package to your R
R CMD INSTALL DTU 
R CMD check DTU
R CMD build DTU
# make the pipeline executable
chmod +x runDTU.sh
# potentially change parameters according to the instructions in runDTU.sh
# make sure Rcript is in your path as it will be called by the bash script
./runDTU.sh
```  

* The analysis is complete when "Sucessfully intersected gene lists" appears in stdout. Results are in `./results/rds/`, direct tool results are in the respective folders in `./results/`.  
* Open `./makeFigures.Rmd` in your R editor and follow the instructions there.


## Content

* Transcript quantifications obtained with [Salmon](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5600148/)
for all samples in the study are in `./rawData/`. Where each subdirectory holds the estimates for the sample identified with the subdir name.  
* For annotation of the transcripts during DTU analysis with DEXSeq and DRIMSeq (when applying tximport::tximport for the tx2gene parameter), we used the reference file provided in `./referenceData`.
(Zipped to be able to fit on github.)  
* Meta data including all covariates used in the model design are supplied in the .csv in `./metaData". The sample IDs match the subdirectory names in `./rawData/`. The cell type estimations listed as covariates ("Microglia\_Genes" and "Oligo\_Genes") have been estimated as described in the paper. The code for the calculation can be found on: https://git.app.uib.no/neuromics/cell-composition-rna-pd. The estimation has be performed on the replication and discovery cohort separately, and all samples of each cohort where included.  
* The log of final run of the pipeline can be found in `./logs`.
* Zipped .rds objects are in `./results/rds`. These hold all results that are necessary to generate tables and figures presented in the paper (except the qPCR wetlab results). They were obtained by running the script `./runDTU.sh`. All arguments and parameters that were used to generate these results can be found in the log file `./log/04-05-20.txt`.  
* The directory `./DTU` contains all helper functions collected into an R package. These are used in both the general pipeline (runDTU.sh) and in the .Rmd file that generates tables and figures. 
* The directory `./scripts` holds Rscripts which are called by `./runDTU.sh` to perferm the general DTU analysis. Basically just wrapper functions around DRIMSeq and DEXSeq functions. The results of these are then assembled by `./scripts/bindToolRes.R` to generate the results found in `./results/rds`. These hold all R objects needed to create figures and tables in `./makeFigures.Rmd`.  
* Code used to generate plots and tables is in `./makeFigures.Rmd`.  

## Built with 

### Bash
5.0.3(1)-release

### R sessionInfo() output  


