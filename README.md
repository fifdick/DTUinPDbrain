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
7z x -so ./results/rds/rds.7z
tar -zxvf ./results/external/dge_results/DESeqOut_CT.tar.gz
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
# make sure Rscript is in your path as it will be called by the bash script
./runDTU.sh
```  

* The analysis is complete when "Sucessfully intersected gene lists" appears in stdout. Results are in `./results/rds/`, direct tool results are in the respective folders in `./results/`.    

```
# extract DGE results, needed to compare and make figures and tables
tar -zxvf ./results/external/dge_results/DESeqOut_CT.tar.gz
```

* Open `./makeFigures.Rmd` in your R editor and follow the instructions there.


## Content

* Transcript quantifications obtained with [Salmon](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5600148/)
for all samples in the study are in `./rawData/`. Where each subdirectory holds the estimates for the sample identified with the subdir name.  
* For annotation of the transcripts during DTU analysis with DEXSeq and DRIMSeq (when applying tximport::tximport for the tx2gene parameter), we used the reference file provided in `./referenceData`.
(Zipped to be able to fit on github.)  
* Meta data including all covariates used in the model design are supplied in the .csv in `./metaData`. The sample IDs match the subdirectory names in `./rawData/`. The cell type estimations listed as covariates ("Microglia\_Genes" and "Oligo\_Genes") have been estimated as described in the paper. The code for the calculation can be found on: https://git.app.uib.no/neuromics/cell-composition-rna-pd. The estimation has be performed on the replication and discovery cohort separately, and all samples of each cohort where included.  
* The log of final run of the pipeline can be found in `./logs`.
* Zipped .rds objects are in `./results/rds`. These hold all results that are necessary to generate tables and figures presented in the paper (except the qPCR wetlab results). They were obtained by running the script `./runDTU.sh`. All arguments and parameters that were used to generate these results can be found in the log file `./log/04-05-20.txt`.  
* The directory `./DTU` contains all helper functions collected into an R package. These are used in both the general pipeline (runDTU.sh) and in the .Rmd file that generates tables and figures. 
* The directory `./scripts` holds Rscripts which are called by `./runDTU.sh` to perferm the general DTU analysis. Basically just wrapper functions around DRIMSeq and DEXSeq functions. The results of these are then assembled by `./scripts/bindToolRes.R` to generate the results found in `./results/rds`. These hold all R objects needed to create figures and tables in `./makeFigures.Rmd`.  
* Differential gene expression results from ["Common gene expression signatures in Parkinson’s disease are driven by changes in cell composition"](https://actaneurocomms.biomedcentral.com/articles/10.1186/s40478-020-00932-7) which are used to compare DTU results to DGE results are in `./results/external/dge_results/`
* qPCR results are summerised in a .csv in `./results/external/qPCR/`
* Code used to generate plots and tables is in `./makeFigures.Rmd`.  

## Built with 

### Bash
5.0.3(1)-release

### R sessionInfo() output  

```
> sessionInfo()
R version 3.6.0 (2019-04-26)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Pop!_OS 19.04

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.8.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.8.
0

locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8
 [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8
 [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C

attached base packages:
 [1] stats4    parallel  grid      stats     graphics
 [6] grDevices utils     datasets  methods   base

other attached packages:
 [1] wiggleplotr_1.8.0
 [2] DEXSeq_1.30.0
 [3] RColorBrewer_1.1-2
 [4] AnnotationDbi_1.46.0
 [5] DESeq2_1.24.0
 [6] SummarizedExperiment_1.14.0
 [7] DelayedArray_0.10.0
 [8] GenomicRanges_1.36.0
 [9] GenomeInfoDb_1.20.0
[10] IRanges_2.18.1
[11] S4Vectors_0.22.0
[12] Biobase_2.44.0
[13] BiocGenerics_0.30.0
[14] BiocParallel_1.18.0
[15] DRIMSeq_1.12.0
[16] reshape2_1.4.3
[17] kableExtra_1.1.0
[18] gridExtra_2.3
[19] aggregation_1.0.1
[20] ggridges_0.5.1
[21] viridis_0.5.1
[22] viridisLite_0.3.0
[23] dplyr_0.8.3
[24] tidyr_1.0.2
[25] readr_1.3.1
[26] Unicode_12.0.0-1
[27] knitr_1.23
[28] STRINGdb_1.24.0
[29] DTU_1.0
[30] ggrepel_0.8.1
[31] matrixStats_0.54.0
[32] cowplot_1.0.0
[33] VennDiagram_1.6.20
[34] futile.logger_1.4.3
[35] ggplotify_0.0.4
[36] patchwork_1.0.0.9000
[37] ggpubr_0.2.1
[38] magrittr_1.5
[39] ggplot2_3.2.1
[40] nvimcom_0.9-83

loaded via a namespace (and not attached):
  [1] backports_1.1.5        Hmisc_4.2-0
  [3] plyr_1.8.4             igraph_1.2.4.1
  [5] lazyeval_0.2.2         splines_3.6.0
  [7] digest_0.6.22          htmltools_0.4.0
  [9] gdata_2.18.0           checkmate_1.9.4
 [11] memoise_1.1.0          cluster_2.0.7-1
 [13] limma_3.40.4           Biostrings_2.52.0
 [15] annotate_1.62.0        prettyunits_1.0.2
 [17] colorspace_1.4-1       blob_1.2.0
 [19] rvest_0.3.4            xfun_0.8
 [21] crayon_1.3.4           RCurl_1.95-4.12
 [23] genefilter_1.66.0      survival_2.43-3
 [25] glue_1.3.1             hash_2.2.6.1
 [27] gtable_0.3.0           zlibbioc_1.30.0
 [29] XVector_0.24.0         webshot_0.5.1
 [31] scales_1.0.0           futile.options_1.0.1
 [33] DBI_1.0.0              edgeR_3.26.5
 [35] Rcpp_1.0.3             plotrix_3.7-6
 [37] xtable_1.8-4           progress_1.2.2
 [39] htmlTable_1.13.1       gridGraphics_0.5-0
 [41] foreign_0.8-71         bit_1.1-14
 [43] Formula_1.2-3          sqldf_0.4-11
 [45] htmlwidgets_1.5.1      httr_1.4.0
 [47] gplots_3.0.1.1         acepack_1.4.1
 [49] pkgconfig_2.0.3        XML_3.98-1.20
 [51] nnet_7.3-12            locfit_1.5-9.1
 [53] tidyselect_1.0.0       rlang_0.4.6
 [55] munsell_0.5.0          tools_3.6.0
 [57] gsubfn_0.7             RSQLite_2.1.2
 [59] evaluate_0.14          stringr_1.4.0
 [61] bit64_0.9-7            caTools_1.17.1.2
 [63] purrr_0.3.2            formatR_1.7
 [65] xml2_1.2.0             biomaRt_2.40.3
 [67] compiler_3.6.0         rstudioapi_0.10
 [69] png_0.1-7              ggsignif_0.5.0
 [71] tibble_2.1.3           statmod_1.4.32
 [73] geneplotter_1.62.0     stringi_1.4.3
 [75] lattice_0.20-38        Matrix_1.2-15
 [77] vctrs_0.2.4            pillar_1.4.2
 [79] lifecycle_0.1.0        BiocManager_1.30.4
 [81] Tmisc_0.1.22           data.table_1.12.2
 [83] bitops_1.0-6           R6_2.4.1
 [85] latticeExtra_0.6-28    hwriter_1.3.2
 [87] KernSmooth_2.23-15     lambda.r_1.2.3
 [89] gtools_3.8.1           assertthat_0.2.1
 [91] chron_2.3-53           proto_1.0.0
 [93] withr_2.1.2            Rsamtools_2.0.0
 [95] GenomeInfoDbData_1.2.1 hms_0.5.0
 [97] rpart_4.1-13           rmarkdown_1.14
 [99] rvcheck_0.1.8          base64enc_0.1-3
```
