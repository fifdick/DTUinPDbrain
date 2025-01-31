% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prep_dtu.R
\name{prep_dtu}
\alias{prep_dtu}
\title{prepare data for DTU analysis}
\usage{
prep_dtu(param_obj)
}
\arguments{
\item{param_obj}{opt parse object parsed from the command line (runDTU.sh)}
}
\value{
tx_lst is a list of tximport results (one for each cohort)
annot is a list of annotation (one for each cohort)
Ds is an object for each cohort that contains the filtered DRIMSeq obj, the pre-filter DRIMSeq obj and the filter parameters used
}
\description{
Used in DTU::run_dex and DTU::run_drim.
}
