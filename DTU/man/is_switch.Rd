% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_switch.R
\name{is_switch}
\alias{is_switch}
\title{Identify isoform switch}
\usage{
is_switch(df, gene)
}
\arguments{
\item{df}{dataframe resulting from DTU::plot_genes_cohort(). (p[[1]]$data)}

\item{gene}{string gene identifier (Gene symbol)}
}
\value{
list with bool specifying whether or not there is an isoform switch, and if yes, the transcript id and biotype of the top (relative) expressed transcript of the control and case group separately
}
\description{
Defined a gene to have an isoform switch, if the top expressed transcript differs between groups. Requires a minimum of 2 DTU events per gene.
}
