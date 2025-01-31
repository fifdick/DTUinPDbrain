% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_gene_info.R
\name{get_gene_info}
\alias{get_gene_info}
\title{Retrieve gene annotation from ensembldb}
\usage{
get_gene_info(
  id_list,
  version = "75",
  g_source = "ensembldb",
  tx = FALSE,
  exons = F,
  uniprot = FALSE,
  symbol = FALSE,
  promoters = FALSE,
  upBP = 2000,
  downBP = 400
)
}
\arguments{
\item{id_list}{character vector of ids used to retrieve gene annotation, if tx==T then it must be ensembl transcript ids that are valid in the ensembl release that is specified with the parameter version. Else it must be ensembl gene ids.}

\item{version}{string specifying ensembl release version (default is 75 which was used within this analysis, the last of grch37)}

\item{g_source}{char string specifing the database used to retrieve annotation (this is now redundant, removed functionality to use ucsc)}

\item{tx}{bool specifying whether the vector provided with id_list are transcript (tx==TRUE) or gene ids}
}
\value{
dataframe of requested gene annotation.
}
\description{
Especially used for making the paper figures. To retrieve transcript biotypes and gene names based on either transcript id or gene id.
}
