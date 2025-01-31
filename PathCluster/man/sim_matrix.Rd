% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_matrix.R
\name{sim_matrix}
\alias{sim_matrix}
\title{Calculate pair-wise similarity between pathways based on their members}
\usage{
sim_matrix(alist)
}
\arguments{
\item{alist}{\code{list()} with one element for each pathway, named with the
corresponding pathway name. Each element on the list is a vector of gene
names/IDs ('\code{character})}
}
\value{
Returns a (square) \code{data.frame} with the numeric similarity
    values (Kappa scores)
}
\description{
From a list of pathways (with gene ID/names as character vectors) outputs a
similarity matrix between each pair of pathways
}
