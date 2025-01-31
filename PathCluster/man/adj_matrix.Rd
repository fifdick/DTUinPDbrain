% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adj_matrix.R
\name{adj_matrix}
\alias{adj_matrix}
\title{Calculate a pair-wise adjacency matrix for a list of pathways}
\usage{
adj_matrix(alist, alpha = 0.05)
}
\arguments{
\item{alist}{\code{list()} with one element for each pathway. Each element
is a binary named vector, where the names are the names/IDs of the genes
that belong to the pathway and the values represent nominal significance
(or lack thereof)}

\item{alpha}{hard threshold to transform the p-values of the test to a
binary adjacency matrix}
}
\value{
Returns a (square) \code{data.frame} with the p-values of the
    pair-wise intersections using a hypergeometric test
}
\description{
From a list of pathways (in which each element is a named vector of genes
    with binary values if the genes were differentially expressed or not),
    calculates a binary pairwise adjacency matrix based on the significance
    of a hypergeometric test of the intersections
}
