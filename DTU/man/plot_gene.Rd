% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_gene.R
\name{plot_gene}
\alias{plot_gene}
\title{Build ggplot obj and dataframe}
\usage{
plot_gene(
  all_vals,
  p_gene_id,
  genes_toplot,
  out = "./",
  plot = T,
  tool,
  only_nom_sig = TRUE,
  selected_samples = NULL,
  want_jitter = TRUE,
  conditions = c("Control", "Case")
)
}
\arguments{
\item{all_vals}{df created in DTU::plot_genes_cohort with fitted and observed fractions and gene level counts}

\item{p_gene_id}{char string identifying gene which should be plotted}

\item{genes_toplot}{gene annotation from which to extract only the transcript annotation specific to the gene in p_gene_id}

\item{tool}{}

\item{only_nom_sig}{}
}
\description{
Build the ggplot object, to plot data extracted with DTU::plot_genes_cohort()i
see plot_genes_cohort to see which inputs the function gets
}
