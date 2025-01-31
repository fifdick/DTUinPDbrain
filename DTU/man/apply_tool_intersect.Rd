% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apply_tool_intersect.R
\name{apply_tool_intersect}
\alias{apply_tool_intersect}
\title{apply tool intersection per cohort}
\usage{
apply_tool_intersect(objs_by_cohort, tool_idx = NULL, tx = T, opt)
}
\arguments{
\item{objs_by_cohort}{: objs listed by cohort and each tool result per cohort (after DTU::split_by_cohorts() was used)}

\item{tx}{: level of list intersection (TRUE || FALSE ) tx vs gene level, default=T}

\item{opt}{opt parse object that was generated in "bindToolRes.R" when called through bash script "runDTU.sh"}

\item{tool_idx:}{which tools to intersect, default:NULL, all tools are being intersected. (Was helpful before when we also used IsoformSwitchAnalyserR. Might become usefull when we try other tools.)}
}
\value{
list with one dataframe per cohort, holding a DTU events (as returned from stageR) that overlap between the tools
}
\description{
Used in "bindToolRes.R" to put together the tool specific result dataframes which were generated by the pipeline
}
