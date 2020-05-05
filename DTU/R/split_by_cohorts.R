#' Reorder nested list with DTU objects 
#'
#' Used in the Rscript "bindToolRes.R". Transfrom listing order from 'per tool' to 'per cohort'.
#' Input obj consists of a list for each tool, and the tool specific list holds gene-lists for each cohort.
#' The function swapps this around such that the outer of the nested list is a list of cohorts within which the two tool specific gene results are listed
#' @param cohort_names char vector with names of cohorts (as used in the pipeline (bash))
#' @param objs the nested list that should be reorderd 
#' @return nested list with the outer list named by cohorts, and a list with specific tool results in each of those
#' @export
split_by_cohorts <- function(cohort_names, objs) {
    objs_by_cohort <- lapply(cohort_names, function(cohort) {
        lapply(objs, function(o) {
            o[[cohort]]
        })
    })
    names(objs_by_cohort) <- cohort_names
    return(objs_by_cohort)
}
