#' Create phenotype dataframe
#' 
#' TODO: REWRITE THIS!!! simple dplyr and split...Create the sample information, read from the metadatafile 
#' @param pheno character string representing path to phenotype csv file
#' @param cohort_names vector of strings containing cohort names
#' @param conditions vector of strings containing names of group conditions (first control then case: c("CT,"CS"))
#' @param batch_vars char vector with column names of covariates that should also be read from file
#' @return list of df with minimum 3 columns: sample_id, condition, rin. And additional covariate columns Length of list == length of \code{cohort_names}
#' @importFrom dplyr %>%
#' @export
create_sample_info <- function(pheno, cohort_names, conditions, batch_vars=NULL, cohort_col_name="procedence") {
 pheno_data <- readr::read_delim(pheno, "\t", escape_double = FALSE, trim_ws = TRUE)
 rownames(pheno_data) <- pheno_data$sample_id
 #remove the non PD-CT samples
 pheno_data <- subset(pheno_data, condition %in% conditions)
 rownames(pheno_data) <- pheno_data$sample_id
 #groups sample_ids by procedence (according to cohort names))
 cohort_ids <- lapply(cohort_names, function(cn) {
  cn_ids <- subset(pheno_data, procedence == cn)$sample_id #as.data.frame(pheno_data[pheno_data[, cohort_col_name] == cn, "sample_id"]) 
  return(cn_ids)
 })
 names(cohort_ids) <- cohort_names
 Y <- lapply(cohort_ids, function(ci) {
  if (length(conditions) == 2) {
   y <- pheno_data[ci, ]$condition
   #TODO MUST BE ADJUSTED TO MODIFYABLE CONDITION ARG
   #select only samples with specified group conditions 
   #the first condition in condition names vector given by user will be 0 reference
   y <- sapply(y, function(cond) (ifelse(grepl(cond, as.character(conditions[1])), 0, 1)))
   y <- as.factor(y)
   print(y)
   y <- stats::relevel(y, ref = "0")
  } else {
   y <- as.numeric(pheno_data[ci, ]$condition)
  }
  names(y) <- ci
  return(y)
 })
 names(Y) <- cohort_names
 info <- lapply(seq(1, length(Y)), function(i) {
  data.frame("sample_id" = cohort_ids[[i]], "condition" = Y[i], "rin" = as.numeric(pheno_data[cohort_ids[[i]], ]$rin))
 })
 info <- lapply(info, function(i) {
  colnames(i) <- c("sample_id", "condition", "rin")
  if (nrow(i) == 0) {
   stop("Something went wrong during sample info creation, no rows in info DF")
  }
  return(i)
 })
 if (!(is.null(batch_vars))) {
  info <- lapply(info, function(i) {
  df <- pheno_data %>% dplyr::as_data_frame(.) %>% dplyr::select(c("sample_id", unlist(batch_vars)))
  i <- i %>% dplyr::inner_join(df)
  return(i)
  })
 }
 names(info) <- cohort_names
 return(info)
}
