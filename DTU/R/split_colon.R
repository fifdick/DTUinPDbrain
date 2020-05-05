#' split string by colon 
#' Really unnecessary function that probably already exists. Used to split the input argument strings from the commandline in runDTU.sh. Arguments like for example cohort_names are specified in the bash script with "discovery:replication". 
#' @param character string that should be split by colon 
#' @return named character vector  
#' @examples 
#' split_colon("Control:Case")
#' split_colon("gender:age:rin")
#' @export
split_colon <- function(string) {
    vec <- unlist(strsplit(string, ":", fixed = T))
    names(vec) <- vec
    return(vec)
}


