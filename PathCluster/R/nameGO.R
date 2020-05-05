#' Fetches GO pathway names with the GO IDs
#'
#' @param clusters object with an \code{items} slot which contains the GO IDs
#'     (e.g. "GO:0000009")
#'
#' @return Returns a \code{data.frame} GO names and titles
#'
#' @export

nameGO <- function(clusters) {
    goterms <- AnnotationDbi::Term(GO.db::GOTERM)[clusters$items]
    dfGo <- data.frame(title=names(goterms), TitleName=goterms)
    df <- merge(clusters, dfGo, by="title", all.x=TRUE, all.y=FALSE)
    names(df)[names(df) == 'title'] <- 'GO_idTitle'
    names(df)[names(df) == 'TitleName'] <- 'title'
    return(df)
}
