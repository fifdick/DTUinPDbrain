#' Wrapper for plotting treemap of the clusters of pathways
#'
#' Plots the treemap of the clusters of pathways
#'
#' @param dtf \code{data.frame()} with the following columns: (1) "items" -
#'     pathway names; (2) "title" - name of the cluster (representative
#'     pathway); (3) "abslog10" - -log10(pvalue) for the pathway (required)
#'
#' @param index names of the columns for the hierarchical mapping
#'
#' @param ... parameters passed to \code{treeplot()} function
#'
#' @return None (only plots)
#'
#' @export

print_map <- function(
    dtf, index=c("title", "items"), vSize="abslog10", type="categorical",
    vColor="title", title="", inflate.labels=FALSE, lowerbound.cex.labels=0,
    bg.labels="#CCCCCCAA", position.legend="none", palette="PiYG",
    fontcolor.labels=c("black","white"), border.col=c("black","white"),
    border.lwds=c(5,2), ...){
    
    if(grepl('^GO', dtf[1,index[1]])){
        clusters <- nameGO(dtf)
    }

    treemap::treemap(
       dtf = dtf,
       index = index, 
       vSize = vSize,
       type = type,
       vColor = vColor,
       title = title,
       inflate.labels = inflate.labels,                # set this to TRUE for space-filling group labels - good for posters
       lowerbound.cex.labels = lowerbound.cex.labels,  # try to draw as many labels as possible (still, some small squares may not get a label)
       bg.labels = bg.labels,                          # define background color of group labels
                                                       # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
       position.legend = position.legend,
       palette = palette,
       fontcolor.labels = fontcolor.labels,
       border.col = border.col,                        # Color of borders of groups, of subgroups, of subsubgroups ....
       border.lwds = border.lwds,
       ...
    )
}
