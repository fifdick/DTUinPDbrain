#' @export
load_objs <- function(paths) {
    e <- new.env()
    for (path in paths) {
        load(path, envir = e)
    }
    obj_names <- ls(e)
    objs <- lapply(obj_names, function(str) {
        # obj <- eval(as.name(paste(str)))
        obj <- get(str, envir = e)
        return(obj)
    })
    names(objs) <- obj_names
    return(objs)
}
