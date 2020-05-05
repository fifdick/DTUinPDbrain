#' Remove redundant pathways
#'
#' Pathways are iteratively removed based on a set of criteria until there are
#' no pairs with a similarity above the given threshold
#'
#' @param alist \code{list()} with one element for each pathway, named with the
#'     corresponding pathway name. Each element on the list is a vector of
#'     genes. Only the names of the elements in the vectors will be used (the
#'     gene ID/names)
#'
#' @param threshold similarity threshold (kappa score). Pairs of pathways (or
#'     clusters of pathways) above this threshold will be iteratively merged
#'     until the clusters cannot be further merged. The higher the threshold,
#'     the less the pathways will be clustered
#'
#' @param lowThreshold minimum pathway size (default=50)
#'
#' @param highThreshold maximum pathway size (default=1000)
#'
#' @param pvalues named numeric vector with the pathway-level p-values. The
#'     names correspond to pathway names
#'
#' @return \code{list()} with two components, "clusterMat" with the cluster
#'     assignments and "simMat" with the final similarity matrix
#'
#' @export

redundant_filter <- function (alist, threshold, pvalues, lowThreshold=50, highThreshold=1000){
    message("Calculating adjacency matrix...")
    adjmat <- adj_matrix(alist) # similarity of pathways based on the significance at the gene-level (0/1)

    message("Calculating similarity matrix...")
    simmat <- sim_matrix(alist) # similarity of pathways based on the genes they share (kappas)
    diag(simmat) <- NA

    pathwayNames <- rownames(simmat)   # pathway names
    nbPathways <- length(pathwayNames) # number of pathways

    cluster <- as.data.frame(matrix(rep(0, nbPathways*nbPathways), nrow=nbPathways, ncol=nbPathways))
    diag(cluster) <- 1
    message(paste0("Number of pathways: ", nbPathways))
    message(paste0("Number of NAs in rownames: ", sum(is.na(pathwayNames))))
    colnames(cluster) <- pathwayNames
    rownames(cluster) <- pathwayNames
  
    while(max(simmat[upper.tri(simmat, diag=FALSE)]) >= threshold){
        index <- which(simmat == max(simmat[upper.tri(simmat, diag = FALSE)]), arr.ind=TRUE)
        t_i = index[1,1] # first row index
        t_j = index[1,2] # first col index
        t_both <- c(t_i, t_j)
        pair_size <- c(length(alist[[rownames(simmat)[t_i]]]), length(alist[[rownames(simmat)[t_j]]]))
        if ( (pair_size[1] <= lowThreshold || pair_size[1] >= highThreshold) && (pair_size[2] <= lowThreshold || pair_size[2] >= highThreshold) ) {
            failed <- rownames(simmat)[t_both[which.max(pair_size)]]
            passed <- rownames(simmat)[t_both[-which.max(pair_size)]]
            simmat <- update_simmat(failed, passed, simmat, alist)
            cluster[failed,passed] <- 1
            cluster[failed,failed] <- -1
        } else if ( (length(alist[[rownames(simmat)[t_i]]]) <= lowThreshold || length(alist[[rownames(simmat)[t_i]]]) >= highThreshold) ) {
            failed <- rownames(simmat)[t_i]
            passed <- rownames(simmat)[t_j]
            simmat <- update_simmat(failed, passed, simmat, alist)
            cluster[failed,passed] <- 1
            cluster[failed,failed] <- -1
        } else if ( (length(alist[[rownames(simmat)[t_j]]]) <= lowThreshold || length(alist[[rownames(simmat)[t_j]]]) >= highThreshold) ) {
            failed <- rownames(simmat)[t_j]
            passed <- rownames(simmat)[t_i]
            simmat <- update_simmat(failed, passed, simmat, alist)
            cluster[failed,passed] <- 1
            cluster[failed,failed] <- -1
        } else {
            if ( pvalues[[rownames(simmat)[t_i]]] != pvalues[[rownames(simmat)[t_j]]] ) {
                less_sig <- which.max(c(pvalues[[rownames(simmat)[t_i]]], pvalues[[rownames(simmat)[t_j]]]))
                failed <- rownames(simmat)[t_both[less_sig]]
                passed <- rownames(simmat)[t_both[-less_sig]]
                simmat <- update_simmat(failed, passed, simmat, alist)
                cluster[failed,passed] <- 1
                cluster[failed,failed] <- -1
            } else if (adjmat[rownames(simmat)[t_i],rownames(simmat)[t_j]] != 0) {
                # reject childterm (more specific)
                child <- which.min(c(sum(alist[[rownames(simmat)[t_i]]]),sum(alist[[rownames(simmat)[t_j]]])))
                failed <- rownames(simmat)[t_both[child]]
                passed <- rownames(simmat)[t_both[-child]]
                simmat <- update_simmat(failed, passed, simmat, alist)
                cluster[failed,passed] <- 1
                cluster[failed,failed] <- -1
                # none of the above filters apply: random choice
            } else {
                failed <- rownames(simmat)[t_i]
                passed <- rownames(simmat)[t_j]
                simmat <- update_simmat(failed, passed, simmat, alist)
                cluster[failed,passed] <- 1
                cluster[failed,failed] <- -1
            }
        }
    }
    return(list("clusterMat"=cluster, "simMat"=simmat))
}


