#' Cluster super-cell data
#'
#' @param D  a dissimilarity matrix or a dist object
#' @param k  number of clusters
#' @param supercell_size  a vector with supercell size (ordered the same way as in D)
#' @param algorithm  which algorithm to use to compute clustering: \code{"hclust"} (default) or \code{"PAM"} (see \link[WeightedCluster]{wcKMedoids})
#' @param method  which method of algorithm to use: \itemize{
#'   \item for \code{"hclust"}: "ward.D", "ward.D2" (default), "single", "complete", "average", "mcquitty", "median" or "centroid", (see \link[stats]{hclust})
#'   \item for \code{"PAM"}: "KMedoids", "PAM" or "PAMonce" (default), (see \link[WeightedCluster]{wcKMedoids})
#' }
#' @param return.hcl whether to return a result of \code{"hclust"} (only for \code{"hclust"} algorithm)
#'
#' @return a list with components
#' \itemize{
#'   \item clustering - vector of clustering assignment of super-cells
#'   \item algo - the algorithm used
#'   \item method - method used with an algorithm
#'   \item hlc - \link[stats]{hclust} result (only for \code{"hclust"} algorithm when \code{return.hcl} is TRUE)
#' }
#' @export




supercell_cluster <- function(D, k = 5, supercell_size = NULL, algorithm = c("hclust", "PAM"), method = NULL, return.hcl = TRUE){
  res <- list()

  if(!(inherits(D, "dist") | inherits(D, "matrix"))){
    stop("D must be a dist object or distance matrix")
  }

  if(is.null(algorithm[1]) | is.na(algorithm[1]) | is.nan(algorithm[1])){
    stop(paste("Please specify algorithm: PAM or hclust"))
  }


  if(algorithm[1] == "PAM"){
    if(is.null(method)) method <- "PAMonce"

    res$clustering <- WeightedCluster::wcKMedoids(diss = D, k = k, weights = supercell_size,  cluster.only = TRUE)
    cl.names       <- 1:k
    names(cl.names)<- sort(unique(res$clustering))
    res$clustering <- unname(cl.names[as.character(res$clustering)])

    res$algo       <- "PAM"
    res$method     <- method

  } else if(algorithm[1] == "hclust"){
    if(is.null(method)) method <- "ward.D2"

    hcl <- stats::hclust(d = D, method = method, members = supercell_size)

    res$clustering <- stats::cutree(hcl, k = k)
    res$algo       <- "hclust"
    res$method     <- hcl$method
    res$hcl        <- hcl
  } else {
    stop(paste("Unknown alhorithm:", algorithm[1], ", please use PAM or hclust!"))
  }

  return(res)
}
