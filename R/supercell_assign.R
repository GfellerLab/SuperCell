#' Assign super-cells to the most aboundant cluster
#'
#'
#' @param clusters a vector of clustering assignment
#' @param supercell_membership a vector of assignment of single-cell data to super-cells (membership field of \code{\link{SCimplify}} function output)
#' @param method method to define the most abuldant cell cluster within super-cells. Available: "jaccard" (default), "relative", "absolute".
#' \itemize{
#'   \item jaccard - assignes super-cell to cluster with the maximum jaccard coefficient (recommended)
#'   \item relative - assignes super-cell to cluster with the maximum relative abundance (normalized by cluster size), may result in assignment of super-cells to poorly represented (small) cluser due to normalizetaion
#'   \item absolute - assignes super-cell to cluster with the maximum absolute abundance within super-cell, may result in disappearence of poorly represented (small) clusters
#' }
#'
#' @return a vector of super-cell assignment to clusters
#'
#' @export
#'


supercell_assign <- function(clusters, supercell_membership, method = c("jaccard", "relative", "absolute")){
  cl.gr            <- table(clusters, supercell_membership)
  cluster.size     <- as.numeric(table(clusters))
  group.size       <- as.numeric(table(supercell_membership))

  if(is.null(method[1]) | is.na(method[1]) | is.nan(method[1])){
    stop(paste("Please specify method: jaccard (recommended), relative or absolute"))
  }

  if(method[1] == "jaccard"){

    cl.gr          <- as.matrix(cl.gr)
    jaccard.mtrx   <- cl.gr

    for(i in rownames(cl.gr)){
      for(j in colnames(cl.gr)){
        jaccard.mtrx[i,j] <- cl.gr[i,j] / (sum(cl.gr[i,]) +  sum(cl.gr[,j]) - cl.gr[i,j])
      }
    }
    res              <- apply(jaccard.mtrx, 2, function(x){names(x)[which.max(x)]})

  } else if(method[1] == "relative"){
    cl.gr            <- sweep(cl.gr, 1, cluster.size, "/")
    res              <- apply(cl.gr, 2, function(x){names(x)[which.max(x)]})

  } else if(method[1] == "absolute"){
    res              <- apply(cl.gr, 2, function(x){names(x)[which.max(x)]})

  } else {
    stop(paste("Unknown value of method (", method[1] , ")", "please, use: jaccard, relative or absolute"))
  }


  return(res)
}
