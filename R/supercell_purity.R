#' Compute purity of super-cells
#'
#'
#' @param clusters vector of clustering assignment (reference assignment)
#' @param supercell_membership vector of assignment of single-cell data to super-cells (membership field of \code{\link{SCimplify}} function output)
#'
#'
#' @return a vector of super-cell purity, which is defined as a proportion of the most abundant cluster within super-cell.
#' With 1 meaning that super-cell consists of single cells from one cluster (reference assignment)
#'
#' @export
#'

supercell_purity <- function(clusters, supercell_membership){
  cl.gr            <- table(clusters, groups)
  cluster.size     <- as.numeric(table(clusters))
  group.size       <- as.numeric(table(groups))

  Ng               <- length(group.size)
  group.max.cl     <- rep(0, Ng)


  cl.gr            <- sweep(cl.gr, 2, group.size, "/")

  res              <- apply(cl.gr, 2, max)
  return(res)
}
