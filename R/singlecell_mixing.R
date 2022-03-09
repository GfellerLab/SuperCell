#' Compute mixing of single-cells within supercell
#'
#' @param SC super-cell object (output of \code{\link{SCimplify}} function)
#' @param clusters vector of clustering assignment (reference assignment)
#'
#'
#' @return a vector of single-cell mixing within super-cell it belongs to, which is defined as:
#' 1 - proportion of cells of the same annotation (e.g., cell type) within the same super-cell
#' With 0 meaning that super-cell consists of single cells from one cluster (reference assignment) and higher values correspond to higher cell type mixing within super-cell
#'
#' @export


sc_mixing_score <- function(
  SC,
  clusters
){


  if("membership" %in% names(SC)){
    membership <- SC$membership
  } else {
    membership <- 1:length(clusters)
  }

  if("cells.use.idx" %in% names(SC)){
    cells.use.idx <- SC[["cells.use.idx"]]
  } else {
    cells.use.idx <- 1:length(membership)
  }

  membership <- membership[cells.use.idx]

  mixing <- t(table(membership, clusters))
  mixing <- sweep(mixing, 2, table(membership), "/")

  N.c    <- length(membership)

  #sc.mixing <- rep(NA, N.c)
  #for(i in 1:N.c){
  #  sc.mixing[i] <- 1 - mixing[clusters[i], membership[i]]
  #}

  sc.mixing <- sapply(1:N.c, function(i){1- mixing[clusters[i], membership[i]]})
  return(sc.mixing)
}
