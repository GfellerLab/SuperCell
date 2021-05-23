#' Rescale supercell object
#'
#' This function recomputes super-cell structure at a different graining level (\code{gamma}) or
#' for a specific number of super-cells (\code{N.SC})
#' @param SC.object super-cell object (an output from \code{\link{SCimplify}} function)
#' @param gamma new grainig level (provide either \code{gamma} or \code{N.SC})
#' @param N.SC new number of super-cells (provide either \code{gamma} or \code{N.SC})
#'
#' @return the same object as \code{\link{SCimplify}} at a new graining level
#' @export


supercell_rescale <- function(SC.object, gamma = NULL, N.SC = NULL){

  if(is.null(SC.object$h_membership) | is.null(SC.object$graph.singlecell)){
    stop("rescaling is impossible if there is no field SC.object$h_membership or SC.object$graph.singlecell,
         please make sure you run function SCimplify with parameters return.singlecell.NW = TRUE and
         return.hierarchical.structure = TRUE")
  }

  if(SC.object$do.approx){
    stop("rescaling is not yet available for an approximate simplification results")
  }

  if(!xor(is.null(gamma), is.null(N.SC))){
    stop("Provide either gamma or N.SC (just 1)")
  }

  N.c <- length(SC.object$membership)

  if(!is.null(N.SC)){
    N.SC <- round(N.c/gamma)
    if(N.SC < 2 | N.SC > N.c)
      stop("N.SC is out of range")
  } else {
    if(gamma >= 1 & gamma < (N.c/2))
      N.SC <- round(N.c/gamma)
    else
      stop("gamma is out of range")
  }

  membership       <- igraph::cut_at(SC.object$h_membership, N.SC)
  supercell_size   <- as.vector(table(membership))

  SC.NW                          <- igraph::contract(SC.object$graph.singlecell, membership)
  SC.NW                          <- igraph::simplify(SC.NW, remove.loops = T, edge.attr.comb="sum")

  igraph::E(SC.NW)$width         <- sqrt(igraph::E(SC.NW)$weight/10)
  igraph::V(SC.NW)$size          <- supercell_size
  igraph::V(SC.NW)$sizesqrt      <- sqrt(igraph::V(SC.NW)$size)


  res <- list(graph.supercells    = SC.NW,
              gamma               = gamma,
              N.SC                = length(membership),
              membership          = membership,
              supercell_size      = supercell_size,
# unchanged fields
              genes.use           = SC.object$genes.use,
              simplification.algo = SC.object$simplification.algo,
              do.approx           = SC.object$do.approx
              #h_membership        = SC.object$h_membership,
              #graph.singlecell    = SC.object$graph.singlecell
              )
}
