#' Run RNAvelocity for super-cells (slightly modified from \link[velocyto.R]{gene.relative.velocity.estimates})
#' Not yet adjusted for super-cell size (not sample-weighted)
#'
#' @param emat spliced (exonic) count matrix
#' @param nmat unspliced (nascent) count matrix
#' @param membership supercell membership ('membership' field of \link{SCimplify})
#' @param supercell_size a vector with supercell size (if emat and nmat provided at super-cell level)
#' @param ... other parameters from \link[velocyto.R]{gene.relative.velocity.estimates}
#'
#' @return reslts of \link[velocyto.R]{gene.relative.velocity.estimates} plus
#' @export


supercell_estimate_velocity <- function(emat, nmat, smat = NULL, membership = NULL, supercell_size = NULL,
                                        do.run.avegaring = (ncol(emat) == length(membership)), kCells = 1, ...){

  N.c <- ncol(emat)

  if(!identical(dim(emat), dim(nmat))){
    stop("emat and nmat have different dimensions")
  }

  if(!is.null(smat) & !identical(dim(smat), dim(nmat))){
    stop("smat has different dimension from emat and nmat")
  }

  ## whether to run averaging of emat & nmat (if nmat provided at a single-cell level)
  if(do.run.avegaring){

    emat <- supercell_GE(ge = emat, groups = membership)
    nmat <- supercell_GE(ge = nmat, groups = membership)

    if(!is.null(smat)){
      smat <- supercell_GE(ge = smat, groups = membership)
    }

    if(!is.null(supercell_size)){
      warning("supercell_size was recomputed from membership")
    }
    supercell_size <- as.vector(table(membership))

    N.SC <- ncol(emat)
  } else { # already averaged data provided as input ?
    if(is.null(supercell_size)){
      supercell_size <- rep(1, N.c)
      warning("supercell_size was replaced with 1")
    } else if(length(supercell_size) != N.c) {
      stop("supercell_size has different length from the number of super-cells")
    }
    N.SC <- ncol(emat)
  }

  rel.velocity <- velocyto.R::gene.relative.velocity.estimates(emat = emat, nmat = nmat, smat = smat,
                                                               kCells = kCells, ...)

  rel.velocity$is_weighted <- FALSE ## for the moment, the weighted version of rna-velocity is under development
  rel.velocity$supercell_size <- supercell_size
  return(rel.velocity)

}

