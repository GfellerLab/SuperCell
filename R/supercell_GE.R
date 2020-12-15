#' Simplification of scRNA-seq dataset
#'
#' This function converts gene-expression matrix of single-cell data into a gene expression
#' matrix of super-cells
#'
#' @param ge gene expression matrix (or any coordinate matrix) with genes as rows and cells as cols
#' @param groups vector of membership (assignment of single-cell to super-cells)
#' @param weights vector of a cell weight (NULL by default), used for computing average gene expression withing cluster of super-cells
#' @param do.median.norm whether to normalize by median value (FALSE by default)
#'
#' @return a matrix of simplified (averaged withing groups) data with ncol equal to number of groups and nrows as in the initial dataset
#' @export

supercell_GE  <- function(ge, groups, weights = NULL, do.median.norm = FALSE){

  if(ncol(ge) != length(groups)){
    stop("Length of the vector groups has to be equal to the number of cols in matrix ge")
  }

  ge        <- as.matrix(ge)
  goups.idx <-  plyr:::split_indices(groups)

  if(is.null(weights)){
    fun <- function(idx){
      Matrix::rowMeans(ge[, idx, drop = FALSE])
    }
  } else {
    if(length(weights) != length(groups))
      stop("weights must be the same length as groups or NULL in case of unweighted averaging")
    fun <- function(idx){
      matrixStats::rowWeightedMeans(ge[, idx, drop = FALSE], w = weights[idx])
    }
  }

  supercell.GE             <- sapply(goups.idx, fun)
  # if(!(TRUE %in% is.na(as.numeric(colnames(supercell.GE))))){
  #   supercell.GE <- supercell.GE[,order(as.numeric(colnames(supercell.GE)))]
  # }

  if(do.median.norm){
    supercell.GE <- (supercell.GE+0.01)/apply(supercell.GE+0.01, 1, median)
  }
  return(supercell.GE)
}


