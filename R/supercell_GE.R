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

supercell_GE <- function(ge, groups, weights = NULL, do.median.norm = FALSE){
  if(ncol(ge) != length(groups)){
    stop("Length of the vector groups has to be equal to the number of cols in matrix ge")
  }
  N.SC <- max(groups)
  supercell_size <- as.vector(table(groups))
  j <- rep(1:N.SC, supercell_size) # column indices of matrix M.AV that, whene GE.SC <- ge %M.AV%

  goups.idx  <- plyr:::split_indices(groups)
  i <- unlist(goups.idx) # row indices of matrix M.AV that, whene GE.SC <- ge %M.AV%

  if(is.null(weights)){
    M.AV <- Matrix::sparseMatrix(i = i, j = j)
    GE.SC <- ge %*% M.AV
    GE.SC <- sweep(GE.SC, 2, supercell_size, "/")
  } else {

    if(length(weights) != length(groups)){
      stop("weights must be the same length as groups or NULL in case of unweighted averaging")
    }
    M.AV <- Matrix::sparseMatrix(i = i, j = j, x = weights[i])
    GE.SC <- ge %*% M.AV

    weighted_supercell_size <- unlist(lapply(goups.idx, FUN = function(x){sum(weights[x])}))
    GE.SC <- sweep(GE.SC, 2, weighted_supercell_size, "/")
  }

  if(do.median.norm){
    GE.SC <- (GE.SC+0.01)/apply(GE.SC+0.01, 1, median)
  }

  return(GE.SC)
}




#' Simplification of scRNA-seq dataset (old version, not used since 12.02.2021)
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


supercell_GE_idx  <- function(ge, groups, weights = NULL, do.median.norm = FALSE){

  if(ncol(ge) != length(groups)){
    stop("Length of the vector groups has to be equal to the number of cols in matrix ge")
  }

  if(ncol(ge) > 200000){
    block.size <- 5000
    N.blocks <- nrow(ge)%/%block.size
    if(nrow(ge)%%block.size > 0) N.blocks <- N.blocks+1
  } else {
    block.size <- nrow(ge)
    N.blocks <- 1
  }

  print("N.blocks:")
  print(N.blocks)

  goups.idx    <- plyr:::split_indices(groups)
  supercell.GE <- c()

  if(N.blocks>0){
    for(i in 1:N.blocks){
      print(i)

      idx.begin <- (i-1)*block.size + 1
      idx.end   <- min(i*block.size,  nrow(ge))

      ge.i      <- as.matrix(ge[idx.begin:idx.end,])

      if(is.null(weights)){
        fun <- function(idx){
          Matrix::rowMeans(ge.i[, idx, drop = FALSE])
        }
      } else {
        if(length(weights) != length(groups))
          stop("weights must be the same length as groups or NULL in case of unweighted averaging")
        fun <- function(idx){
          matrixStats::rowWeightedMeans(ge.i[, idx, drop = FALSE], w = weights[idx])
        }
      }

      supercell.GE             <- rbind(supercell.GE, sapply(goups.idx, fun))
    }
    # if(!(TRUE %in% is.na(as.numeric(colnames(supercell.GE))))){
    #   supercell.GE <- supercell.GE[,order(as.numeric(colnames(supercell.GE)))]
    # }

  }
  if(do.median.norm){
    supercell.GE <- (supercell.GE+0.01)/apply(supercell.GE+0.01, 1, median)
  }
  return(supercell.GE)
}


