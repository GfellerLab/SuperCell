#' Simplification of scRNA-seq dataset
#'
#' This function converts gene-expression matrix of single-cell data into a gene expression
#' matrix of super-cells
#'
#' @param ge gene expression matrix (or any coordinate matrix) with genes as rows and cells as cols
#' @param groups vector of membership (assignment of single-cell to super-cells)
#' @param do.median.norm whether to normalize
#'
#' @return a matrix of simplified (averaged withing groups) data with ncol equal to number of groups and nrows as in the initial dataset
#' @export

supercell_GE  <- function(ge, groups, do.median.norm = FALSE){
  goups.idx <-  plyr:::split_indices(groups)

  fun <- function(idx){
    Matrix::rowMeans(ge[, idx, drop = FALSE])
  }
  supercell.GE             <- sapply(goups.idx, fun)

  if(!(TRUE %in% is.na(as.numeric(colnames(supercell.GE))))){
    supercell.GE <- supercell.GE[,order(as.numeric(colnames(supercell.GE)))]
  }

  if(do.median.norm){
    supercell.GE <- (supercell.GE+0.01)/apply(supercell.GE+0.01, 1, median)
  }
  return(supercell.GE)
}


###### rest of the functions were used to identify the fastest way to average gene expression within super-cells

supercell_GE_loop  <- function(ge, groups, do.median.norm = FALSE){
  u.groups        <- unique(groups)
  u.groups        <- u.groups[!is.na(u.groups)]
  N.groups        <- length(u.groups)
  N.genes         <- nrow(ge)
  supercell.GE    <- matrix(0, nrow = N.genes, ncol = N.groups)
  rownames(supercell.GE) <- rownames(ge)
  colnames(supercell.GE) <- u.groups
  for(g in as.character(u.groups)){
    idxs             <- which(groups == g)
    #print(idxs)
    if(length(idxs) < 2){
      supercell.GE[,g] <- ge[,idxs]
    } else{
      if(N.genes > 1){
        supercell.GE[,g] <- apply(ge[,idxs], 1, mean)
      } else {
        supercell.GE[,g] <- mean(ge[,idxs])
      }
    }
  }
  if(!(TRUE %in% is.na(as.numeric(colnames(supercell.GE))))){
    supercell.GE <- supercell.GE[,order(as.numeric(colnames(supercell.GE)))]
  }

  if(do.median.norm){
    supercell.GE <- (supercell.GE+0.01)/apply(supercell.GE+0.01, 1, median)
  }
  return(supercell.GE)
}



supercell_GE_sapply  <- function(ge, groups, do.median.norm = FALSE){
  goups.idx <-  plyr:::split_indices(groups)

  fun <- function(idx){
      Matrix::rowMeans(ge[, idx, drop = FALSE])
  }
  supercell.GE             <- sapply(goups.idx, fun)

  if(do.median.norm){
    supercell.GE <- (supercell.GE+0.01)/apply(supercell.GE+0.01, 1, median)
  }
  return(supercell.GE)
}



supercell_GE_aggr  <- function(ge, groups, do.median.norm = FALSE){
  u.groups        <- unique(groups)
  u.groups        <- u.groups[!is.na(u.groups)]
  N.groups        <- length(u.groups)
  N.genes         <- nrow(ge)

  supercell.GE             <- aggregate(ge, list(groups), mean)

  if(do.median.norm){
    supercell.GE <- (supercell.GE+0.01)/apply(supercell.GE+0.01, 1, median)
  }
  return(supercell.GE)
}


