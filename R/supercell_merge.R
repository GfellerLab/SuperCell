#' Merging independent SuperCell objects
#'
#' This function merges independent SuperCell objects
#'
#' @param SCs list of SuperCell objects (resuls of \link[SCimplify])
#' @param fields which additional fields (e.g., metadata) of the the SuperCell objects to keep when merging
#'
#' @return a list with components
#' \itemize{
#'   \item membership - assigmnent of each single cell to a particular metacell
#'   \item cell.ids - the original ids of single-cells
#'   \item supercell_size - size of metacells (former super-cells)
#'   \item gamma -  graining level of the merged object (astimated as an average size of metacells as the independent SuperCell objects might have different graining levels)
#'   \item N.SC - number of obtained metacells
#' }
#'
#'
#' @examples
#' \dontrun{
#' data(cell_lines) # list with GE - gene expression matrix (logcounts), meta - cell meta data
#' GE <- cell_lines$GE
#' cell.meta <- cell_lines$meta
#'
#' cell.idx.HCC827 <- which(cell.meta == "HCC827")
#' cell.idx.H838   <- which(cell.meta == "H838")
#'
#' SC.HCC827 <- SCimplify(GE[,cell.idx.HCC827],  # log-normalized gene expression matrix
#'                 gamma = 20, # graining level
#'                 n.var.genes = 1000, # number of top varible genes to use for the dimensionality reduction, the list of genes can be provided instead (with 'genes.use')
#'                 k.knn = 5, # k for kNN algorithm
#'                 n.pc = 10) # number of proncipal components to use
#' SC.HCC827$cell.line <- supercell_assign(cell.meta[cell.idx.HCC827], supercell_membership = SC.HCC827$membership)
#'
#' SC.H838 <- SCimplify(GE[,cell.idx.H838],  # log-normalized gene expression matrix
#'                 gamma = 30, # graining level
#'                 n.var.genes = 1000, # number of top varible genes to use for the dimensionality reduction, the list of genes can be provided instead (with 'genes.use')
#'                 k.knn = 5, # k for kNN algorithm
#'                 n.pc = 15) # number of proncipal components to use
#' SC.H838$cell.line <- supercell_assign(cell.meta[cell.idx.H838], supercell_membership = SC.H838$membership)
#'
#' SC.merged <- supercell_merge(list(SC.HCC827, SC.H838), fields = c("cell.line"))
#'
#' # compute metacell gene expression for SC.HCC827
#' SC.GE.HCC827 <- supercell_GE(GE[, cell.idx.HCC827], groups = SC.HCC827$membership)
#' # compute metacell gene expression for SC.H838
#' SC.GE.H838 <- supercell_GE(GE[, cell.idx.H838], groups = SC.H838$membership)
#' # merge GE matricies
#' SC.GE.merged <- supercell_mergeGE(list(SC.GE.HCC827, SC.GE.H838))
#'
#' }
#' @export

supercell_merge <- function(
  SCs,
  fields = c()
){
  N <- length(SCs)

  # update membership vector
  membership <- c()
  max.membership <- 0

  supercell_size <- c()
  gamma <- c()

  for(i in 1:N){
    cur.SC <- SCs[[i]]

    cur.membership <- cur.SC$membership
    cur.membership <- cur.membership + max.membership
    membership     <- c(membership, cur.membership)
    max.membership <- max(membership)

    supercell_size <- c(supercell_size, as.vector(table(cur.membership)))

    gamma          <- c(gamma, cur.SC$gamma)
  }


  res <- list(
    membership        = membership,
    cell.ids          = names(membership),
    supercell_size    = supercell_size,
    N.SC              = length(unique(membership)),
    gamma             = round(mean(supercell_size))
  )

  meta.fields <- list()
  for(f in fields){
    meta.fields[[f]] <- c()
    for(i in 1:N){

      cur.SC <- SCs[[i]]

      if(!(f %in% names(cur.SC))){
        warning(paste0("Field '", f, "' was not found in a SuperCell object, replaced with NAs \n"))
        cur.meta.f <- rep(NA, cur.SC$N.SC)
      } else {
        cur.meta.f <- cur.SC[[f]]
      }

      meta.fields[[f]] <- c(meta.fields[[f]], cur.meta.f)
    }
  }

  res <- c(res, meta.fields)

  return(res)
}


#' Merging metacell gene expression matrices from several independent SuperCell objects
#'
#' This function merges independent SuperCell objects
#'
#' @param SC.GEs list of metacell gene expression matricies (result of \link[supercell_GE]), make sure the order of the gene expression metricies is the same as in the call of \link[supercell_merge]
#'
#' @return a merged matrix of gene expression
#'
#' @examples see examples in \link[supercell_merge]
#'
#' @export

supercell_mergeGE <- function(
  SC.GEs
){

  N <- length(SC.GEs)
  common.genes <- rownames(SC.GEs[[1]])

  for(i in 2:N){
    cur.GE <- SC.GEs[[i]]
    common.genes <- intersect(common.genes, rownames(cur.GE))
  }

  cat(paste(length(common.genes), "of common genes has been detected"))

  res <- SC.GEs[[1]][common.genes,]
  for(i in 2:N){
    res <- cbind(res, SC.GEs[[i]][common.genes,])
  }

  return(res)

}
