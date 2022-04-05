#' Merging independent SuperCell objects
#'
#' This function merges independent SuperCell objects
#'
#' @param SCs list of SuperCell objects (resuls of \link[SCimplify])
#'
#' @return a list with components
#' \itemize{
#'   \item membership - assigmnent of each single cell to a particular metacell
#'   \item cell.ids - the original ids of single-cells
#'   \item supercell_size - size of metacells (former super-cells)
#'   \item gamma -  graining level of the merged object (astimated as an average size of metacells as the independent SuperCell objects might have different graining levels)
#'   \item N.SC - number of obtained metacells
#'}
#'
#' @export

supercell_merge <- function(
  SCs
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

  return(res)
}


#' Merging metacell gene expression matrices from several independent SuperCell objects
#'
#' This function merges independent SuperCell objects
#'
#' @param SC.GEs list of metacell gene expression matricies (result of \link[supercell_GE]), make sure the order of the gene expression metricies is the same as in the call of \link[supercell_merge]
#'
#' @return a merged matrix of gene expression
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
