#' Compute purity of super-cells
#'
#'
#' @param clusters vector of clustering assignment (reference assignment)
#' @param supercell_membership vector of assignment of single-cell data to super-cells (membership field of \code{\link{SCimplify}} function output)
#' @param method method to compute super-cell purity.
#' \code{"max_proportion"} if the purity is defined as a proportion of the most abundant cluster (cell type) within super-cell or
#' \code{"entropy"} if the purity is defined as the Shanon entropy of the cell types super-cell consists of.
#'
#' @return a vector of super-cell purity, which is defined as:
#' - proportion of the most abundant cluster within super-cell for \code{method = "max_proportion"} or
#' - Shanon entropy for \code{method = "entropy"}.
#' With 1 meaning that super-cell consists of single cells from one cluster (reference assignment)
#'
#' @export
#'

supercell_purity <- function(
  clusters,
  supercell_membership,
  method = c("max_proportion", "entropy")[1]
){

  if(!(method %in% c("max_proportion", "entropy"))){
    stop(paste("Method", method, "is not known. The available methods are:", paste(method, collapse = ",")))
  }

  cl.gr            <- table(clusters, supercell_membership)

  switch(method,

         entropy = {
           res <- apply(cl.gr, 2, entropy::entropy)
         },

         max_proportion = {
           cluster.size     <- as.numeric(table(clusters))
           group.size       <- as.numeric(table(supercell_membership))

           Ng               <- length(group.size)
           group.max.cl     <- rep(0, Ng)


           cl.gr            <- sweep(cl.gr, 2, group.size, "/")

           res              <- apply(cl.gr, 2, max)
         },

         {
           stop(paste("Method", method, "is not known. The available methods are:", paste(method, collapse = ",")))
         }
  )


  return(res)
}
