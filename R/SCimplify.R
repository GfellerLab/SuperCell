#' Detection of super-cells
#'
#' This function detects super-cells from single-cell gene expression matrix
#'
#'
#' @param X gene expression matrix with rows to be genes and cols to be cells
#' @param genes.use a list of genes used to compute PCA
#' @param n.var.genes if \code{"genes.use"} is not provided, \code{"n.var.genes"} genes with the largest variation are used
#' @param k.knn parameter to compute single-cell kNN network
#' @param gamma graining level of data (proportion of number of single cells in the initial dataset to the number of super-cells in the final dataset)
#' @param do.scale whether to scale gene expression matrix to compure PCA
#' @param n.pc number of principal components to use for construction of single-cell kNN network
#' @param do.approx compute approximate kNN in case of a large dataset (>20000)
#' @param approx.N number of cells to subsample for an approximate approach
#' @param seed seed to use to subsample cells for an approximate approach
#' @param igraph.clustering clustering method to identify super-cells (available methods "walktrap" (default) and "louvain" (not recommended, gamma is ignored)).
#' @param return.singlecell.NW whether return single-cell network (which consists of approx.N if \code{"do.approx"} or all cells otherwise)
#'
#' @return a list with components
#' \itemize{
#'   \item graph.supercells - igraph object of a simplified network (number of nodes corresponds to number of super-cells)
#'   \item membership - assigmnent of each single cell to a particular super-cell
#'   \item graph.singlecells - igraph object (kNN network) of single-cell data
#' }
#' @export
#'
#'
#'


SCimplify <- function(X,
                      genes.use = NULL,
                      n.var.genes = min(1000, nrow(X)),
                      k.knn = 5,
                      gamma = 10,
                      do.scale = TRUE,
                      n.pc = 10,
                      do.approx = FALSE,
                      approx.N = 20000,
                      seed = 12345,
                      igraph.clustering = c("walktrap", "louvain"),
                      return.singlecell.NW = TRUE){

  if(is.null(genes.use)){
    n.var.genes <- min(n.var.genes, nrow(X))
    gene.var    <- apply(X, 1, var)
    genes.use   <- names(sort(gene.var, decreasing = TRUE))[1:n.var.genes]
  }

  if(do.approx){
    set.seed(seed)
    approx.N            <- min(approx.N, ncol(X))
    presample           <- sample(1:ncol(X), size = approx.N, replace = FALSE)
    presampled.cell.ids <- colnames(X)[sort(presample)]
    rest.cell.ids       <- setdiff(colnames(X), presampled.cell.ids)
  } else {
    presampled.cell.ids <- colnames(X)
    rest.cell.ids       <- c()
  }

  X.for.pca            <- t(X[genes.use, presampled.cell.ids])
  if(do.scale){ X.for.pca            <- scale(X.for.pca) }

  PCA.presampled        <- prcomp(X.for.pca, rank. = n.pc, scale. = F, center = F)

  sc.nw <- build_knn_graph(X = PCA.presampled$x, k = k.knn, from = "coordinates", use.nn2 = TRUE, dist_method = "euclidean")

  #simplify
  N.c <- ncol(X)
  k   <- round(N.c/gamma)

  if(igraph.clustering[1] == "walktrap"){
    g.s              <- igraph::cluster_walktrap(sc.nw$graph.knn)
    g.s$membership   <- igraph::cut_at(g.s, k)

  } else if(igraph.clustering[1] == "louvain") {
    warning(paste("igraph.clustering =", igraph.clustering, ", gamma is ignored"))
    g.s    <- igraph::cluster_louvain(sc.nw$graph.knn)

  } else {
    stop(paste("Unknown clustering method (", igraph.clustering, "), please use louvain or walkrtap"))
  }

  membership.presampled        <- g.s$membership
  names(membership.presampled) <- presampled.cell.ids

  SC.NW                        <- igraph::contract(sc.nw$graph.knn, membership.presampled)
  SC.NW                        <- igraph::simplify(SC.NW, remove.loops = T, edge.attr.comb="sum")

  if(do.approx){

    PCA.averaged.SC      <- t(supercell_GE(t(PCA.presampled$x), groups = membership.presampled))

    X.for.roration       <- t(X[genes.use, rest.cell.ids])
    if(do.scale){ X.for.roration <- scale(X.for.roration) }

    PCA.ommited          <- X.for.roration %*% PCA.presampled$rotation[, 1:n.pc]

    D.omitted.subsampled <- proxy::dist(PCA.ommited, PCA.averaged.SC)

    membership.omitted   <- apply(D.omitted.subsampled, 1, which.min)
    names(membership.omitted) <- rest.cell.ids

    membership.all       <- c(membership.presampled, membership.omitted)
    membership.all       <- membership.all[colnames(X)]
  } else {
    membership.all       <- membership.presampled[colnames(X)]
  }

  membership       <- membership.all

  supercell_size   <- as.vector(table(membership))

  E(SC.NW)$width         <- E(SC.NW)$weight/10
  V(SC.NW)$size          <- supercell_size
  V(SC.NW)$sizesqrt      <- sqrt(V(SC.NW)$size)

  res <- list(graph.supercells = SC.NW,
              membership = membership,
              supercell_size = supercell_size,
              genes.use = genes.use)

  if(return.singlecell.NW){res$graph.singlecell <- sc.nw$graph.knn}
  return(res)
}

