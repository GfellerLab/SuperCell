#' Detection of super-cells
#'
#' This function detects super-cells from single-cell gene expression matrix
#'
#'
#' @param X log-normalized gene expression matrix with rows to be genes and cols to be cells
#' @param genes.use a vector of genes used to compute PCA
#' @param genes.exclude a vector of genes to be excluded when computing PCA
#' @param n.var.genes if \code{"genes.use"} is not provided, \code{"n.var.genes"} genes with the largest variation are used
#' @param gamma graining level of data (proportion of number of single cells in the initial dataset to the number of super-cells in the final dataset)
#' @param k.knn parameter to compute single-cell kNN network
#' @param do.scale whether to scale gene expression matrix when computing PCA
#' @param n.pc number of principal components to use for construction of single-cell kNN network
#' @param fast.pca use irlba::irlba as a faster version of prcomp (one used in Seurat package)
#' @param do.approx compute approximate kNN in case of a large dataset (>50'000)
#' @param directed whether to build a directed graph
#' @param approx.N number of cells to subsample for an approximate approach
#' @param seed seed to use to subsample cells for an approximate approach
#' @param igraph.clustering clustering method to identify super-cells (available methods "walktrap" (default) and "louvain" (not recommended, gamma is ignored)).
#' @param return.singlecell.NW whether return single-cell network (which consists of approx.N if \code{"do.approx"} or all cells otherwise)
#' @param return.hierarchical.structure whether return hierarchical structure of super-cell
#'
#' @return a list with components
#' \itemize{
#'   \item graph.supercells - igraph object of a simplified network (number of nodes corresponds to number of super-cells)
#'   \item membership - assigmnent of each single cell to a particular super-cell
#'   \item graph.singlecells - igraph object (kNN network) of single-cell data
#'   \item supercell_size - size of super-cells
#' }
#'
#' @examples
#' \dontrun{
#' data(cell_lines) # list with GE - gene expression matrix (logcounts), meta - cell meta data
#' GE <- cell_lines$GE
#'
#' SC <- SCimplify(GE,  # log-normalized gene expression matrix
#'                 gamma = 20, # graining level
#'                 n.var.genes = 1000, # number of top varible genes to use for the dimensionality reduction, the list of genes can be provided instead (with 'genes.use')
#'                 k.knn = 5, # k for kNN algorithm
#'                 n.pc = 10, # number of proncipal components to use
#'                 do.approx) # whether to run an approxiate coarse-graining (for large daatasets with N_cells > 50'000)
#'
#' }
#' @export
#'
#'
#'


SCimplify <- function(X,
                      genes.use = NULL,
                      genes.exclude = NULL,
                      n.var.genes = min(1000, nrow(X)),
                      gamma = 10,
                      k.knn = 5,
                      do.scale = TRUE,
                      n.pc = 10,
                      fast.pca = TRUE,
                      do.approx = FALSE,
                      approx.N = 20000,
                      directed = FALSE,
                      use.nn2 = TRUE,
                      seed = 12345,
                      igraph.clustering = c("walktrap", "louvain"),
                      return.singlecell.NW = TRUE,
                      return.hierarchical.structure = TRUE,
                      block.size = 10000){

  N.c <- ncol(X)

  if(is.null(rownames(X))){
     if(!(is.null(genes.use) | is.null(genes.exclude))){
        stop("rownames(X) is Null \nGene expression matrix X is expected to have genes as rownames")
     } else {
       warning("colnames(X) is Null, \nGene expression matrix X is expected to have genes as rownames! \ngenes will be created automatically in a form 'gene_i' ")
       rownames(X) <- paste("gene", 1:nrow(X), sep = "_")
     }
}

  if(is.null(colnames(X))){
    warning("colnames(X) is Null, \nGene expression matrix X is expected to have cellIDs as colnames! \nCellIDs will be created automatically in a form 'cell_i' ")
    colnames(X) <- paste("cell", 1:N.c, sep = "_")
  }

  keep.genes    <- setdiff(rownames(X), genes.exclude)
  X             <- X[keep.genes,]


  if(is.null(genes.use)){
    n.var.genes <- min(n.var.genes, nrow(X))
    if(N.c > 50000){
      set.seed(seed)
      idx         <- sample(N.c, 50000)
      gene.var    <- apply(X[,idx], 1, var)
    } else {
      gene.var    <- apply(X, 1, var)
    }

    genes.use   <- names(sort(gene.var, decreasing = TRUE))[1:n.var.genes]
  }

  if(length(intersect(genes.use, genes.exclude)) > 0){
    stop("Sets of genes.use and genes.exclude have non-empty intersection")
  }

  genes.use <- genes.use[genes.use %in% rownames(X)]
  X <- X[genes.use,]

  if(do.approx & approx.N >= N.c){
    do.approx <- FALSE
    warning("N.approx is larger or equal to the number of single cells, thus, an exact simplification will be performed")
  }

  if(do.approx & (approx.N < round(N.c/gamma))){
    approx.N <- round(N.c/gamma)
    warning(paste("N.approx is set to N.SC", approx.N))
  }

  if(do.approx & ((N.c/gamma) > (approx.N/3))){
    warning("N.approx is not much larger than desired number of super-cells, so an approximate simplification may take londer than an exact one!")
  }

  if(do.approx){
    set.seed(seed)
    approx.N            <- min(approx.N, N.c)
    presample           <- sample(1:N.c, size = approx.N, replace = FALSE)
    presampled.cell.ids <- colnames(X)[sort(presample)]
    rest.cell.ids       <- setdiff(colnames(X), presampled.cell.ids)
  } else {
    presampled.cell.ids <- colnames(X)
    rest.cell.ids       <- c()
  }

  X.for.pca            <- Matrix::t(X[genes.use, presampled.cell.ids])
  if(do.scale){ X.for.pca            <- scale(X.for.pca) }
  X.for.pca[is.na(X.for.pca)] <- 0

  if(is.null(n.pc[1]) | min(n.pc) < 1){stop("Please, provide a range or a number of components to use: n.pc")}
  if(length(n.pc)==1) n.pc <- 1:n.pc

  if(fast.pca & (N.c < 1000)){
    warning("Normal pca is computed because number of cell is low for irlba::irlba()")
    fast.pca <- FALSE
  }

  if(!fast.pca){
    PCA.presampled          <- prcomp(X.for.pca, rank. = max(n.pc), scale. = F, center = F)
  } else {
    PCA.presampled          <- irlba::irlba(X.for.pca, nv = max(n.pc, 25))
    PCA.presampled$x        <- PCA.presampled$u %*% diag(PCA.presampled$d)
    PCA.presampled$rotation <- PCA.presampled$v
  }

  sc.nw <- build_knn_graph(X = PCA.presampled$x[,n.pc], k = k.knn, from = "coordinates", use.nn2 = use.nn2, dist_method = "euclidean", directed = directed)

  #simplify

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

    PCA.averaged.SC      <- as.matrix(Matrix::t(supercell_GE(t(PCA.presampled$x[,n.pc]), groups = membership.presampled)))
    X.for.roration       <- Matrix::t(X[genes.use, rest.cell.ids])



    if(do.scale){ X.for.roration <- scale(X.for.roration) }
    X.for.roration[is.na(X.for.roration)] <- 0


    membership.omitted   <- c()
    if(is.null(block.size) | is.na(block.size)) block.size <- 10000

    N.blocks <- length(rest.cell.ids)%/%block.size
    if(length(rest.cell.ids)%%block.size > 0) N.blocks <- N.blocks+1


    if(N.blocks>0){
      for(i in 1:N.blocks){ # compute knn by blocks
        idx.begin <- (i-1)*block.size + 1
        idx.end   <- min(i*block.size,  length(rest.cell.ids))

        cur.rest.cell.ids    <- rest.cell.ids[idx.begin:idx.end]

        PCA.ommited          <- X.for.roration[cur.rest.cell.ids,] %*% PCA.presampled$rotation[, n.pc] ###

        D.omitted.subsampled <- proxy::dist(PCA.ommited, PCA.averaged.SC) ###

        membership.omitted.cur        <- apply(D.omitted.subsampled, 1, which.min) ###
        names(membership.omitted.cur) <- cur.rest.cell.ids ###

        membership.omitted   <- c(membership.omitted, membership.omitted.cur)
      }
    }

    membership.all       <- c(membership.presampled, membership.omitted)
    membership.all       <- membership.all[colnames(X)]
  } else {
    membership.all       <- membership.presampled[colnames(X)]
  }

  membership       <- membership.all

  supercell_size   <- as.vector(table(membership))

  igraph::E(SC.NW)$width         <- sqrt(igraph::E(SC.NW)$weight/10)
  igraph::V(SC.NW)$size          <- supercell_size
  igraph::V(SC.NW)$sizesqrt      <- sqrt(igraph::V(SC.NW)$size)

  res <- list(graph.supercells = SC.NW,
              gamma = gamma,
              N.SC = length(membership),
              membership = membership,
              supercell_size = supercell_size,
              genes.use = genes.use,
              simplification.algo = igraph.clustering[1],
              do.approx = do.approx,
              n.pc = n.pc,
              k.knn = k.knn)

  if(return.singlecell.NW){res$graph.singlecell <- sc.nw$graph.knn}

  if(igraph.clustering[1] == "walktrap" & return.hierarchical.structure)  res$h_membership <- g.s

  return(res)
}

