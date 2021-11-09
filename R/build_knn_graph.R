#' Build kNN graph
#'
#' Build kNN graph either from distance (from == "dist") or from coordinates (from == "coordinates")
#'
#' @param X either distance or matrix of coordinates (rows are samples and cols are coordinates)
#' @param k kNN parameter
#' @param from from which data type to build kNN network: "dist" if X is a distance (dissimilarity) or "coordinates" if X is a matrix with coordinates as cols and cells as rows
#' @param use.nn2 whether use RANN::nn2 method to buid kNN network faster (avaivable only for "coordinates" option)
#' @param return_neighbors_order whether return order of neighbors (not available for nn2 option)
#' @param dist_method method to compute dist (if X is a matrix of coordinates) available: c("cor", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
#' @param cor_method if distance is computed as correlation (dist_method == "cor), which type of correlation to use (available: "pearson", "kendall", "spearman")
#' @param p p param in \code{"dist"} function
#' @param directed whether to build a directed graph
#'
#' @return a list with components
#' \itemize{
#'   \item graph.knn - igraph object
#'   \item order - Nxk matrix with indecies of k nearest neighbors ordered by relevance (from 1st to k-th)
#' }
#'
#' @export
#'

build_knn_graph <- function(
  X,
  k = 5,
  from = c("dist", "coordinates"),
  use.nn2 = TRUE,
  return_neighbors_order = F,
  dist_method = "euclidean",
  cor_method = "pearson",
  p = 2,
  directed = FALSE,
  DoSNN = FALSE,
  which.snn = c("bluster", "dbscan"),
  pruning = NULL,
  kmin = 0,
  ...
  )
{
  av.methods <- c("dist", "coordinates")
  method <-  pmatch(from[1], av.methods)
  if(is.na(method)){
    stop(paste("Unlnown method", from, "Available methods are", paste(av.methods, collapse = ", ")))
  }


  if (method == 2){ # from coordinates
    if(use.nn2){
      if(dist_method != "euclidean"){
        stop(paste0("Fast nn2 function from RANN package is used, so",
                    dist_method, "distnce is not acceptable.
                  To use nn2 method, please, choose eucleadian distance.
                  If you want to use", dist_method, "distance, please set parameter use.nn2 to FALSE"))}
      mode <- ifelse(directed, 'out', 'all')
      return(build_knn_graph_nn2(X = X, k = k, mode = mode, DoSNN = DoSNN, pruning = pruning, which.snn = which.snn, kmin = kmin, ...))
    } else {
      av.dist      <- c("cor", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
      dist_method_ <-  pmatch(dist_method, av.dist)
      if(is.na(dist_method_)){
        stop(paste("Unknown distance method:", dist_method, "Available dist methods are", paste(av.dist, collapse = ", ") ))
      }
      if(dist_method_ == 1){
        #print("cor")
        #print(cor_method)
        av.cor_methods <- c("pearson", "kendall", "spearman")
        cor_method_    <- pmatch(cor_method, av.cor_methods)
        if(is.na(cor_method_)){
          stop(paste("Unknown cor method:", cor_method, "Available cor methods are", paste(av.cor_methods)))
        }
        X <- as.dist(as.matrix(1 - cor(t(X), method = cor_method)))
      } else {
        X <- dist(X, method = dist_method)
      }
    }
  } else {
    if(use.nn2 == TRUE){
      stop("Method nn2 cannot be applied to distance, to use fast nn2 method, please provide coordinates rather than distance
             and set parameter from to coordinates")
    }
    return(knn_graph_from_dist(D = X, k = k, return_neighbors_order = return_neighbors_order))
  }


  ### now X is distance in any case
  return(knn_graph_from_dist(D = X, k = k, return_neighbors_order = return_neighbors_order))

}


#' Build kNN graph using RANN::nn2
#' (used in \code{"build_knn_graph"})
#'
#' @param X matrix of coordinates (rows are samples and cols are coordinates)
#' @param k kNN parameter
#' @param return_neighbors_order whether return order of neighbors
#' @param mode mode of \link[igraph]{graph_from_adj_list} ('all' -- undirected graph, 'out' -- directed graph)
#'

#' @return a list with components
#' \itemize{
#'   \item graph.knn - igraph object
#' }
#'

build_knn_graph_nn2 <- function(
  X,
  k = min(5, ncol(X)),
  mode = 'all',
  DoSNN = FALSE,
  which.snn = c("bluster", "dbscan"),
  pruning = NULL,
  kmin = 0,
  ...
){


  if(is.numeric(pruning)){
    if(pruning < 0 | pruning > 1){
      stop(paste("Pruning has to be numeric from 0 to 1 or NULL"))
    }
  } else {
    if(!is.null(pruning)){
      warning("Pruning has to be numeric from 0 to 1 or NULL. pruning was set to NULL")
    }
    pruning <- NULL
  }


  if(DoSNN){

    which.snn <- which.snn[1]
    if(!(which.snn %in% c("bluster", "dbscan"))){
      stop(paste("Which SNN:", which.snn, "is unknown, available parameters are:", paste(c("bluster", "dbscan"), collapse = ", ")))
    }

    if(which.snn == "bluster"){

      nn2.res <- RANN::nn2(data = X, k = k) # it is faster and memory more efficietn to run RANN::nn2() + bluster::neighborsToSNNGraph() instead of bluster::makeSNNGraph()
      snn <- bluster::neighborsToSNNGraph(nn2.res$nn.idx, ...)

      if(!is.null(pruning)){
        pca_edje_dist      <- apply(igraph::get.edgelist(snn), 1, function(i){dist(X[i,])})
        igraph::E(snn)$pca_edje_dist <- pca_edje_dist
        pruning_cutoff     <- quantile(pca_edje_dist, pruning)
        edges_to_remove    <- which(pca_edje_dist > pruning_cutoff)
        snn <- igraph::delete_edges(snn, edges_to_remove)
      }

      graph.knn     <- igraph::simplify(snn, remove.multiple = T)
      #igraph::E(snn)$weight <- 1
      return(res <- list(graph.knn = graph.knn))

    } else if (which.snn == "dbscan"){

      kt          <- max(1, round(k/3))

      dbsnn       <- dbscan::sNN(x = X, k = 2*(k-1), ...) # does not create self-loop, so k = k-1

      if(!is.null(pruning)){
        pca_edje_dist      <- as.vector(dbsnn$dist)
        pruning_cutoff     <- quantile(pca_edje_dist, pruning)

        remove_edges       <- dbsnn$dist > pruning_cutoff
        # keep at least (kmin) edges
        if(kmin > 0){
          remove_edges[,1:kmin] <- FALSE
        }

        nn.idx             <- dbsnn$id

        dbsnn$id[remove_edges]     <- NA
        dbsnn$dist[remove_edges]   <- NA
        dbsnn$shared[remove_edges] <- NA

      }

      adj.knn       <- split(dbsnn$id, rep(1:nrow(dbsnn$id), times = ncol(dbsnn$id))) # get adj list
      adj.knn       <- lapply(adj.knn, function(x){x[!is.na(x)]}) # remove NA

      edge.dist   <- split(dbsnn$dist, rep(1:nrow(dbsnn$dist), times = ncol(dbsnn$dist))) # get edje weight in a format of adj list
      edge.dist   <- lapply(edge.dist, function(x){x[!is.na(x)]}) # remove NA

      edge.weight   <- split(dbsnn$shared, rep(1:nrow(dbsnn$shared), times = ncol(dbsnn$shared))) # get edje weight in a format of adj list
      edge.weight   <- lapply(edge.weight, function(x){x[!is.na(x)]}) # remove NA


      graph.knn                             <- igraph::graph_from_adj_list(adj.knn, duplicate = F, mode = mode)
      print(paste("N edges =", igraph::ecount(graph.knn)))
      print(paste("Length of edge.weight =", length(unlist(edge.weight))))
      print(paste("Length of edge.weight =", length(unlist(edge.weight))))

      igraph::E(graph.knn)$weight           <- unlist(edge.weight)
      igraph::E(graph.knn)$pca_edje_dist    <- unlist(edge.dist)

      return(res <- list(graph.knn = graph.knn))

    } else {
      stop(paste("Which SNN:", which.snn, "is unknown, available parameters are:", paste(c("bluster", "dbscan"), collapse = ", ")))
    }
  }


  # run just kNN
  nn2.res <- RANN::nn2(data = X, k = k)
  if(!is.null(pruning)){

    nn2.res$nn.dists <- nn2.res$nn.dists[,-1] #remove self loops
    nn2.res$nn.idx   <- nn2.res$nn.idx[,-1]

    pca_edje_dist      <- as.vector(nn2.res$nn.dists)
    pruning_cutoff     <- quantile(pca_edje_dist, pruning)

    remove_edges       <- nn2.res$nn.dists > pruning_cutoff
    # keep at least (kmin) edges
    if(kmin > 0){
      remove_edges[,1:kmin] <- FALSE
    }

    nn.idx                   <- nn2.res$nn.idx
    nn.idx[remove_edges]     <- NA

  } else {
    nn.idx <- nn2.res$nn.idx
  }

  adj.knn       <- split(nn.idx, rep(1:nrow(nn.idx), times = ncol(nn.idx))) # get adj list
  adj.knn       <- lapply(adj.knn, function(x){x[!is.na(x)]}) # remove NA

  graph.knn     <- igraph::graph_from_adj_list(adj.knn, duplicate = F, mode = mode)

  graph.knn     <- igraph::simplify(graph.knn, remove.multiple = T)
  igraph::E(graph.knn)$weight <- 1

  pca_edje_dist      <- apply(igraph::get.edgelist(graph.knn), 1, function(i){dist(X[i,])})
  igraph::E(graph.knn)$pca_edje_dist <- pca_edje_dist

  return(res <- list(graph.knn = graph.knn))

}


#' Build kNN graph from distance
#' (used in \code{"build_knn_graph"})
#'
#' @param X dist matrix or dist object (preferentially)
#' @param k kNN parameter
#' @param mode mode of \link[igraph]{graph_from_adj_list} ('all' -- undirected graph, 'out' -- directed graph)
#'
#' @return a list with components
#' \itemize{
#'   \item graph.knn - igraph object
#'   \item order - Nxk matrix with indecies of k nearest neighbors ordered by relevance (from 1st to k-th)
#' }
#'

knn_graph_from_dist <- function(D, k = 5, return_neighbors_order = T, mode = 'all'){

  ##print("Start knn_graph_from_dist")
  if(!is.matrix(D) & class(D) != "dist"){
    stop("D (X) mast be a matrix or dist!")
  }

  if(class(D) != "dist"){
    D <- as.dist(D)
  }



  N        <- (1 + sqrt(1+8*length(D)))/2 # number of cells

  if (k >= N)
    stop("Not enought neighbors in data set!")
  if (k < 1)
    stop("Invalid number of nearest neighbors, k must be >= 1!")

  row <- function(i, N){
    return(c(if(i>1) D[(i-1)+c(0:(i-2))*(N - 1 - c(1:(i-1))/2)],
             NA,
             if(i < N) D[((i-1)*(N-1) - ((i-1)*(i-2)/2) + 1) : (((i-1)*(N-1) - ((i-1)*(i-2)/2) + 1) + N-i-1)]))
  }
  neighbors <- t(sapply(1:N, function(i) {order(row(i,N))[1:k]}))

  adj.knn <- split(neighbors, rep(1:nrow(neighbors), times = ncol(neighbors)))


  graph.knn     <- igraph::graph_from_adj_list(adj.knn,  duplicate = F, mode = mode)
  graph.knn     <- igraph::simplify(graph.knn, remove.multiple = T)
  igraph::E(graph.knn)$weight <- 1

  if(return_neighbors_order){
    res <- list(graph.knn = graph.knn,
                order = neighbors)
  } else {res <- list(graph.knn = graph.knn)}

  return(res)
}




