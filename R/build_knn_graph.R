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
  directed = TRUE
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
      return(build_knn_graph_nn2(X = X, k = k, mode = mode))
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

build_knn_graph_nn2 <- function(X, k = min(5, ncol(X)), mode = 'all'){
  nn2.res <- RANN::nn2(data = X, k = k)
  nn2.res <- nn2.res$nn.idx

  adj.knn       <- split(nn2.res, rep(1:nrow(nn2.res), times = ncol(nn2.res))) # get adj list

  graph.knn     <- igraph::graph_from_adj_list(adj.knn,  duplicate = F, mode = mode)

  graph.knn     <- igraph::simplify(graph.knn, remove.multiple = T)
  igraph::E(graph.knn)$weight <- 1

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




