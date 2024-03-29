% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_knn_graph.R
\name{build_knn_graph_nn2}
\alias{build_knn_graph_nn2}
\title{Build kNN graph using RANN::nn2
(used in \code{"build_knn_graph"})}
\usage{
build_knn_graph_nn2(
  X,
  k = min(5, ncol(X)),
  mode = "all",
  DoSNN = FALSE,
  which.snn = c("bluster", "dbscan"),
  pruning = NULL,
  kmin = 0,
  ...
)
}
\arguments{
\item{X}{matrix of coordinates (rows are samples and cols are coordinates)}

\item{k}{kNN parameter}

\item{mode}{mode of \link[igraph]{graph_from_adj_list} ('all' -- undirected graph, 'out' -- directed graph)}

\item{DoSNN}{whether to apply shared nearest neighbors (default is \code{FALSE})}

\item{which.snn}{whether to use \link[bluster]{neighborsToSNNGraph} or \link[dbscan]{sNN} for sNN graph construction}

\item{pruning}{quantile to perform edge pruning (default is \code{NULL} - no pruning applied) based on PCA distance distribution}

\item{kmin}{keep at least \code{kmin} edges in single-cell graph when pruning applied (idnored if \code{is.null(pruning)})}

\item{...}{other parameters of \link[bluster]{neighborsToSNNGraph} or \link[dbscan]{sNN}}
}
\value{
a list with components
\itemize{
  \item graph.knn - igraph object
}
}
\description{
Build kNN graph using RANN::nn2
(used in \code{"build_knn_graph"})
}
