% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/supercell_UMAP.R
\name{supercell_UMAP}
\alias{supercell_UMAP}
\title{Compute UMAP of super-cells}
\usage{
supercell_UMAP(SC, PCA_name = "SC_PCA", n.comp = NULL, n_neighbors = 15, ...)
}
\arguments{
\item{SC}{super-cell structure (output of \link{SCimplify}) with a field \code{PCA_name} containing PCA result}

\item{PCA_name}{name of \code{SC} field containing result of \link{supercell_prcomp}}

\item{n.comp}{number of vector of principal components to use for computing UMAP}

\item{n_neighbors}{number of neighbors (parameter of  \link[umap]{umap})}

\item{...}{other parameters of \link[umap]{umap}}
}
\value{
\link[umap]{umap} result
}
\description{
Computes UMAP of super-cells
}
