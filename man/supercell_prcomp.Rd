% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/supercell_prcomp.R
\name{supercell_prcomp}
\alias{supercell_prcomp}
\title{compute PCA for super-cell data (sample-weighted data)}
\usage{
supercell_prcomp(
  X,
  supercell_size = NULL,
  k = 20,
  do.scale = TRUE,
  do.center = TRUE,
  seed = 12345
)
}
\arguments{
\item{X}{an averaged gene expression matrix (! where rows represent super-cells and cols represent genes)}

\item{supercell_size}{a vector with supercell size (ordered the same way as in X)}

\item{k}{number of components to compute}

\item{seed}{a seed to use for \code{set.seed}}
}
\value{
the same object as \link[stats]{prcomp} result
}
\description{
compute PCA for super-cell data (sample-weighted data)
}