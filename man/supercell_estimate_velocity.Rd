% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/supercell_estimate_velocity.R
\name{supercell_estimate_velocity}
\alias{supercell_estimate_velocity}
\title{Run RNAvelocity for super-cells (slightly modified from \link[velocyto.R]{gene.relative.velocity.estimates})
Not yet adjusted for super-cell size (not sample-weighted)}
\usage{
supercell_estimate_velocity(
  emat,
  nmat,
  smat = NULL,
  membership = NULL,
  supercell_size = NULL,
  do.run.avegaring = (ncol(emat) == length(membership)),
  kCells = 10,
  ...
)
}
\arguments{
\item{emat}{spliced (exonic) count matrix (see \link[velocyto.R]{gene.relative.velocity.estimates})}

\item{nmat}{unspliced (nascent) count matrix (\link[velocyto.R]{gene.relative.velocity.estimates})}

\item{smat}{optional spanning read matrix (used in offset calculations) (\link[velocyto.R]{gene.relative.velocity.estimates})}

\item{membership}{supercell membership ('membership' field of \link{SCimplify})}

\item{supercell_size}{a vector with supercell size (if emat and nmat provided at super-cell level)}

\item{do.run.avegaring}{whether to run averaging of emat & nmat (if nmat provided at a single-cell level)}

\item{kCells}{number of k nearest neighbors (NN) to use in slope calculation smoothing (see \link[velocyto.R]{gene.relative.velocity.estimates})}

\item{...}{other parameters from \link[velocyto.R]{gene.relative.velocity.estimates}}
}
\value{
results of \link[velocyto.R]{gene.relative.velocity.estimates} plus metacell size vector
}
\description{
Run RNAvelocity for super-cells (slightly modified from \link[velocyto.R]{gene.relative.velocity.estimates})
Not yet adjusted for super-cell size (not sample-weighted)
}
