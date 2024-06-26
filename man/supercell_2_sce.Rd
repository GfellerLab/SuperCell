% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/supercell_2_sce.R
\name{supercell_2_sce}
\alias{supercell_2_sce}
\title{Super-cells to SingleCellExperiment object}
\usage{
supercell_2_sce(
  SC.GE,
  SC,
  fields = c(),
  var.genes = NULL,
  do.preproc = TRUE,
  is.log.normalized = TRUE,
  do.center = TRUE,
  do.scale = TRUE,
  ncomponents = 50
)
}
\arguments{
\item{SC.GE}{gene expression matrix with genes as rows and cells as columns}

\item{SC}{super-cell (output of \link{SCimplify} function)}

\item{fields}{which fields of \code{SC} to use as cell metadata}

\item{var.genes}{set of genes used as a set of variable features of SingleCellExperiment (by default is the set of genes used to generate super-cells)}

\item{do.preproc}{whether to do prepocessing, including data normalization, scaling, HVG, PCA, nearest neighbors, \code{TRUE} by default, change to \code{FALSE} to speed up conversion}

\item{is.log.normalized}{whether \code{SC.GE} is log-normalized counts. If yes, then SingleCellExperiment field \code{assay name = 'logcounts'} else \code{assay name = 'counts'}}

\item{do.center}{whether to center gene expression matrix to compute PCA}

\item{do.scale}{whether to scale gene expression matrix to compute PCA}

\item{ncomponents}{number of principal components to compute}
}
\value{
\link[SingleCellExperiment]{SingleCellExperiment} object
}
\description{
This function transforms super-cell gene expression and super-cell partition into \link[SingleCellExperiment]{SingleCellExperiment} object
}
\examples{
\dontrun{
data(cell_lines)
SC           <- SCimplify(cell_lines$GE, gamma = 20)
SC$ident     <- supercell_assign(clusters = cell_lines$meta, supercell_membership = SC$membership)
SC.GE        <- supercell_GE(cell_lines$GE, SC$membership)
sce          <- supercell_2_sce(SC.GE = SC.GE, SC = SC, fields = c("ident"))
}
}
