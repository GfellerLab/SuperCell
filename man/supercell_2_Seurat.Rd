% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/supercell_2_Seurat.R
\name{supercell_2_Seurat}
\alias{supercell_2_Seurat}
\title{Super-cells to Seurat object}
\usage{
supercell_2_Seurat(
  SC.GE,
  SC,
  fields = c(),
  var.genes = NULL,
  do.preproc = TRUE,
  is.log.normalized = TRUE,
  do.center = TRUE,
  do.scale = TRUE,
  N.comp = NULL,
  output.assay.version = "v4"
)
}
\arguments{
\item{SC.GE}{gene expression matrix with genes as rows and cells as columns}

\item{SC}{super-cell (output of \link{SCimplify} function)}

\item{fields}{which fields of \code{SC} to use as cell metadata}

\item{var.genes}{set of genes used as a set of variable features of Seurat (by default is the set of genes used to generate super-cells), ignored if \code{!do.preproc}}

\item{do.preproc}{whether to do prepocessing, including data normalization, scaling, HVG, PCA, nearest neighbors, \code{TRUE} by default, change to \code{FALSE} to speed up conversion}

\item{is.log.normalized}{whether \code{SC.GE} is log-normalized counts. If yes, then Seurat field \code{data} is replaced with \code{counts} after normalization (see 'Details' section), ignored if \code{!do.preproc}}

\item{do.center}{whether to center gene expression matrix to compute PCA, ignored if \code{!do.preproc}}

\item{do.scale}{whether to scale gene expression matrix to compute PCA, ignored if \code{!do.preproc}}

\item{N.comp}{number of principal components to use for construction of single-cell kNN network, ignored if \code{!do.preproc}}

\item{output.assay.version}{version of the seurat assay in output, \code{`"v4"`} by default, \code{`"v5"`} requires \link[Seurat]{Seurat} v5 installed.}
}
\value{
\link[Seurat]{Seurat} object
}
\description{
This function transforms super-cell gene expression and super-cell partition into \link[Seurat]{Seurat} object
}
\details{
Since the input of \link[Seurat]{CreateSeuratObject} should be unnormalized count matrix (UMIs or TPMs, see \link[Seurat]{CreateSeuratObject}).
Thus, we manually set field \code{`assays$RNA@data`} to \code{SC.GE} if \code{is.log.normalized == TRUE}.
Avoid running \link[Seurat]{NormalizeData} for the obtained Seurat object, otherwise this will overwrite field \code{`assays$RNA@data`}.
If you have run \link[Seurat]{NormalizeData}, then make sure to replace \code{`assays$RNA@data`} with correct matrix by running
\code{`your_seurat@assays$RNA@data <- your_seurat@assays$RNA@counts`}.

Since super-cells have different size (consist of different number of single cells), we use sample-weighted algorithms for all
 possible steps of the downstream analysis, including scaling and dimensionality reduction. Thus, generated Seurat object  comes
 with the results of sample-wighted scaling (available as \code{`your_seurat@assays$RNA@scale.data`} or
  \code{`your_seurat@assays$RNA@misc[["scale.data.weighted"]]`} to reproduce if the first one has been overwritten) and PCA (available as
\code{`your_seurat@reductions$pca`} or \code{`your_seurat@reductions$pca_weighted`} to reproduce if the first one has been overwritten).
}
\examples{
\dontrun{
data(cell_lines)
SC           <- SCimplify(cell_lines$GE, gamma = 20)
SC$ident     <- supercell_assign(clusters = cell_lines$meta, supercell_membership = SC$membership)
SC.GE        <- supercell_GE(cell_lines$GE, SC$membership)
m.seurat     <- supercell_2_Seurat(SC.GE = SC.GE, SC = SC, fields = c("ident"))
}
}
