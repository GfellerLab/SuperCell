#' Super-cells to Seurat object
#'
#' This function transforms super-cell gene expression and super-cell partition into \link[Seurat]{Seurat} object
#'
#'
#' @param SC.GE gene expression matrix with genes as rows and cells as columns
#' @param SC super-cell (output of \code{\link{SCimplify}} function)
#' @param fields which fields of \code{SC} to use as cell metadata
#' @param var.genes set of genes used as a set of variable features of Seurat (by default is the set of genes used to generate super-cells)
#' @param is.log.normalized whether \code{SC.GE} is log-normalized counts. If yes, then Seurat field \code{data} is replaced with \code{counts} after normalization (see 'Details' section)
#' @param do.center whether to center gene expression matrix to compute PCA
#' @param do.scale whether to scale gene expression matrix to compute PCA
#'
#'
#' @details
#' Since the input of \link[Seurat]{CreateSeuratObject} should be unnormalized count matrix (UMIs or TPMs, see \link[Seurat]{CreateSeuratObject}).
#' Thus, we manually set field \code{`assays$RNA@data`} to \code{SC.GE} if \code{is.log.normalized == TRUE}.
#' Avoid running \link[Seurat]{NormilizeData} for the obtained Seurat object, otherwise this will overwtite fieid \code{`assays$RNA@data`}.
#' If you have run \link[Seurat]{NormilizeData}, then make sure to replace \code{`assays$RNA@data`} with correct matrix by running
#' \code{`your_seurat@assays$RNA@data <- your_seurat@assays$RNA@counts`}.
#'
#' Since super-cells have different size (consist of different number of single cells), we use sample-weighted algorithms for all
#'  possible steps of the downstream analysis, including scaling and dimensionality reduction. Thus, generated Seurat object  comes
#'  with the results of sample-wighted scaling (available as \code{`your_seurat@assays$RNA@scale.data`} or
#'   \code{`your_seurat@assays$RNA@misc[["scale.data.weighted"]]`} to reproduce if the first one has been overwritten) and PCA (available as
#' \code{`your_seurat@reductions$pca`} or \code{`your_seurat@reductions$pca_weighted`} to reproduce if the first one has been overwritten).
#'
#'
#'
#' @return \link[Seurat]{Seurat} object
#'
#'@examples
#'\dontrun{
#' data(cell_lines)
#' SC           <- SCimplify(cell_lines$GE, gamma = 20)
#' SC$ident     <- supercell_assign(clusters = cell_lines$meta, supercell_membership = SC$membership)
#' SC.GE        <- supercell_GE(cell_lines$GE, SC$membership)
#' m.seurat     <- supercell_2_Seurat(SC.GE = SC.GE, SC = SC, fields = c("ident"))
#'}
#' @export


supercell_2_Seurat <- function(SC.GE, SC, fields = c(),
                               var.genes = NULL,
                               is.log.normalized = TRUE,
                               do.center = TRUE,
                               do.scale = TRUE){
  N.c <- ncol(SC.GE)
  if(is.null(SC$supercell_size)){
    warning(paste0("supercell_size field of SC is missing, size of all super-cells set to 1"))
    supercell_size <- rep(1, N.c)
  } else {
    supercell_size <- SC$supercell_size
  }

  if(length(supercell_size) != N.c){
    stop(paste0("length of SC$supercell_size has to be the same as number of super-cells ", N.c))
  }

  ## Name all cells to create Seurat Object
  if(is.null(colnames(SC.GE))){
    colnames(SC.GE) <- as.character(1:N.c)
  }
  counts <- SC.GE

  ## If fields is numerical, map them to names
  if(is.numeric(fields)){
    fields <- names(SC)[fields]
  }

  ## Keep only available fiedls
  fields <- intersect(fields, names(SC))

  if(length(fields) > 0){
    SC.fields <- SC[fields]
  } else {
    SC.fields <- NULL
  }

  ## Keep only fields that are specific to cells
  SC.field.length <- lapply(SC.fields, length)
  SC.fields       <- SC.fields[which(SC.field.length == N.c)]

  meta     <- data.frame(size = SC$supercell_size, row.names = colnames(SC.GE), stringsAsFactors = FALSE)

  if(length(SC.fields) > 0){
    meta <- cbind(meta, SC.fields)
  }
  m.seurat <- Seurat::CreateSeuratObject(counts = SC.GE, meta.data = meta)

  ## Normalize data, so Seurat does not generate warning
  m.seurat <- Seurat::NormalizeData(m.seurat)

  ## If SC.GE is log-normalized gene expression, than field data has to be rewritten
  if(is.log.normalized){
    m.seurat@assays$RNA@data <- m.seurat@assays$RNA@counts
  }

  ## Sample-weighted scaling
  m.seurat@assays$RNA@scale.data <- t(as.matrix(corpcor::wt.scale(Matrix::t((m.seurat@assays$RNA@data)),
                                                                w = meta$size,
                                                                center = do.center,
                                                                scale = do.scale)))


  m.seurat@assays$RNA@misc[["scale.data.weighted"]] <- m.seurat@assays$RNA@scale.data

  if(is.null(var.genes)){
    var.genes <- sort(SC$genes.use)
  }

  VariableFeatures(m.seurat) <- var.genes
  m.seurat <- RunPCA(m.seurat, verbose = F)
  m.seurat@reductions$pca_seurat <- m.seurat@reductions$pca

  my_pca <- supercell_prcomp(X = Matrix::t(SC.GE), genes.use = var.genes,
                             fast.pca = TRUE,
                             supercell_size = meta$supercell_size,
                             k = dim(m.seurat@reductions$pca@cell.embeddings)[2],
                             do.scale = do.scale, do.center = do.center,
                             double.centering = FALSE)

  dimnames(my_pca$x) <- dimnames(m.seurat@reductions$pca_seurat)
  m.seurat@reductions$pca@cell.embeddings  <- my_pca$x
  m.seurat@reductions$pca@feature.loadings <- my_pca$rotation
  m.seurat@reductions$pca@stdev            <- my_pca$sdev

  m.seurat@reductions$pca_weighted         <- m.seurat@reductions$pca

  ## Super-cell netwok:
  ## 1) create graph field
  m.seurat            <- FindNeighbors(m.seurat, compute.SNN = TRUE, verbose = FALSE)

  ## 2) add self-loops to our super-cell graph to indicate super-cell size (does not work, as Seurat removes loops...)
 # SC$graph.supercells <- igraph::add_edges(SC$graph.supercells, edges = rep(1:N.c, each = 2), weight = supercell_size)
  adj.mtx             <- igraph::get.adjacency(SC$graph.supercells, attr = "weight")

  ## 3) replace generated Seurat network with the super-cell network
  m.seurat@graphs$RNA_nn@i                <- adj.mtx@i
  m.seurat@graphs$RNA_nn@p                <- adj.mtx@p
  m.seurat@graphs$RNA_nn@Dim              <- adj.mtx@Dim
  m.seurat@graphs$RNA_nn@x                <- adj.mtx@x
  m.seurat@graphs$RNA_nn@factors          <- adj.mtx@factors

  m.seurat@graphs$RNA_super_cells         <- m.seurat@graphs$RNA_nn
  return(m.seurat)
}

