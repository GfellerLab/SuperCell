#' Super-cells to SingleCellExperiment object
#'
#' This function transforms super-cell gene expression and super-cell partition into \link[SingleCellExperiment]{SingleCellExperiment} object
#'
#' @param SC.GE gene expression matrix with genes as rows and cells as columns
#' @param SC super-cell (output of \code{\link{SCimplify}} function)
#' @param fields which fields of \code{SC} to use as cell metadata
#' @param var.genes set of genes used as a set of variable features of SingleCellExperiment (by default is the set of genes used to generate super-cells)
#' @param is.log.normalized whether \code{SC.GE} is log-normalized counts. If yes, then SingleCellExperiment field \code{assay name = 'logcounts'} else \code{assay name = 'counts'}
#' @param do.center whether to center gene expression matrix to compute PCA
#' @param do.scale whether to scale gene expression matrix to compute PCA
#' @param ncomponents number of principal components to compute
#'
#' @details
#'
#' @return \link[SingleCellExperiment]{SingleCellExperiment} object
#'
#'@examples
#'\dontrun{
#' data(cell_lines)
#' SC           <- SCimplify(cell_lines$GE, gamma = 20)
#' SC$ident     <- supercell_assign(clusters = cell_lines$meta, supercell_membership = SC$membership)
#' SC.GE        <- supercell_GE(cell_lines$GE, SC$membership)
#' sce          <- supercell_2_sce(SC.GE = SC.GE, SC = SC, fields = c("ident"))
#'}
#' @export

supercell_2_sce <- function(SC.GE, SC, fields = c(),
                            var.genes = NULL,
                            is.log.normalized = TRUE,
                            do.center = TRUE,
                            do.scale = TRUE,
                            ncomponents = 50){

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

  ## If fields is numerical (not recommended), map them to names
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


  ## If SC.GE is log-normalized gene expression, than field data has to be rewritten
  if(is.log.normalized){
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts=SC.GE), # load as normalized data
                                                      colData = meta,
                                                      rowData = data.frame(gene_names = SC.GE@Dimnames[[1]]))
  } else {
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts=SC.GE),
                                                      colData = meta,
                                                      rowData = data.frame(gene_names = SC.GE@Dimnames[[1]]))
  }

  ## Sample-weighted scaling
  assay(sce, "scale.data") <- t(as.matrix(corpcor::wt.scale(Matrix::t(SC.GE),
                                                            w = meta$size,
                                                            center = do.center,
                                                            scale = do.scale)))


  if(is.null(var.genes)){
    var.genes <- sort(SC$genes.use)
  }

  metadata(sce) <- list(var.genes = var.genes)


  my_pca <- supercell_prcomp(X = Matrix::t(SC.GE), genes.use = var.genes,
                             fast.pca = TRUE,
                             supercell_size = meta$supercell_size,
                             k = ncomponents,
                             do.scale = do.scale, do.center = do.center,
                             double.centering = FALSE)

  colnames(my_pca$x) <- paste0('PC', 1:ncol(my_pca$x))
  rownames(my_pca$x) <- NULL
  SingleCellExperiment::reducedDim(sce, "PCA") <- my_pca$x
  SingleCellExperiment::reducedDim(sce, "PCA_weighted") <- my_pca$x

  return(sce)
}
