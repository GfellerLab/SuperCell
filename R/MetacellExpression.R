#' MetacellExpression
#'
#' \code{MetacellExpression} 
#' Compute metacells from a Seurat single-cell object. 
#' @param object A Seurat single-cell object. It has to be preprocessed (eg. latent space computed) for the assay(s) used to identify metacells
#' @return  A Seurat object with aggregated data
#' @import Seurat
#' @import SeuratObject
#' @export
#' 

MetacellExpression <- function(object, pb.method = "aggregate", assays = NULL, features = NULL,
                               return.seurat = TRUE, group.by = "ident", add.ident = NULL,
                               layer = "counts", verbose = TRUE, ...)
{
  SeuratObject::CheckDots(..., fxns = "CreateSeuratObject")
  if (!is.null(x = add.ident)) {
    .Deprecated(msg = "'add.ident' is a deprecated argument, please use the 'group.by' argument instead")
    group.by <- c("ident", add.ident)
  }
  if (!(pb.method %in% c("average", "aggregate"))) {
    stop("'pb.method' must be either 'average' or 'aggregate'")
  }
  object.assays <- .FilterObjects(object = object, classes.keep = c("Assay",
                                                                    "Assay5"))
  assays <- assays %||% object.assays
  if (!all(assays %in% object.assays)) {
    assays <- assays[assays %in% object.assays]
    if (length(x = assays) == 0) {
      stop("None of the requested assays are present in the object")
    }
    else {
      warning("Requested assays that do not exist in object. Proceeding with existing assays only.")
    }
  }
  if (length(x = layer) == 1) {
    layer <- rep_len(x = layer, length.out = length(x = assays))
  }
  else if (length(x = layer) != length(x = assays)) {
    stop("Number of layers provided does not match number of assays")
  }
  data <- FetchData(object = object, vars = rev(x = group.by))
  data <- data[which(rowSums(x = is.na(x = data)) == 0), ,
               drop = F]
  if (nrow(x = data) < ncol(x = object)) {
    message("Removing cells with NA for 1 or more grouping variables")
    object <- subset(x = object, cells = rownames(x = data))
  }
  for (i in 1:ncol(x = data)) {
    data[, i] <- as.factor(x = data[, i])
  }
  num.levels <- sapply(X = 1:ncol(x = data), FUN = function(i) {
    length(x = levels(x = data[, i]))
  })
  if (any(num.levels == 1)) {
    message(paste0("The following grouping variables have 1 value and will be ignored: ",
                   paste0(colnames(x = data)[which(num.levels <= 1)],
                          collapse = ", ")))
    group.by <- colnames(x = data)[which(num.levels > 1)]
    data <- data[, which(num.levels > 1), drop = F]
  }
  if (ncol(x = data) == 0) {
    message("All grouping variables have 1 value only. Computing across all cells.")
    category.matrix <- matrix(data = 1, nrow = ncol(x = object),
                              dimnames = list(Cells(x = object), "all"))
    if (pb.method == "average") {
      category.matrix <- category.matrix/sum(category.matrix)
    }
  }
  else {
    category.matrix <- Matrix::sparse.model.matrix(object = as.formula(object = paste0("~0+",
                                                                                       paste0("data[,", 1:length(x = group.by), "]", collapse = ":"))))
    #print(dim(category.matrix))
    colsums <- Matrix::colSums(x = category.matrix)
    category.matrix <- category.matrix[, colsums > 0]
    colsums <- colsums[colsums > 0]
    if (pb.method == "average") {
      category.matrix <- Seurat:::Sweep(x = category.matrix, MARGIN = 2,
                                        STATS = colsums, FUN = "/")
    }
    colnames(x = category.matrix) <- sapply(X = colnames(x = category.matrix),
                                            FUN = function(name) {
                                              name <- gsub(pattern = "data\\[, [1-9]*\\]",
                                                           replacement = "", x = name)
                                              return(paste0(rev(x = unlist(x = strsplit(x = name,
                                                                                        split = ":"))), collapse = "_"))
                                            })
  }
  data.return <- list()
  for (i in 1:length(x = assays)) {
    data.use <- GetAssayData(object = object, assay = assays[i],
                             layer = layer[i])
    features.to.avg <- features %||% rownames(x = data.use)
    if (inherits(x = features, what = "list")) {
      features.to.avg <- features[i]
    }
    if (IsMatrixEmpty(x = data.use)) {
      warning("The ", layer[i], " layer for the ", assays[i],
              " assay is empty. Skipping assay.", immediate. = TRUE,
              call. = FALSE)
      next
    }
    bad.features <- setdiff(x = features.to.avg, y = rownames(x = data.use))
    if (length(x = bad.features) > 0) {
      warning("The following ", length(x = bad.features),
              " features were not found in the ", assays[i],
              " assay: ", paste(bad.features, collapse = ", "),
              call. = FALSE, immediate. = TRUE)
    }
    features.assay <- intersect(x = features.to.avg, y = rownames(x = data.use))
    if (length(x = features.assay) > 0) {
      data.use <- data.use[features.assay, ]
    }
    else {
      warning("None of the features specified were found in the ",
              assays[i], " assay.", call. = FALSE, immediate. = TRUE)
      next
    }
    data.return[[i]] <- as.sparse(x = (data.use %*% category.matrix))
    colnames(data.return[[i]]) <- paste0("Metacell_", c(1:ncol(data.return[[i]])))
    names(x = data.return)[i] <- assays[[i]]
  }
  if (return.seurat) {
    toRet <- CreateSeuratObject(counts = data.return[[1]],
                                project = if (pb.method == "average")
                                  "Average"
                                else "Aggregate", assay = names(x = data.return)[1],
                                ...)
    if (length(x = data.return) > 1) {
      for (i in 2:length(x = data.return)) {
        toRet[[names(x = data.return)[i]]] <- CreateAssay5Object(counts = data.return[[i]])
      }
    }
    if (DefaultAssay(object = object) %in% names(x = data.return)) {
      DefaultAssay(object = toRet) <- DefaultAssay(object = object)
    }
    if ("ident" %in% group.by) {
      first.cells <- c()
      for (i in 1:ncol(x = category.matrix)) {
        first.cells <- c(first.cells, Position(x = category.matrix[,
                                                                   i], f = function(x) {
                                                                     x > 0
                                                                   }))
      }
      Idents(object = toRet) <- Idents(object = object)[first.cells]
    }
    return(toRet)
  }
  else {
    return(data.return)
  }
}
