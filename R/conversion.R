#' Convert Metacells (Metacell-2) to Super-cell like object
#' @param adata anndata object of metacells (the output of \href{https://metacells.readthedocs.io/en/latest/_modules/metacells/pipeline/collect.html#collect_metacells}{\code{collect_metacells()}})
#' @param obs.sc a dataframe of the single-cell anndata object used to compute metacells (anndata after applying \href{https://metacells.readthedocs.io/en/latest/_modules/metacells/pipeline/divide_and_conquer.html#divide_and_conquer_pipeline}{\code{divide_and_conquer_pipeline()}} function)
#' @return a list of super-cell like object (similar to the output of \link{SCimplify})
#'
#' @export

metacell2_anndata_2_supercell <- function(
  adata,
  obs.sc
){
  res <- list()

  N.c                <- nrow(obs.sc)
  cell.ids           <- obs.sc$X
  cells.use.ids      <- cell.ids[obs.sc$metacell > -1] # '-1' corresponds to outliers and '-2' for non-clear cells

  membership         <- obs.sc$metacell + 1
  names(membership)  <- cell.ids
  membership[membership<=0] <- NA

  cells.use.idx      <- pmatch(cells.use.ids, cell.ids)
  supercell_size     <- as.vector(table(membership))

  gamma              <- adata$uns$gamma

  N.SC               <- length(supercell_size)
  gamma.actual       <- round(N.c/N.SC)

  genes.use          <- rownames(adata$var)[adata$var$feature_gene == 1]


  obs                <- adata$obs
  rownames(obs)      <- as.character(1:nrow(obs))

  counts <- Matrix::t(adata$chunk_X())
  rownames(counts) <- rownames(adata$var)
  if(is.null(colnames(counts))){
    colnames(counts) <- as.character(1:ncol(counts))
  }


  res <- list(
    membership            = membership,
    supercell_size        = supercell_size,
    genes.use             = genes.use,
    simplification.algo   = "Metacell2",
    gamma                 = gamma,
    gamma.actual          = gamma.actual,
    cells.use.ids         = cells.use.ids,
    cells.use.idx         = cells.use.idx,
    SC.counts             = counts,
    SC.meta               = obs
  )

  return(res)
}


#' Convert Anndata metacell object (Metacell-2 or SEACells) to Super-cell like object
#' @param adata anndata object of metacells
#' (for example, the output of \href{https://metacells.readthedocs.io/en/latest/_modules/metacells/pipeline/collect.html#collect_metacells}{\code{collect_metacells()}} for Metacells
#' or the output of \href{https://github.com/dpeerlab/SEACells/blob/406ebd3a7ee200788358464b345fa0c8d88fda69/SEACells/core.py#L532}{\code{SEACells.core.summarize_by_SEACell}})
#'
#' Please, **make sure**, adata has `uns['sc.obs']` field containing observation information of single-cell data, in particular,
#' a column 'membership' (single-cell assignemnt to metacells)
#'
#' @param simplification.algo metacell construction algorithm (i.e., Metacell2 or SEACells)
#'
#' @return a list of super-cell like object (similar to the output of \link{SCimplify})
#'
#'
#' @export

anndata_2_supercell <- function(
  adata,
  simplification.algo = "unknown"
){
  res <- list()
  if(!("uns" %in% names(adata))){
    stop("adata has to contain 'uns' field!")
  }
  if(!("sc.obs" %in% names(adata[["uns"]]))){
    stop("adata[['uns']] has to contain 'sc.obs' field!")
  }

  obs.sc             <- adata[["uns"]][["sc.obs"]]

  if(!("membership" %in% colnames(obs.sc))){
    stop("adata[['uns']][[sc.obs]] has to contain 'membership' column!")
  }

  N.c                <- nrow(obs.sc)
  cell.ids           <- obs.sc$X
  cells.use.ids      <- cell.ids[!is.na(obs.sc$membership)]

  membership         <- obs.sc$membership
  names(membership)  <- cell.ids

  cells.use.idx      <- pmatch(cells.use.ids, cell.ids)
  supercell_size     <- as.vector(table(membership))

  gamma              <- adata$uns$gamma

  N.SC               <- length(supercell_size)
  gamma.actual       <- round(N.c/N.SC)

  genes.use          <- rownames(adata$var)[adata$var$feature_gene == 1]


  obs                <- adata$obs
  rownames(obs)      <- as.character(1:nrow(obs))

  counts <- Matrix::t(adata$chunk_X())
  rownames(counts) <- rownames(adata$var)
  if(is.null(colnames(counts))){
    colnames(counts) <- as.character(1:ncol(counts))
  }


  res <- list(
    membership            = membership,
    supercell_size        = supercell_size,
    genes.use             = genes.use,
    simplification.algo   = simplification.algo,
    gamma                 = gamma,
    gamma.actual          = gamma.actual,
    cells.use.ids         = cells.use.ids,
    cells.use.idx         = cells.use.idx,
    SC.counts             = counts,
    SC.meta               = obs
  )

  return(res)
}
