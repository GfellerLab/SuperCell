#' Differential expression analysis of supep-cell data. Most of the parameters are the same as in Seurat \link[Seurat]{FindMarkers} (for simplicity)
#'
#'
#' @param ge gene expression matrix for super-cells (rows - genes, cols - super-cells)
#' @param supercell_size a vector with supercell size (ordered the same way as in \code{ge})
#' @param clusters a vector with clustering information (ordered the same way as in \code{ge})
#' @param ident.1 name(s) of cluster for which markers are computed
#' @param ident.2 name(s) of clusters for comparison. If \code{NULL} (defauld), then all the other clusters used
#' @param genes.use set of genes to test. Defeult -- all genes in \code{ge}
#' @param logfc.threshold log fold change threshold for genes to be considered in the further analysis
#' @param min.expr minimal expression (default 0)
#' @param min.pct remove genes with lower percentage of detection from the set of genes which will be tested
#' @param seed random seed to use
#' @param only.pos whether to compute only positive (upregulated) markers
#' @param return.extra.info whether to return extra information about test and its statistics. Default is FALSE.
#'
#' @return a matrix with a test name (t-test), statisctics, adjusted p-values, logFC, percenrage of detection in eacg ident and mean expresiion
#'
#' @export
#'



supercell_FindMarkers <- function(ge, supercell_size = NULL, clusters, ident.1, ident.2 = NULL, genes.use = NULL,
                                  logfc.threshold = 0.25, min.expr = 0., min.pct = 0.1, seed = 12345,
                                  only.pos = FALSE, return.extra.info = FALSE){

  total.number.of.genes <- nrow(ge)

  if(sum(is.na( pmatch(ident.1, unique(clusters)))) > 0){
    stop(paste("Unknown ident.1", ident.1))
  }
  if(!is.null(ident.2)){
    if(sum(is.na(pmatch(ident.2, unique(clusters)))) > 0){
      stop(paste("Unknown ident.2", ident.2))
    }
  } else { ## ident.2 == NULL
    ident.2 <- setdiff(unique(clusters), ident.1)
  }

  if(is.null(supercell_size[1])){
    supercell_size <- rep(1, ncol(ge))
  }

  cell.idx.ident.1 <- which(clusters %in% ident.1)
  cell.weigth.1    <- supercell_size[cell.idx.ident.1]

  cell.idx.ident.2 <- which(clusters %in% ident.2)
  cell.weigth.2    <- supercell_size[cell.idx.ident.2]

  if(min(length(cell.idx.ident.1), length(cell.idx.ident.2)) < 2){
    return(res = list(adj.p.value = NA, logFC = NA))
  }

  if(is.null(genes.use)){
    genes.use <- rownames(ge)
  } else {
    genes.use <- intersect(genes.use, rownames(ge))
  }

  if(length(genes.use) < 2){
    return(res = list(adj.p.value = NA, logFC = NA))
  }

  # compute percentage of cells expressing a gene and filter out poorly expressed genes
  #print("Compute percentage of genes expressed in super cells (detection rate)")

  pct.1 <- apply(ge[genes.use, cell.idx.ident.1], 1, function(x){sum(cell.weigth.1[x>min.expr])/sum(cell.weigth.1)})
  pct.2 <- apply(ge[genes.use, cell.idx.ident.2], 1, function(x){sum(cell.weigth.2[x>min.expr])/sum(cell.weigth.2)})
  max.pct.1.2 <- apply(cbind(pct.1, pct.2), 1, max)

  #print(min.pct.1.2)
  genes.min.pct <- genes.use[which(max.pct.1.2 > min.pct)] # filter out genes with low detection rate
  genes.use <- intersect(genes.use, genes.min.pct)

  if(length(genes.use) < 2){
    return(res = list(adj.p.value = NA, logFC = NA))
  }

  #print("compute average LogFC")
  w.mean.1 <- apply(ge[genes.use, cell.idx.ident.1 ], 1, function(x){Hmisc::wtd.mean(x, weights = cell.weigth.1)})
  w.mean.2 <- apply(ge[genes.use, cell.idx.ident.2 ], 1, function(x){Hmisc::wtd.mean(x, weights = cell.weigth.2)})
  ## exp of normalized expression is used to compute av.logFC
  w.mean.exp.1 <- apply(ge[genes.use, cell.idx.ident.1 ], 1, function(x){Hmisc::wtd.mean(x = expm1(x = x), weights = cell.weigth.1)}) # take true expression value (as seurat does), not log-transformed
  w.mean.exp.2 <- apply(ge[genes.use, cell.idx.ident.2 ], 1, function(x){Hmisc::wtd.mean(x = expm1(x = x), weights = cell.weigth.2)})

  logFC    <- log(w.mean.exp.1 + 1.01) - log(w.mean.exp.2 + 1.01)

  if(!only.pos){
    genes.logFC <- genes.use[which(abs(logFC) > logfc.threshold)]
  } else {
    genes.logFC <- genes.use[which(logFC > logfc.threshold)]
  }
  genes.use <- intersect(genes.use, genes.logFC)

  if(length(genes.use) < 2){
    return(res = list(adj.p.value = NA, logFC = NA))
  }


  #print("Filtered genes")



  w.t.test <- apply(ge[genes.use,], 1, function(x){weights::wtd.t.test(x = x[cell.idx.ident.1], y = x[cell.idx.ident.2],
                                                                       weight = cell.weigth.1, weighty  = cell.weigth.2,
                                                                       mean1 = FALSE, samedata = FALSE)}) # mean1 == T results in less degreas of freedom (lower p value), thus, df = number of super cells, not total number of cells


  if(length((w.t.test))>1){
    dd            <- data.frame(matrix(unlist(w.t.test), ncol = 8, byrow = T))
    colnames(dd)  <- c("test", "t.value", "df", "p.value", "difference", "mean_x", "mean_y", "std_err")
    rownames(dd)  <- names(w.t.test)

    adj.p.value <- p.adjust(as.vector(dd$p.value), method = p.adjust.methods[4], n = total.number.of.genes)

    pct.1 <- pct.1[rownames(dd)]
    pct.2 <- pct.2[rownames(dd)]
    logFC <- logFC[rownames(dd)]
    w.mean.1 <- w.mean.1[rownames(dd)]
    w.mean.2 <- w.mean.2[rownames(dd)]

    res <- cbind(dd, adj.p.value, pct.1, pct.2, logFC, w.mean.1, w.mean.2)

    for(j in colnames(res)){
      if(j != "test")
        res[,j] <- as.numeric(as.vector(res[,j]))
    }

    order <- order(adj.p.value, 1/(abs(logFC)+1), decreasing = FALSE)
    res   <- res[order,]

    if(!return.extra.info) res <- res[,-c(1:3, 5:8)]

  } else { res = list(adj.p.value = NA, logFC = NA)}



  return(res)

}




#' Differential expression analysis of supep-cell data. Most of the parameters are the same as in Seurat \link[Seurat]{FindAllMarkers} (for simplicity)
#'
#'
#' @param ge gene expression matrix for super-cells (rows - genes, cols - super-cells)
#' @param supercell_size a vector with supercell size (ordered the same way as in \code{ge})
#' @param clusters a vector with clustering information (ordered the same way as in \code{ge})
#' @param genes.use set of genes to test. Defeult -- all genes in \code{ge}
#' @param logfc.threshold log fold change threshold for genes to be considered in the further analysis
#' @param min.expr minimal expression (default 0)
#' @param min.pct remove genes with lower percentage of detection from the set of genes which will be tested
#' @param seed random seed to use
#' @param only.pos whether to compute only positive (upregulated) markers
#' @param return.extra.info whether to return extra information about test and its statistics. Default is FALSE.
#'
#' @return list of results of \link[SCimple]{supercell_FindMarkers}
#'
#' @export
#'

supercell_FindAllMarkers <- function(ge, supercell_size = NULL, clusters, genes.use = NULL, logfc.threshold = 0.25,
                                     min.expr = 0., min.pct = 0.1, seed = 12345, only.pos = FALSE, return.extra.info = FALSE){
  res <- list()
  unique.ident <- sort(unique(clusters))
  for(ident.1 in unique.ident){
    res[[ident.1]] <- supercell_FindMarkers(ge = ge,
                                            supercell_size = supercell_size,
                                            clusters = clusters,
                                            ident.1 = ident.1,
                                            ident.2 = NULL,
                                            genes.use = genes.use,
                                            logfc.threshold = logfc.threshold,
                                            min.expr = min.expr,
                                            min.pct = min.pct,
                                            seed = seed,
                                            only.pos = only.pos,
                                            return.extra.info = return.extra.info)

  }
  return(res)
}
