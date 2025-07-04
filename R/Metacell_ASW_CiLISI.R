#' MetacellBioSingleASW
#'
#' \code{MetacellBioSingleASW} 
#' Compute silhouette of bio labels of the metacells at the single cell levels. (relies on scIntegrationMetrics and Seurat R packages). 
#' @param sobj.mc A metacell suerat object
#' @param meta.data.sc A meta data dataframe of single cells..
#' @param reduction.name Metacell latent space to use to compute the silhouettes. 
#' @param label.col label column in metacell meta data to use to compute the silhouettes.
#' @param return.all.res whether to return all silhouettes coefficients or only the average (default).
#' @param return.all.res dimensions to use to compute the silhouettes (default all computed present in the Seurat object).
#' @param kept.labels labels to keep to compute the silhouettes (default all). Useful to discard annotated cells.
#' @return  A vector containing the compactness of each metacell.
#' @examples
#' bioASW <- MetacellBioSingleASW(sobj.mc = sobj.mc,
#'                          reduction.name = "integrated.rpca",
#'                          meta.data.sc = meta.data.sc,
#'                          dims = 1:50,
#'                          label.col = "label")
#' @import Seurat
#' @export

MetacellBioSingleASW <- function(sobj.mc,
                                 meta.data.sc,
                                 reduction.name,
                                 label.col,
                                 return.all.res = FALSE,
                                 dims = NULL,
                                 kept.labels = NULL) {
  
  if (!is.null(kept.labels)) {
    sobj.mc <- sobj.mc[,sobj.mc[[label.col]][,1] %in% kept.labels]
    sobj.mc@misc$membership <- sobj.mc@misc$membership[sobj.mc@misc$membership %in% Cells(sobj.mc)]
    meta.data.sc <- meta.data.sc[names(sobj.mc@misc$membership),]
  }
  if (is.null(dims)) {
    dims <- 1:dim(sobj.mc[[reduction.name]])[2]
  }
  metacells.embeddings <- Seurat::FetchData(sobj.mc,vars = names(sobj.mc[[reduction.name]])[dims]) 
  expanded.embeddings <- metacells.embeddings[sobj.mc@misc$membership[rownames(meta.data.sc)],]
  rownames(expanded.embeddings) <- rownames(meta.data.sc)
  
  metacells.meta.data <- FetchData(sobj.mc,vars = c("orig.ident",label.col))
  expanded.meta.data <- metacells.meta.data[sobj.mc@misc$membership[rownames(meta.data.sc)],]
  rownames(expanded.meta.data) <- rownames(meta.data.sc)
  res <- scIntegrationMetrics::compute_silhouette(as.matrix(expanded.embeddings),
                                                  expanded.meta.data,
                                                  label_colnames = c(label.col))
  if(!return.all.res) {
    res <- mean(res[,label.col])
  }
  
  return(res)
} 


#' MetacellBioSplitASW
#'
#' \code{MetacellBioSplitASW} 
#' Compute weigted silhouette of bio labels of the metacells after splitting. (relie on Seurat R packages). 
#' @param sobj.mc A metacell suerat object
#' @param meta.data.sc A meta data dataframe of single cells..
#' @param reduction.name Metacell latent space to use to compute the silhouettes. 
#' @param label.col label column in metacell meta data to use to compute the silhouettes.
#' @param return.all.res whether to return all silhouettes coefficients or only the average (default).
#' @param return.all.res dimensions to use to compute the silhouettes (default all computed present in the Seurat object).
#' @param kept.labels labels to keep to compute the silhouettes (default all). Useful to discard annotated cells.
#' @return  A vector containing the compactness of each metacell.
#' @examples
#' bioASW <- MetacellBioSplitASW(sobj.mc = sobj.mc,
#'                          reduction.name = "integrated.rpca",
#'                          meta.data.sc = meta.data.sc,
#'                          dims = 1:50,
#'                          label.col = "label")
#' @import Seurat
#' @export

MetacellBioSplitASW <- function(sobj.mc,
                                meta.data.sc,
                                reduction.name,
                                label.col,
                                return.all.res = FALSE,
                                dims = NULL,
                                kept.labels = NULL) {
  
  if (!is.null(kept.labels)) {
    # discard single cells of non kept labels
    sobj.mc@misc$membership <- sobj.mc@misc$membership[meta.data.sc[names(sobj.mc@misc$membership),label.col] %in% kept.labels]
    # discard metacells of non kept labels
    sobj.mc <- sobj.mc[,sobj.mc[[label.col]][,1] %in% kept.labels]
    sobj.mc@misc$membership <- sobj.mc@misc$membership[sobj.mc@misc$membership %in% Cells(sobj.mc)]
    meta.data.sc <- meta.data.sc[names(sobj.mc@misc$membership),]
    
  }
  if (is.null(dims)) {
    dims <- 1:dim(sobj.mc[[reduction.name]])[2]
  }
  splitted.size <- table(interaction(sobj.mc@misc$membership,meta.data.sc[names(sobj.mc@misc$membership),label.col], drop = TRUE,sep = ";;"))
  
  splitted.membership <- names(splitted.size)
  
  x <- stringr::str_split_fixed(splitted.membership,";;",n=2)[,2]
  metacells <- stringr::str_split_fixed(splitted.membership,";;",n=2)[,1]
  splitted.embeddings <- Embeddings(sobj.mc[[reduction.name]])[metacells,dims]
  
  rownames(splitted.embeddings) <- splitted.membership
  
  dist.metacells <- dist(splitted.embeddings)
  res <- supercell_silhouette(x, dist.metacells, supercell_size = splitted.size)
  #print(res)
  if(!return.all.res) {
    res <- res$avg.width
  }
  
  return(res)
  
}

#' MetacellBatchSplitASW
#'
#' \code{MetacellBatchSplitASW} 
#' Compute weigted silhouette of batch labels of the metacells after splitting. (relie on Seurat R packages). 
#' Compute weighted average silhouette for batch label in each cell type and take the mean (based on batch_ASW proposed by Luecken et al 2022) 
#' @param sobj.mc A metacell seurat object
#' @param meta.data.sc A meta data dataframe of single cells..
#' @param reduction.name Metacell latent space to use to compute the silhouettes. 
#' @param label.col column name of bio labels in metacell meta data .
#' @param batch.col column name of batch labels in metacell meta data to use to compute the silhouettes.
#' @param return.all.res whether to return all silhouettes coefficients or only the average (default).
#' @param return.all.res dimensions to use to compute the silhouettes (default all computed present in the Seurat object).
#' @param kept.labels labels to keep to compute the silhouettes (default all). Useful to discard annotated cells.
#' @return  A vector containing the compactness of each metacell.
#' @examples
#' batch <- MetacellBatchSplitASW(sobj.mc = sobj.mc,
#'                          reduction.name = "integrated.rpca",
#'                          meta.data.sc = meta.data.sc,
#'                          dims = 1:50,
#'                          batch.label.col = "batch",
#'                          label.col = "label")
#' @import Seurat
#' @export


MetacellBatchSplitASW <- function(sobj.mc,
                                  reduction.name,
                                  meta.data.sc,
                                  bio.label.col,
                                  batch.label.col,
                                  kept.labels = NULL,
                                  dims = NULL,
                                  return.res.all = F,
                                  min.obs = 10) {
  
  if (!is.null(kept.labels)) {
    sobj.mc <- sobj.mc[,sobj.mc[[bio.label.col]][,1] %in% kept.labels]
    sobj.mc@misc$membership <- sobj.mc@misc$membership[sobj.mc@misc$membership %in% Cells(sobj.mc)]
    meta.data.sc <- meta.data.sc[names(sobj.mc@misc$membership),]
    
    meta.data.sc <- meta.data.sc[meta.data.sc[[bio.label.col]] %in% kept.labels,]
    sobj.mc@misc$membership <- sobj.mc@misc$membership[rownames(meta.data.sc)]
    sobj.mc <- sobj.mc[,Cells(sobj.mc) %in% sobj.mc@misc$membership]
  } else {
    kept.labels <- unique(sobj.mc[[bio.label.col]][,1])
  }
  if (is.null(dims)) {
    dims <- 1:dim(sobj.mc[[reduction.name]])[2]
  }
  
  kept.labels <- names(which(table(sobj.mc[[bio.label.col]][,1]) > min.obs & names(table(sobj.mc[[bio.label.col]][,1])) %in% kept.labels))
  print(kept.labels)
  splitted.size <- table(interaction(sobj.mc@misc$membership,
                                     meta.data.sc[names(sobj.mc@misc$membership),bio.label.col],
                                     meta.data.sc[names(sobj.mc@misc$membership),batch.label.col], 
                                     drop = TRUE,sep = ";;"))
  
  splitted.membership <- names(splitted.size)
  
  bio <- stringr::str_split_fixed(splitted.membership,";;",n=3)[,2]
  batches <-  stringr::str_split_fixed(splitted.membership,";;",n=3)[,3]
  
  metacells <- stringr::str_split_fixed(splitted.membership,";;",n=3)[,1]
  
  splitted.embeddings <- Embeddings(sobj.mc[[reduction.name]])[metacells,dims]
  rownames(splitted.embeddings) <- splitted.membership
  
  
  
  
  asw.batch.ct <- lapply(kept.labels, function(X) {
    sub.metacells <- splitted.membership[bio == X]
    dist.metacells <- dist(splitted.embeddings[sub.metacells,])
    sub.splitted.size <- splitted.size[bio == X]
    weighted.sil.res <- supercell_silhouette(x = batches[bio == X],
                                             dist.metacells, 
                                             supercell_size = sub.splitted.size)
    abs.s <- abs(as.numeric(weighted.sil.res$s[,"silhouette width"]))
    abs.s.batch <-  1 - abs.s
    res <- sum(abs.s.batch*sub.splitted.size)/sum(sub.splitted.size)
    res.all <- list("abs.s.batch" = abs.s.batch,"abs.s.batch.avg" = res)
    return(res.all)
  }
  )
  names(asw.batch.ct) <- kept.labels
  if (return.res.all) {
    return(asw.batch.ct)
  } else {
    asw.batch.ct <- lapply(asw.batch.ct,function(x){return(x$abs.s.batch)})
    asw.batch.ct <- mean(do.call(c,asw.batch.ct))
  }
  return(asw.batch.ct)
}



MetacellSplitCiLISI <- function(sobj.mc,
                                reduction.name,
                                meta.data.sc,
                                bio.label.col,
                                batch.label.col,
                                kept.labels = NULL,
                                dims = NULL,
                                return.res.all = F,
                                min.obs = 10) {
  if (!is.null(kept.labels)) {
    # discard single cells of non kept labels
    sobj.mc@misc$membership <- sobj.mc@misc$membership[meta.data.sc[names(sobj.mc@misc$membership),bio.label.col] %in% kept.labels]
    # discard metacells of non kept labels
    sobj.mc <- sobj.mc[,sobj.mc[[bio.label.col]][,1] %in% kept.labels]
    sobj.mc@misc$membership <- sobj.mc@misc$membership[sobj.mc@misc$membership %in% Cells(sobj.mc)]
    meta.data.sc <- meta.data.sc[names(sobj.mc@misc$membership),]
  } else {
    kept.labels <- unique(sobj.mc[[bio.label.col]][,1])
  }
  if (is.null(dims)) {
    dims <- 1:dim(sobj.mc[[reduction.name]])[2]
  }
  
  kept.labels <- names(which(table(sobj.mc[[bio.label.col]][,1]) > min.obs & names(table(sobj.mc[[bio.label.col]][,1])) %in% kept.labels))
  
  splitted.size <- table(interaction(sobj.mc@misc$membership,
                                     meta.data.sc[names(sobj.mc@misc$membership),bio.label.col],
                                     meta.data.sc[names(sobj.mc@misc$membership),batch.label.col], 
                                     drop = TRUE,sep = ";;"))
  
  splitted.membership <- names(splitted.size)
  
  bio <- stringr::str_split_fixed(splitted.membership,";;",n=3)[,2]
  batches <-  stringr::str_split_fixed(splitted.membership,";;",n=3)[,3]
  
  metacells <- stringr::str_split_fixed(splitted.membership,";;",n=3)[,1]
  
  splitted.embeddings <- Embeddings(sobj.mc[[reduction.name]])[metacells,dims]
  rownames(splitted.embeddings) <- splitted.membership
  
  metacell.meta_data <- data.frame("label_colnames" = batches,
                                   "split_by_colname" = bio,
                                   "size" = as.numeric(splitted.size),
                                   "size_raw" = splitted.size)
  
  # print(head(metacell.meta_data$size))
  # print(head(splitted.size))
  
  rownames(metacell.meta_data) <- splitted.membership
  
  lisi_splitByCelltype <- scIntegrationMetrics::compute_lisi_splitBy(X = splitted.embeddings, 
                                                                     meta_data = metacell.meta_data, 
                                                                     label_colnames = "label_colnames", 
                                                                     perplexity = 30, 
                                                                     split_by_colname = "split_by_colname", 
                                                                     metricsLabels = kept.labels, 
                                                                     normalize = T)
  
  res <- do.call(rbind,lisi_splitByCelltype)
  rownames(metacell.meta_data) <- paste0(metacell.meta_data$split_by_colname,".",rownames(metacell.meta_data))
  res$size <- metacell.meta_data[rownames(res),"size"]
  #head(res)
  
  ciLISI <- sum((res$label_colnames * res$size ))/sum(res$size)
  
  return(ciLISI)
}


getIntegrationMetricsMetacells <- function(sobj.mc,
                                           meta.data.sc,
                                           bio.label.col,
                                           batch.label.col,
                                           reduction.name,
                                           metricsLabels = NULL,
                                           dims = NULL,
                                           metrics = c("celltype_ASW","CiLISI")) {
  
  res <- list()
  
  if ("celltype_ASW" %in% metrics) {
    res$celltype_ASW <- MetacellBioSplitASW(sobj.mc = sobj.mc,
                                            label.col = bio.label.col,
                                            meta.data.sc =meta.data.sc,
                                            dims = dims,
                                            reduction.name = reduction.name,
                                            kept.labels = metricsLabels)
  }
  
  if ("CiLISI" %in% metrics) {
    res$CiLISI <- MetacellSplitCiLISI(sobj.mc = sobj.mc,
                                      bio.label.col = bio.label.col,
                                      batch.label.col = batch.label.col,
                                      meta.data.sc =meta.data.sc,
                                      dims = dims,
                                      kept.labels = metricsLabels,
                                      reduction.name = reduction.name)
    
  } 
  
  return(res)
  
}

# supercell_silhouette(x = sobj.mc$Sample,
#                     dist(Embeddings(sobj.mc[[reduction.name]])[,dims]), 
#                      supercell_size = rep(1,length(sobj.mc$Sample)))
# 
# cluster_sil <- cluster::silhouette(as.numeric(as.factor(sobj.mc$Sample)),dist(Embeddings(sobj.mc[[reduction.name]])[,dims]))
# mean(cluster_sil[,3])
# 
# mean(scIntegrationMetrics::compute_silhouette(as.matrix(Embeddings(sobj.mc[[reduction.name]])[,dims]),
#                                          sobj.mc@meta.data,
#                                          label_colnames = c(batch.label.col))[,1])

