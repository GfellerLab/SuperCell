
#' Weighted mean function
#' @export
#'
weighted.mean.fxn <- function (x, weights) {
  rowSums(x %*% Matrix::Diagonal(x = weights))/sum(weights)
}

#' Multimodal markers for all identity classes
#'
#' Function identifying markers in 2 modalities using FindAllMarkers.SuperCell
#'
#'
#' @param object Seurat metacell object with a metadata size column
#' @param assay1 Seurat assay to consider as the first modality
#' @param assay2 Seurat assay to consider as the second modality
#' @param test.use (default is "survey_weighted_t") Test to use in FindAllMarkers.SuperCell
#' @param fc.name1 fc.name parameter that FindAllMarkers.SuperCell will consider for the first modality
#' @param fc.name2 fc.name parameter that FindAllMarkers.SuperCell will consider for the second modality
#' @param logfc.threshold1 logfc.threshold parameter that FindAllMarkers.SuperCell will consider for the first modality
#' @param logfc.threshold2 logfc.threshold parameter that FindAllMarkers.SuperCell will consider for the second modality
#' @param mean.fxn1 mean.fxn parameter that FindAllMarkers.SuperCell will consider for the first modality
#' @param mean.fxn2 mean.fxn parameter that FindAllMarkers.SuperCell will consider for the second modality
#' @inheritParams FindAllMarkers.SuperCell
#' @return A list of 3 data.frames: (i) a data frame containing a ranked list of putative multimodal markers as rows, and associated statistics as columns (p-values, t statistic, degree of freedom (df) of the weighted t-test).
#' avg_logFC, pct.1, pc.2 are computed taking into account metacell sizes, (ii) a data frame containing the output of FindAllMarkers.SuperCell run on the first modality, and (iii) a data frame containing the output of FindAllMarkers.SuperCell run onthe second modality.
#' @import Seurat
#' @export

FindMultimodalMarkers.SuperCell <- function(seurat.obj, group.by = 'celltype',
                                 assay1 = 'RNA', assay2 = "chromvar",
                                 min.cells.feature = 0, min.cells.group = 0, min.pct = 0.01,
                                 padj.cutoff = 0.05, base = 2, only.pos = T, return.thresh = 1,
                                 logfc.threshold1 = 0.1, logfc.threshold2 = 0.1,
                                 test.use = "survey_weighted_t",
                                 fc.name1 = "avg_log2FC", fc.name2 = "avg_diff",
                                 mean.fxn1 = NULL, mean.fxn2 = weighted.mean.fxn, ...) {

  Idents(seurat.obj) <- group.by
  markers_mod1 <- FindAllMarkers.SuperCell(
    object = seurat.obj, assay = assay1,
    min.cells.feature = min.cells.feature,
    min.cells.group = min.cells.group,
    base = base,
    min.pct = min.pct, fc.name = fc.name1,
    logfc.threshold = logfc.threshold1,
    only.pos = only.pos, mean.fxn = mean.fxn1, test.use = test.use,
    return.thresh = return.thresh, ...
  )
  colnames(markers_mod1) <- paste0(assay1, ".", colnames(markers_mod1))

  DefaultAssay(seurat.obj) <- assay2

  markers_mod2 <- FindAllMarkers.SuperCell(
    object = seurat.obj,
    assay = assay2,
    min.cells.feature = min.cells.feature,
    min.cells.group = min.cells.group,
    base = base,
    min.pct = min.pct,
    logfc.threshold = logfc.threshold2, test.use = test.use,
    only.pos = only.pos, fc.name = fc.name2, mean.fxn = weighted.mean.fxn, #SuperCell:::weighted.mean.fxn
    return.thresh = return.thresh, ...
  )
  colnames(markers_mod2) <- paste0(assay2, ".", colnames(markers_mod2))

  # markers_mod1$gene <- markers_mod1[, paste0(assay1, ".gene")]

  if(assay1 == "chromvar"){
    DefaultAssay(seurat.obj) <- "ATAC"
    markers_mod1$gene <- ConvertMotifID(seurat.obj, id = markers_mod1[, paste0(assay1, ".gene")])
  }else{
    markers_mod1$gene <- markers_mod1[, paste0(assay1, ".gene")]
  }


  if(assay2 == "chromvar"){
    DefaultAssay(seurat.obj) <- "ATAC"
    markers_mod2$gene <- ConvertMotifID(seurat.obj, id = markers_mod2[, paste0(assay2, ".gene")])
  }else{
    markers_mod2$gene <- markers_mod2[, paste0(assay2, ".gene")]
  }


  markers.all <- vector()
  for(celltype in unique(seurat.obj@meta.data[,group.by])){
    # print(celltype)
    ctmarkers_mod1 <- dplyr::filter(
      markers_mod1,
      !!dplyr::sym(paste0(assay1, ".cluster")) == celltype,
      !!dplyr::sym(paste0(assay1, ".p_val_adj")) < padj.cutoff,
      !!dplyr::sym(paste0(assay1, ".", fc.name1)) > 0) %>%
      dplyr::arrange(-!!dplyr::sym(paste0(assay1, ".", fc.name1)))

    ctmarkers_mod2 <- dplyr::filter(
      markers_mod2,
      !!dplyr::sym(paste0(assay2, ".cluster")) == celltype,
      !!dplyr::sym(paste0(assay2, ".p_val_adj"))  < padj.cutoff,
      !!dplyr::sym(paste0(assay2, ".", fc.name2)) > 0) %>%
      dplyr::arrange(-!!dplyr::sym(paste0(assay2, ".", fc.name2)))

    multimodal.markers <- dplyr::inner_join(
      x = ctmarkers_mod1,
      y = ctmarkers_mod2,
      by = "gene"
    )

    X <- as(GetAssayData(seurat.obj, assay = assay1), "dgCMatrix")
    y <- factor(ifelse(seurat.obj@meta.data[, group.by] == celltype, celltype, "other"))
    group.size <- as.numeric(table(y))
    n1n2 <- group.size * (ncol(X) - group.size)
    rank_res <- presto::rank_matrix(Matrix::t(X[multimodal.markers[, paste0(assay1, ".gene")],]))
    ustat <- presto:::compute_ustat(rank_res$X_ranked, y, n1n2, group.size)
    auc.rna <- t(ustat/n1n2)[,1]

    X <- as(GetAssayData(seurat.obj, assay = assay2), "dgCMatrix")
    rank_res <- presto::rank_matrix(Matrix::t(X[multimodal.markers[, paste0(assay2, ".gene")],]))
    ustat <- presto:::compute_ustat(rank_res$X_ranked, y, n1n2, group.size)
    auc.motif <- t(ustat/n1n2)[,1]

    multimodal.markers[, paste0(assay1, ".auc")] <- auc.rna
    multimodal.markers[, paste0(assay2, ".auc")]<- auc.motif #unlist(lapply(auc.res, function(x) x[,"motif.auc"]))
    multimodal.markers$auc.mean <- rowMeans(data.frame(auc.rna, auc.motif))

    multimodal.markers[, paste0(assay1, ".rpb")] <- multimodal.markers[, paste0(assay1, ".t_value")] / sqrt(multimodal.markers[, paste0(assay1, ".t_value")]^2 + multimodal.markers[, paste0(assay1, ".df")])
    multimodal.markers[, paste0(assay2, ".rpb")] <- multimodal.markers[, paste0(assay2, ".t_value")] / sqrt(multimodal.markers[, paste0(assay2, ".t_value")]^2 + multimodal.markers[, paste0(assay2, ".df")])
    multimodal.markers$mean.rpb <- rowMeans(multimodal.markers[, c(paste0(assay1, ".rpb"), paste0(assay2, ".rpb"))])

    multimodal.markers <- dplyr::arrange(multimodal.markers, -mean.rpb)

    if(nrow(multimodal.markers) != 0){
      markers.all <- rbind(markers.all, multimodal.markers)
    }
  }

  return(
    setNames(
      list(
        markers.all,
        markers_mod1,
        markers_mod2
      ), nm = c("MultimodalMarkers", paste0("markers.", assay1), paste0("markers.", assay2))
    )
  )
}
