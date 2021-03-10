#' Violin plots
#'
#' Violin plots (similar to \link[Seurat]{VlnPlot} with some changes for super-cells)
#'
#' @param ge a gene expression matrix (ncol same as number of super-cells)
#' @param supercell_size a vector with supercell size (ordered the same way as in \code{ge})
#' @param clusters a vector with clustering information (ordered the same way as in \code{ge})
#' @param features name of geneы of for which gene expression is plotted
#' @param idents idents (clusters) to plot (default all)
#' @param color.use colors for idents
#' @param pt.size point size (0 by default)
#' @param pch shape of jitter dots
#' @param y.max max of y axis
#' @param same.y.lims same y axis for all plots
#' @param adjust param of geom_violin
#' @param ncol number of colums in combined plot
#' @param combine combine plots into a single \link[patchwork]{patchwork}ed ggplot object. If FALSE, return a list of ggplot
#'
#' @return combined ggplot or list of ggplots if \code{combine = TRUE}
#' @export

supercell_VlnPlot <- function(ge,
                              supercell_size = NULL,
                              clusters,
                              features = NULL,
                              idents = NULL,
                              color.use = NULL,
                              pt.size = 0,
                              pch = "o",
                              y.max = NULL,
                              same.y.lims = FALSE,
                              adjust = 1,
                              ncol = NULL,
                              combine = TRUE){

  N.SC <- ncol(ge) # number of super-cells

  if(is.null(clusters)) clusters <- 1

  if((length(clusters) != N.SC) & length(clusters) != 1){
    stop(paste0("clusters has to be a vector of the same lenght as ge1 (", N.SC, ") or 1, not ", length(clusters)))
  }
  if(length(clusters) == 1){
    clusters <- rep(clusters, N.SC)
  }


  if(is.null(supercell_size)) supercell_size <- rep(1, N.SC)

  if((length(supercell_size) != N.SC) & length(supercell_size) != 1){
    stop(paste0("supercell_size has to be a vector of the same lenght as ge1 (", N.SC, ") or 1, not ", length(supercell_size)))
  }
  if(length(supercell_size) == 1){
    supercell_size <- rep(supercell_size, N.SC)
  }

  if(!is.null(color.use)){
    if(length(color.use) < length(idents)){
      warning(paste0("Length of color.use (", length(color.use), ") is smaller than number of idents (",
                     length(idents),"), color.use will not be used"))
      color.use <- NULL
    }
  }

  if(is.null(idents)) idents <- sort(unique(clusters))

  ids.keep.idents <- which(clusters %in% idents)
  features.set <- features[features %in% rownames(ge)]

  if(same.y.lims & is.null(x = y.max)){
    y.max <- max(ge[features.set, ids.keep.idents])
  }


  p.list <- list()
  for(gene.i in features.set){
    ge.i <- ge[gene.i,]
    p.list[[gene.i]] <- supercell_VlnPlot_single(ge1 = ge.i[ids.keep.idents],
                                                 supercell_size = supercell_size[ids.keep.idents],
                                                 clusters = clusters[ids.keep.idents],
                                                 feature = gene.i,
                                                 color.use = color.use,
                                                 pt.size = pt.size,
                                                 pch = pch,
                                                 y.max = y.max,
                                                 adjust = adjust)
  }


  if(combine) {
    p.list <- patchwork::wrap_plots(p.list, ncol = ncol, guides = "collect")
  }
  return(p.list)
}


#' Plot  Violin plot for 1 feature
#'
#' Used for supercell_VlnPlot
#'
#' @param ge1 a gene expression vector (same length as number of super-cells)
#' @param supercell_size a vector with supercell size (ordered the same way as in \code{ge})
#' @param clusters a vector with clustering information (ordered the same way as in \code{ge})
#' @param color.use colors for idents
#' @param pt.size point size (0 by default)
#' @param pch shape of jitter dots
#' @param y.max max of y axis
#' @param adjust param of geom_violin
#'
#' @importFrom ggplot2 ggplot aes_string geom_point scale_size scale_radius geom_violin scale_color_manual
#' theme element_blank labs scale_color_identity scale_color_distiller geom_jitter aes
#' scale_color_gradient guides guide_legend guide_colorbar scale_y_continuous scale_x_continuous element_text
#' geom_text margin
#'

supercell_VlnPlot_single <- function(ge1,
                               supercell_size = NULL,
                               clusters,
                               feature = NULL,
                               color.use = NULL,
                               pt.size = 0,
                               pch = "o",
                               y.max = NULL,
                               adjust = 1){


  N.SC <- length(ge1)

  sc.expansion <- rep(1:length(clusters), supercell_size)
  plot.df <- data.frame(x = clusters[sc.expansion],
                        y = ge1[sc.expansion],
                        group = clusters[sc.expansion],
                        size = supercell_size[sc.expansion])

  plot.df.sc <- data.frame(x = clusters,
                           y = ge1,
                           group = clusters,
                           size = supercell_size)

  if(is.null(y.max)) y.max <- max(ge1)

  g <- ggplot(data = plot.df, aes(x = x, y = y, fill = factor(group))) +
    geom_violin(scale="width", trim = TRUE, adjust = adjust) +
    scale_y_continuous(limits = c(0, y.max)) +
    labs(x = "Ident", y = "Expression Level", title = feature) +
    cowplot::theme_cowplot()

  if(pt.size > 0){
    g <- g + geom_jitter(data = plot.df.sc, mapping = aes(x = x, y = y, size = pt.size*size),
                         pch = pch,
                         position = position_jitterdodge(jitter.width = 1., dodge.width = 0.25))
  }

  if(!is.null(color.use)){
    g <- g + scale_fill_manual(values = color.use)
  }
  g <- g + theme(legend.position = "none",
                 plot.title = element_text(hjust = 0.5),
                 axis.text.x = element_text(angle = 45, hjust = 1))

  return(g)
}



