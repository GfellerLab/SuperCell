#' Violin plots
#'
#' Violin plots (similar to \link[Seurat]{VlnPlot} with some changes for super-cells)
#'
#' @param ge a gene expression matrix (ncol same as number of super-cells)
#' @param supercell_size a vector with supercell size (ordered the same way as in \code{ge})
#' @param clusters a vector with clustering information (ordered the same way as in \code{ge})
#' @param features name of genes of for which gene expression is plotted
#' @param idents idents (clusters) to plot (default all)
#' @param color.use colors for idents
#' @param pt.size point size (0 by default)
#' @param pch shape of jitter dots
#' @param y.max max of y axis
#' @param y.min min of y axis
#' @param same.y.lims same y axis for all plots
#' @param adjust param of geom_violin
#' @param ncol number of columns in combined plot
#' @param combine combine plots into a single \link[patchwork]{patchwork}ed ggplot object. If FALSE, return a list of ggplot
#' @param angle.text.x rotation of x text
#' @param angle.text.y rotation of y text
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
                              y.min = NULL,
                              same.y.lims = FALSE,
                              adjust = 1,
                              ncol = NULL,
                              combine = TRUE,
                              angle.text.y = 90,
                              angle.text.x = 45){

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
  if(is.null(x = y.min)){
    y.min <- min(0,min(ge[features.set, ids.keep.idents]))
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
                                                 y.min = y.min,
                                                 adjust = adjust,
                                                 angle.text.y = angle.text.y,
                                                 angle.text.x = angle.text.x)
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
#' @param feature gene to plot
#' @param color.use colors for idents
#' @param pt.size point size (0 by default)
#' @param pch shape of jitter dots
#' @param y.max max of y axis
#' @param y.min min of y axis
#' @param adjust param of geom_violin
#' @param angle.text.x rotation of x text
#' @param angle.text.y rotation of y text
#'
#' @importFrom rlang .data
#'

supercell_VlnPlot_single <- function(ge1,
                               supercell_size = NULL,
                               clusters,
                               feature = NULL,
                               color.use = NULL,
                               pt.size = 0,
                               pch = "o",
                               y.max = NULL,
                               y.min = NULL,
                               adjust = 1,
                               angle.text.y = 90,
                               angle.text.x = 45){


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

  g <- ggplot2::ggplot(data = plot.df, ggplot2::aes(x = .data$x, y = .data$y, fill = factor(.data$group))) +
    ggplot2::geom_violin(scale="width", trim = TRUE, adjust = adjust) +
    ggplot2::scale_y_continuous(limits = c(y.min, y.max)) +
    ggplot2::labs(x = "Ident", y = "Expression Level", title = feature) +
    cowplot::theme_cowplot()

  if(pt.size > 0){
    g <- g + ggplot2::geom_jitter(data = plot.df.sc, mapping = ggplot2::aes(x = .data$x, y = .data$y, size = pt.size*.data$size),
                         pch = pch,
                         position = ggplot2::position_jitterdodge(jitter.width = 1., dodge.width = 0.25))
  }

  if(!is.null(color.use)){
    g <- g + ggplot2::scale_fill_manual(values = color.use)
  }
  g <- g + ggplot2::theme(legend.position = "none",
                 plot.title = ggplot2::element_text(hjust = 0.5),
                 axis.text.x = ggplot2::element_text(angle = angle.text.x, hjust = 1),
                 axis.text.y = ggplot2::element_text(angle = angle.text.y, hjust = 0.5))

  return(g)
}



