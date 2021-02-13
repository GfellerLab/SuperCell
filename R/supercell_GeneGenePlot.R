#' Gene-gene correlation plot
#'
#' Plots gene-gene expression and computes their correaltion
#'
#' @param ge a gene expression matrix of super-cells (ncol same as number of super-cells)
#' @param gene_x gene or vector of genes (if vector, has to be the same lenght as gene_y)
#' @param gene_y gene or vector of genes (if vector, has to be the same lenght as gene_x)
#' @param supercell_size a vector with supercell size (ordered the same way as in \code{ge})
#' @param clusters a vector with clustering information (ordered the same way as in \code{ge})
#' @param color.use colors for idents
#' @param idents idents (clusters) to plot (default all)
#' @param pt.size point size (if supercells have identical sizes)
#' @param alpha transparency
#' @param x.max max of x axis
#' @param y.max max of y axis
#' @param same.x.lims same x axis for all plots
#' @param same.y.lims same y axis for all plots
#' @param ncol number of colums in combined plot
#' @param combine combine plots into a single \link[patchwork]{patchwork}ed ggplot object. If FALSE, return a list of ggplot
#' @param sort.by.corr whether to sort plots by absolute value of correlation (fist plot genes with largest (anti-)correlation)
#'
#'
#'@return a list with components
#' \itemize{
#'   \item p - is a combined ggplot or list of ggplots if \code{combine = TRUE}
#'   \item w.cor - weighted correlation between genes
#' }
#' @return a list, where
#'
#' @export

supercell_GeneGenePlot <- function(ge,
                                   gene_x,
                                   gene_y,
                                   supercell_size = NULL,
                                   clusters = NULL,
                                   color.use = NULL,
                                   idents = NULL,
                                   pt.size = 1,
                                   alpha = 0.9,
                                   x.max = NULL,
                                   y.max = NULL,
                                   same.x.lims = FALSE,
                                   same.y.lims = FALSE,
                                   ncol = NULL,
                                   combine = TRUE,
                                   sort.by.corr = TRUE){

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


  if(is.null(idents)) idents <- sort(unique(clusters))

  ids.keep.idents <- which(clusters %in% idents)


  if(!is.null(color.use)){
    if(length(color.use) < length(idents)){
      warning(paste0("Length of color.use (", length(color.use), ") is smaller than number of idents (",
                     length(idents),"), color.use will not be used"))
      color.use <- NULL
    }
  }


  if(length(gene_x) != length(gene_y)){
    if(length(gene_x) == 1){
      gene_x <- rep(gene_x, length(gene_y))
    } else {
      if(length(gene_y) == 1){
        gene_y <- rep(gene_y, length(gene_x))
      } else{
        stop("Vectors gene_x and gene_y need to have the same length or one of them has to be a vector of length 1")
      }
    }
  }

  # keep genes that are present in the gene expression dataset
  features.set_x <- gene_x[gene_x %in% rownames(ge) & gene_y %in% rownames(ge)]
  features.set_y <- gene_y[gene_x %in% rownames(ge) & gene_y %in% rownames(ge)]


  if(same.y.lims & is.null(x = y.max)){
    y.max <- max(ge[features.set_y, ids.keep.idents])
  }

  if(same.x.lims & is.null(x.max)){
    x.max <- max(ge[features.set_x, ids.keep.idents])
  }

  p.list <- list()



  for(i in 1:length(features.set_x)){

    genes.i <- paste(features.set_x[i], features.set_y[i], sep = "_")

    p.list[[genes.i]] <- supercell_GeneGenePlot_single(ge_x = ge[features.set_x[i], ids.keep.idents],
                                                       ge_y = ge[features.set_y[i], ids.keep.idents],
                                                       gene_x_name = features.set_x[i],
                                                       gene_y_name = features.set_y[i],
                                                       supercell_size = supercell_size[ids.keep.idents],
                                                       clusters = clusters[ids.keep.idents],
                                                       color.use = color.use,
                                                       x.max = x.max,
                                                       y.max = y.max,
                                                       pt.size = pt.size,
                                                       alpha = alpha)
  }




  if(sort.by.corr){ # sort plots by absolute value of correkation
    p.list <- p.list[names(sort(abs(unlist(lapply(p.list, FUN = function(x){x$w.cor[1]}))), decreasing = T))]
  }


  w.cor.list <- lapply(p.list, FUN = function(x){x$w.cor})
  p.list <-lapply(p.list, FUN = function(x){x$g})

  if(combine) {
    p.list <- patchwork::wrap_plots(p.list, ncol = ncol, guides = "collect")
  }
  return(list(p = p.list, w.cor = w.cor.list))
}







#' Plot  Gene-gene correlation plot for 1 feature
#'
#' Used for \link{supercell_GeneGenePlot}
#'
#' @param ge_x first gene expression vector (same length as number of super-cells)
#' @param ge_y second gene expression vector (same length as number of super-cells)
#' @param gene_x_name name of gene x
#' @param gene_y_name name of gene y
#' @param supercell_size a vector with supercell size (ordered the same way as in \code{ge})
#' @param clusters a vector with clustering information (ordered the same way as in \code{ge})
#' @param color.use colors for idents
#' @param x.max max of x axis
#' @param y.max max of y axis
#' @param pt.size point size (0 by default)
#' @param alpha transparency of dots
#'
#' @importFrom ggplot2 ggplot aes_string geom_point scale_size scale_radius geom_violin scale_color_manual
#' theme element_blank labs scale_color_identity scale_color_distiller geom_jitter aes
#' scale_color_gradient guides guide_legend guide_colorbar scale_y_continuous

supercell_GeneGenePlot_single <- function(ge_x,
                                          ge_y,
                                          gene_x_name,
                                          gene_y_name,
                                          supercell_size = NULL,
                                          clusters = NULL,
                                          color.use = NULL,
                                          x.max = NULL,
                                          y.max = NULL,
                                          pt.size = 1,
                                          alpha = 0.9){



  N.SC <- length(ge_x)

  plot.df <- data.frame(x = ge_x,
                        y = ge_y,
                        identity = factor(clusters),
                        size = supercell_size)


  w.cor   <- weights::wtd.cor(plot.df$x, plot.df$y, weight = plot.df$size)


  if(is.null(x.max)) x.max <- NA
  if(is.null(y.max)) y.max <- NA


  if(length(unique(plot.df$size)) == 1){

     g <- ggplot(plot.df, aes(x = x, y = y, color = identity)) +
       geom_point(alpha = alpha, size = pt.size)
  } else {
    g <- ggplot(plot.df, aes(x = x, y = y, color = identity, size = size)) +
      geom_point(alpha = alpha)
  }

  g <- g +
    scale_x_continuous(limits = c(0, x.max)) +
    scale_y_continuous(limits = c(0, y.max)) +
    theme_classic() + theme(asp = 1) + #, legend.position = "none"
    labs(x = gene_x_name,
         y = gene_y_name,
         title = paste0("\nw_cor = ", signif(w.cor[1], 2), ", pval = ", signif(w.cor[4])))


  if(!is.null(color.use)){
    g <- g + scale_color_manual(values = color.use)
  }

  res <- list(g = g, w.cor = w.cor)
  return(res)
}