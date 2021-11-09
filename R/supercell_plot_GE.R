#' Plot super-cell NW colored by an expression of a gene (gradient color)
#'
#' @param SC.nw a super-cell network (a field \code{supercell_network} of the output of \link{SCimplify})
#' @param ge a gene expression vector (same length as number of super-cells)
#' @param color.use colors of gradient
#' @param n.color.gradient number of bins of the gradient, default is 10
#' @param legend.side a side parameter of \link[plotfunctions]{gradientLegend} fucntion (default is 4)
#' @param gene.name name of gene of for which gene expression is plotted
#' @param ... rest of the parameters of \link{supercell_plot} function

#'
#' @return plot of a super-cell network with color representing an expression level
#' @export

supercell_plot_GE <- function(SC.nw,
                              ge,
                              color.use = c("gray", "blue"),
                              n.color.gradient = 10,
                              lay.method = c("nicely", "fr", "components", "drl", "graphopt", "dh"),
                              lay = NULL,
                              alpha = 0,
                              seed = 12345,
                              main = NA,
                              do.frames = TRUE,
                              do.extra.log.rescale = FALSE,
                              log.base = 2,
                              do.extra.sqtr.rescale = FALSE,
                              frame.color = "black",
                              weights = NULL,
                              min.cell.size = 0,
                              legend.side = 4,
                              gene.name = NULL){

  if(length(ge) != igraph::vcount(SC.nw)){
    stop("ge must be the same length as number of super-cells")
  }

  if(is.null(color.use[1]) | is.na(color.use[1]) | length(color.use) < 2){
    color.use <- c("gray", "blue")
  }

  areColors <- function(x) {
    sapply(x, function(X) {
      tryCatch(is.matrix(col2rgb(X)),
               error = function(e) FALSE)
    })
  }

  if(sum(areColors(color.use)) != length(areColors(color.use))){
    stop(paste("color.use contains non color variables:", paste(color.use[!areColors(color.use)], collapse = ", ")))
  }


  if(is.null(n.color.gradient) | is.na(n.color.gradient) | is.nan(n.color.gradient)){
    n.color.gradient <- 10
  }
  if(n.color.gradient < 3)
    warning(paste0("Color gradient length is very small ", n.color.gradient ))

  color.gradient <- grDevices::colorRampPalette(color.use)(n.color.gradient)
  names(color.gradient) <- as.character(1:length(color.gradient))

  brks           <- seq(floor(min(c(ge, 0))), ceiling(max(ge)), length.out = n.color.gradient+1)
  #print(brks)
  if(length(unique(brks)) == length(brks)){
    group          <- as.character(cut(ge, breaks=brks,  include.lowest=TRUE, labels = FALSE))
  } else {
    group <- as.character(rep(1, length(ge)))
  }

  p <- supercell_plot(SC.nw = SC.nw,
                   group = group,
                   color.use = color.gradient[unique(group)],
                   lay.method = lay.method,
                   lay = lay,
                   alpha = alpha,
                   seed = seed,
                   main = paste(main, gene.name),
                   do.frames = do.frames,
                   do.extra.log.rescale = do.extra.log.rescale,
                   log.base = log.base,
                   do.extra.sqtr.rescale = do.extra.sqtr.rescale,
                   frame.color = frame.color,
                   weights = weights,
                   min.cell.size = min.cell.size,
                   return.meta = TRUE)


  p$p
  plotfunctions::gradientLegend(valRange = c(min(brks), max(brks)), color = color.gradient, side = legend.side)

  p.res <- grDevices::recordPlot()

  return(p.res)
}



