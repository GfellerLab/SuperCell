#' Plot super-cell NW
#'
#' @param SC.nw a super-cell network (a field \code{supercell_network} of the output of \link{SCimplify})
#' @param group an assigment of super cell to any group (for ploting in different colors)
#' @param color.use colros to use for groups, if \code{NULL}, an automatic palette of colors will be applied
#' @param lay.method method to compute layout of the network (for the moment there several available: "nicely"
#' for \link[igraph]{layout_licely} and "fr" for \link[igraph]{layout_with_fr}, "components" for \link[igraph]{layout_components})
#' @param lay a particular layout of a graph to plot (in is not \code{NULL}, \code{lay.method} is ignored and new layout is not computed)
#' @param seed a random seed used to compute graph layout
#' @param main a title of a plot
#' @param do.frames whether to keep vertex.frames in the plot
#'
#' @return plot of a super-cell network
#' @export


supercell_plot <- function(SC.nw, group = NULL, color.use = NULL, lay.method = c("nicely", "fr", "components"), lay = NULL, seed = 12345, main = NA, do.frames = TRUE){
  N.SC <- igraph::vcount(SC.nw)
  if(is.null(igraph::V(SC.nw)$size)) igraph::V(SC.nw)$size <- 1
  vsize <- sqrt(igraph::V(SC.nw)$size)

  if(is.null(group)) group <- rep(1, N.SC)
  if(length(group) != N.SC) stop(paste("Vector groups has to be the same length as number of super-cells:", N.SC))

  N.groups <- length(unique(group))

  if(is.null(color.use)) color.use <- scales::hue_pal()(N.groups)
  if(length(color.use) != N.groups) stop(paste("Vector color.use has to be the same length as number of groups:", N.groups))

  names(color.use) <- sort(unique(group))


  if(is.null(lay)){ # compute new layout
    lay.method <- lay.method[1]
    if(is.na(lay.method) | is.nan(lay.method) | is.null(lay.method)) lay.method = "nicely"

    lay.method <- match(lay.method, c("nicely", "fr", "components"))

    if(is.na(lay.method)) stop(paste("Unknown layout!"))

    set.seed(seed)
    if(lay.method == 1)
      lay <- igraph::layout_nicely(SC.nw)
    else if (lay.method == 2)
      lay <- igraph::layout_with_fr(SC.nw)
    else if (lay.method == 3)
      lay <- igraph::layout_components(SC.nw)
    else stop(paste("Unknown lay.method", lay.method))
  }


  v.colors <- color.use[group]
  v.frame.colors <- ifelse(do.frames, "black", v.colors)
  p <- plot(SC.nw, vertex.label = NA, vertex.color = v.colors, vertex.frame.color = v.colors,
            layout = lay, main = main, vertex.size = vsize)

}
