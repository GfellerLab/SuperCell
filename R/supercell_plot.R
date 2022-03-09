#' Plot super-cell NW
#'
#' @param SC.nw a super-cell network (a field \code{supercell_network} of the output of \link{SCimplify})
#' @param group an assigment of super cell to any group (for ploting in different colors)
#' @param color.use colros to use for groups, if \code{NULL}, an automatic palette of colors will be applied
#' @param lay.method method to compute layout of the network (for the moment there several available: "nicely"
#' for \link[igraph]{layout_licely} and "fr" for \link[igraph]{layout_with_fr}, "components" for \link[igraph]{layout_components},
#' "drl" for \link[igraph]{layout_with_drl}, "graphopt" for \link[igraph]{layout_with_graphopt}). If your dataset has clear clusters, use "components"
#' @param lay a particular layout of a graph to plot (in is not \code{NULL}, \code{lay.method} is ignored and new layout is not computed)
#' @param alpha a rotation of the layout (either provided or computed)
#' @param seed a random seed used to compute graph layout
#' @param main a title of a plot
#' @param do.frames whether to keep vertex.frames in the plot
#' @param do.extra.log.rescale whether to log-scale node size (to balance plot if some super-cells are large and covers smaller super-cells)
#' @param do.directed whether to plot edge direction
#' @param log.base base with thich to log-scale node size
#' @param do.extra.sqtr.rescale  whether to sqrt-scale node size (to balance plot if some super-cells are large and covers smaller super-cells)
#' @param frame.color color of node frames, black by default
#' @param weights edge weights used for some layout algorithms
#' @param min.cell.size do not plot cells with smaller size
#' @param return.meta whether to return all the meta data
#'
#' @return plot of a super-cell network
#'
#' @examples
#' \dontrun{
#' data(cell_lines) # list with GE - gene expression matrix (logcounts), meta - cell meta data
#' GE <- cell_lines$GE
#' cell.meta <- cell_lines$meta
#'
#' SC <- SCimplify(GE,  # gene expression matrix
#'                 gamma = 20) # graining level
#'
#' # Assign super-cell to a cell line
#' SC2cellline  <- supercell_assign(clusters = cell.meta, # single-cell assigment to cell lines
#'                                  supercell_membership = SC$membership) # single-cell assignment to super-cells
#'
#' # Plot super-cell network colored by cell line
#' supercell_plot(SC$graph.supercells, # network
#'                group = SC2cellline, # group assignment
#'                main = "Super-cell colored by cell line assignment",
#'                lay.method = 'nicely')
#' }
#'
#' @export


supercell_plot <- function(SC.nw,
                           group = NULL,
                           color.use = NULL,
                           lay.method = c("nicely", "fr", "components", "drl", "graphopt"),
                           lay = NULL,
                           alpha = 0,
                           seed = 12345,
                           main = NA,
                           do.frames = TRUE,
                           do.extra.log.rescale = FALSE,
                           do.directed = FALSE,
                           log.base = 2,
                           do.extra.sqtr.rescale = FALSE,
                           frame.color = "black",
                           weights = NULL,
                           min.cell.size = 0,
                           return.meta = FALSE){



  N.SC <- igraph::vcount(SC.nw)
  if(is.null(igraph::V(SC.nw)$size)) igraph::V(SC.nw)$size <- 1
  vsize <- igraph::V(SC.nw)$size

  cells.to.remove <- which(vsize < min.cell.size)
  cells.to.use    <- setdiff(1:N.SC, cells.to.remove)

  vsize <- sqrt(vsize)
  if(do.extra.log.rescale)
    vsize <- log(vsize+1, base = log.base)
  if(do.extra.sqtr.rescale)
    vsize <- sqrt(vsize)

  if(is.null(group)) group <- rep(1, N.SC)
  if(is.numeric(group)){
    group <- group + (1 - min(group)) # min group to 1
  }
  if(!is.character(group)) group <- as.character(group)
  if(length(group) != N.SC) stop(paste("Vector groups has to be the same length as number of super-cells:", N.SC))

  N.groups <- length(unique(group))

  if(is.null(color.use)) color.use <- scales::hue_pal()(N.groups)
  if(length(color.use) != N.groups) stop(paste("Vector color.use has to be the same length as number of groups:", N.groups, "\n"))

  if(is.null(names(color.use))){
    names(color.use) <- sort(unique(group))
  }

  if(!is.null(weights) & length(cells.to.remove)>0){
    weights <- weights[-sort(unique(unlist(igraph::incident_edges(SC.nw, cells.to.remove))))]
  }

  if(length(cells.to.remove)>0){
    SC.nw <- igraph::delete_vertices(SC.nw, cells.to.remove)
  }



  if(is.null(lay)){ # compute new layout
    lay.method <- lay.method[1]
    if(is.na(lay.method) | is.nan(lay.method) | is.null(lay.method)) lay.method = "nicely"

    lay.method <- match(lay.method,  c("nicely", "fr", "components", "drl", "graphopt", "dh"))

    if(is.na(lay.method)) stop(paste("Unknown layout!\n"))

    if(lay.method == "components" & igraph::vcount(SC.nw) > 5000){
      warning("Number of super-cells is too large for the 'components' layout, lay.method was changed to 'nicely'!\n")
    }

    set.seed(seed)
    if(lay.method == 1)
      lay <- igraph::layout_nicely(SC.nw, weights = weights)
    else if (lay.method == 2)
      lay <- igraph::layout_with_fr(SC.nw, weights = weights)
    else if (lay.method == 3)
      lay <- igraph::layout_components(SC.nw)
    else if (lay.method == 4)
      lay <- igraph::layout_with_drl(SC.nw, weights = weights)
    else if (lay.method == 5)
      lay <- igraph::layout_with_graphopt(SC.nw)
    else stop(paste("Unknown lay.method", lay.method, "\n"))
  }

  # do rotation
  if(is.numeric(alpha)){
    if((alpha %% (2*pi)) != 0){
      M.rotation  <- matrix(c(cos(alpha), sin(alpha),
                              -sin(alpha), cos(alpha)), byrow = T, ncol = 2)
      lay.rotated <- lay %*% M.rotation
    } else {
      lay.rotated <- as.matrix(lay)
    }
  } else {
    stop(paste("Alpha has to be numeric! \n"))
  }

  v.colors <- color.use[group]
  ifelse(do.frames, v.frame.colors <- rep(frame.color, length(v.colors)), v.frame.colors <- v.colors)

  plot(igraph::as.undirected(SC.nw), vertex.label = NA, vertex.color = v.colors[cells.to.use],
       vertex.frame.color = v.frame.colors[cells.to.use],
       layout = lay.rotated, main = main, vertex.size = vsize[cells.to.use])

  if(return.meta){
    p <- grDevices::recordPlot()
    res <- list(p = p,
                edgelist = edgelist,
                lay = lay,
                lay.rotated = lay.rotated,
                alpha = alpha,
                group = group[cells.to.use],
                cells.to.use = cells.to.use,
                vsize.used = vsize[cells.to.use],
                vsize.actual = igraph::V(SC.nw)$size,
                v.colors = v.colors)
    return(res)
  }
}


