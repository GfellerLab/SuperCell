
#' FeatureFeaturePlot from SuperCell seurat object
#'
#' FeatureFeaturePlot function using metacell sizes
#'
#' @param seurat.obj A Seurat object or Seurat metacell object
#' @param assays (default is c("RNA", "RNA")) Vector of 2 assays from which `feature.x` and `feature.y` will be selected
#' @param feature.x Vector of features from the first assay to represent on the x axis
#' @param feature.y Vector of features from the second assay to represent on the y axis
#' @param cluster Names of the meta.data column to use to color the cells or metacells
#' @param is.normalized (default is FALSE) Bolean indicating whether to normalize the assays in `assays`
#' @param normalization.method.1 Method for normalization for the first assay in `assays`
#' @param normalization.method.2 Method for normalization for the second assay in `assays`
#' @param norm.margin.1 (default is 1) If performing CLR normalization for the first assay in `assays`, normalize across features (1) or cells (2)
#' @param norm.margin.2 (default is 1) If performing CLR normalization for the second assay in `assays`, normalize across features (1) or cells (2)
#' @param use.size Consider metacell sizes in the scatter plot (dot sizes) and in the correlation test between the features in feature.x and feature.y
#' @param plot (default is FALSE) Bolean indicating whether to generate and save the scatter plots or only retrieve correlation values
#' @return A list of SingleCorPlot scatter plot for single-cells or metacells
#' @import Seurat
#' @export

FeatureFeaturePlot.SuperCell <- function (seurat.obj, feature.x, feature.y, cluster = NULL,
                                          assays = c("RNA", "RNA"), #nb.cores = 2,
                                          norm.margin.1 = 1, normalization.method.1 = "LogNormalize",
                                          norm.margin.2 = 1, normalization.method.2 = "LogNormalize",
                                          is.normalized = F, plot = F, color.use = NULL, use.size = T, add.key = T,...)
{
  Seurat::DefaultAssay(seurat.obj) <- assays[1]
  if (!is.normalized) {
    seurat.obj <- Seurat::NormalizeData(seurat.obj, normalization.method = normalization.method.1,
                                        margin = norm.margin.1)
  }
  fe1 <- Seurat::GetAssayData(seurat.obj, slot = "data", assay = assays[1])[feature.x, , drop= F ]
  feature.x <- paste0(gsub("_", "", tolower(assays[1])), "_", feature.x)
  rownames(fe1) <- feature.x

  Seurat::DefaultAssay(seurat.obj) <- assays[2]
  if (!is.normalized) {
    seurat.obj <- Seurat::NormalizeData(seurat.obj, normalization.method = normalization.method.2,
                                        margin = norm.margin.2)
  }
  fe2 <- Seurat::GetAssayData(seurat.obj, slot = "data", assay = assays[2])[feature.y, , drop= F]
  feature.y <- paste0(gsub("_", "", tolower(assays[2])), "_", feature.y)
  rownames(fe2) <- feature.y
  fe <- rbind(fe1, fe2)

  if (use.size) {
    sizes <- as.numeric(seurat.obj$size)
  } else {
    sizes <- rep(1, ncol(seurat.obj))
  }

  res <- list()
  cor.res <- lapply(1:length(feature.x), function(i) {
    tryCatch({
      weights::wtd.cor(fe[feature.x[i], ], fe[feature.y[i], ], weight = sizes)
    }, error = function(e) {
      data.frame(correlation = NA, p.value = NA)
    })
  })

  res[["w.cor"]] <- lapply(cor.res, function(x) x[,"correlation"])
  res[["w.pval"]] <- lapply(cor.res, function(x) x[, "p.value"])
  if(plot){
    seurat.obj$size <- sizes
    res[["p"]] <- lapply(1:length(feature.x), function(i){
      if(add.key){
        FeatureScatter.SuperCell(object = seurat.obj,
                                 size.by = "size",
                                 feature1 = feature.x[i],
                                 feature2 = feature.y[i],
                                 group.by = cluster,
                                 cols = color.use,
                                 plot.cor = plot,...)
      }else{
        FeatureScatter.SuperCell(object = seurat.obj,
                                 size.by = "size",
                                 feature1 = feature.x[i],
                                 feature2 = feature.y[i],
                                 group.by = cluster,
                                 cols = color.use,
                                 plot.cor = plot,...) +
          xlab(gsub(paste0(gsub("_", "", tolower(assays[1])),"_"), "", feature.x[i])) +
          ylab(gsub(paste0(gsub("_", "", tolower(assays[2])),"_"), "", feature.y[i]))
      }
    } )
    names(res[["p"]]) <- feature.x
  }

  w.cor <- data.frame(row.names = paste(feature.x, feature.y, sep = "_"),
                      feature1 = feature.x,
                      feature2 = feature.y,
                      w.cor = as.numeric(res$w.cor),
                      w.pval = as.numeric(res$w.pval),
                      w.qval = p.adjust(as.numeric(res$w.pval), method = "bonferroni"))

  if (plot) {
    return(list(cor.res = w.cor, plots = res$p))
  }else{
    return(w.cor)
  }
}

#' FeatureScatter plot with SuperCell
#'
#' FeatureScatter Seurat function using metacell sizes
#'
#' @param object A Seurat metacell object
#' @inheritParams Seurat::FeatureScatter
#' @return A FeatureScatter sctter plot for metacells
#' @import Seurat
#' @export

FeatureScatter.SuperCell <- function (object, feature1, feature2, cells = NULL, shuffle = FALSE,
                                      seed = 1, group.by = NULL, split.by = NULL, cols = NULL,
                                      size.by = "size", pt.size = 0.5, shape.by = NULL, span = NULL,
                                      smooth = FALSE, combine = TRUE, slot = "data", plot.cor = TRUE,
                                      ncol = NULL, raster = NULL, raster.dpi = c(512, 512), jitter = FALSE,
                                      log = FALSE)
{
  cells <- cells %||% colnames(x = object)
  if (isTRUE(x = shuffle)) {
    set.seed(seed = seed)
    cells <- sample(x = cells)
  }
  group.by <- group.by %||% "ident"
  data <- Seurat::FetchData(object = object, vars = c(feature1, feature2, size.by, group.by), cells = cells, slot = slot)
  if (!grepl(pattern = feature1, x = names(x = data)[1])) {
    rlang::abort(message = paste("Feature 1", sQuote(x = feature1),
                          "not found"))
  }
  if (!grepl(pattern = feature2, x = names(x = data)[2])) {
    rlang::abort(message = paste("Feature 2", sQuote(x = feature2),
                          "not found"))
  }
  feature1 <- names(x = data)[1]
  feature2 <- names(x = data)[2]
  group.by <- intersect(x = group.by, y = names(x = data)[3:ncol(x = data)])
  for (group in group.by) {
    if (!is.factor(x = data[, group])) {
      data[, group] <- factor(x = data[, group])
    }
  }
  if (!is.null(x = split.by)) {
    split <- Seurat::FetchData(object = object, vars = split.by,
                               clean = TRUE)[split.by]
    data <- data[rownames(split), ]
    data[, split.by] <- split
  }
  plots <- lapply(X = group.by, FUN = function(x) {
    plot <- SingleCorPlot.SuperCell(data = data[, c(feature1, feature2, split.by, size.by)], col.by = data[, x],
                                    cols = cols, pt.size = pt.size, smooth = smooth, size.by = size.by,
                                    legend.title = "Identity", span = span, plot.cor = plot.cor,
                                    raster = raster, raster.dpi = raster.dpi, jitter = jitter)
    if (!is.null(x = split.by)) {
      plot <- plot + Seurat:::FacetTheme() + ggplot2::facet_wrap(facets = dplyr::vars(!!rlang::sym(x = split.by)),
                                                                 ncol = if (length(x = group.by) > 1 || is.null(x = ncol)) {
                                                                   length(x = unique(x = data[, split.by]))
                                                                 }
                                                                 else {
                                                                   ncol
                                                                 })
    }
    if (log) {
      plot <- plot + ggplot2::scale_x_log10() + ggplot2::scale_y_log10()
    }
    plot
  })
  if (isTRUE(x = length(x = plots) == 1)) {
    return(plots[[1]])
  }
  if (isTRUE(x = combine)) {
    plots <- patchwork::wrap_plots(plots, ncol = length(x = group.by))
  }
  return(plots)
}

#' SingleCorPlot plot with SuperCell
#'
#' SingleCorPlot Seurat function using metacell sizes
#'
#' @param object A Seurat metacell object
#' @inheritParams Seurat::SingleCorPlot
#' @return A SingleCorPlot scatter plot for metacells
#' @import Seurat
#' @export
#'
SingleCorPlot.SuperCell <- function (data, col.by = NULL, cols = NULL, pt.size = NULL,
                                      smooth = FALSE, rows.highlight = NULL, legend.title = NULL,
                                      na.value = "grey50", size.by, span = NULL, raster = NULL,
                                      raster.dpi = NULL, plot.cor = FALSE, jitter = TRUE)
{
  pt.size <- pt.size %||% AutoPointSize(data = data, raster = raster)
  if ((nrow(x = data) > 1e+05) & is.null(x = raster)) {
    message("Rasterizing points since number of points exceeds 100,000.",
            "\nTo disable this behavior set `raster=FALSE`")
  }
  raster <- raster %||% (nrow(x = data) > 1e+05)
  if (!is.null(x = raster.dpi)) {
    if (!is.numeric(x = raster.dpi) || length(x = raster.dpi) != 2)
      stop("'raster.dpi' must be a two-length numeric vector")
  }
  orig.names <- colnames(x = data)
  names.plot <- colnames(x = data) <- gsub(pattern = "-",
                                           replacement = ".", x = colnames(x = data), fixed = TRUE)
  names.plot <- colnames(x = data) <- gsub(pattern = ":",
                                           replacement = ".", x = colnames(x = data), fixed = TRUE)
  names.plot <- colnames(x = data) <- gsub(pattern = " ",
                                           replacement = ".", x = colnames(x = data), fixed = TRUE)
  if (ncol(x = data) < 2) {
    msg <- "Too few variables passed"
    if (ncol(x = data) == 1) {
      msg <- paste0(msg, ", only have ", colnames(x = data)[1])
    }
    stop(msg, call. = FALSE)
  }
  plot.cor <- if (isTRUE(x = plot.cor)) {
    # round(x = weights::wtd.cors(x = data[, 1], y = data[, 2],
    #                             weight = data[, size.by]), digits = 2)
    test <- weights::wtd.cor(x = data[, 1], y = data[, 2],
                             weight = data[, size.by])
    paste0("cor = ", round(test[,1],2) ,"\n p.val ",
           ifelse(test[,4] < 10^-16,
                  ifelse(test[,4] == 0, "= 0", "< e-16"),
                  paste0("= ", sprintf("%.2e", test[,4]))) )
  }else ("")
  if (!is.null(x = rows.highlight)) {
    highlight.info <- Seurat:::SetHighlight(cells.highlight = rows.highlight,
                                            cells.all = rownames(x = data), sizes.highlight = pt.size,
                                            cols.highlight = "red", col.base = "black", pt.size = pt.size,
                                            raster = raster)
    cols <- highlight.info$color
    col.by <- factor(x = highlight.info$highlight, levels = rev(x = highlight.info$plot.order))
    plot.order <- order(col.by)
    data <- data[plot.order, ]
    col.by <- col.by[plot.order]
  }
  if (!is.null(x = col.by)) {
    data$colors <- col.by
  }
  plot <- ggplot2::ggplot(data = data, mapping = ggplot2::aes_string(size = size.by,  x = names.plot[1], y = names.plot[2])) +
    ggplot2::labs(x = orig.names[1], y = orig.names[2], title = plot.cor, fill = legend.title)
  if (smooth) {
    plot <- plot + ggplot2::stat_density2d(mapping = aes(fill = ..density..^0.25),
                                           geom = "tile", contour = FALSE, n = 200,
                                           h = Bandwidth(data = data[, names.plot])) +
      ggplot2::scale_fill_continuous(low = "white", high = "dodgerblue4") + ggplot2::guides(fill = FALSE)
  }
  position <- NULL
  if (jitter) {
    position <- "jitter"
  }else {
    position <- "identity"
  }
  if (!is.null(x = col.by)) {
    if (raster) {
      plot <- plot + scattermore::geom_scattermore(mapping = ggplot2::aes_string(fill = "colors"),
                                                   position = position, pointsize = pt.size, pixels = raster.dpi)
    }else {
      if(length(unique(data[,size.by]))==1){
        plot <- plot + ggplot2::geom_point(mapping = ggplot2::aes_string(fill = "colors", color = "colors"),
                                           pch = 21, #color = "black",
                                           position = position, size = pt.size)
      }else{
        plot <- plot + ggplot2::geom_point(mapping = ggplot2::aes_string(fill = "colors",
                                                                         size = size.by), pch = 21, color = "black",
                                           position = position)
      }

    }
  }else {
    if (raster) {
      plot <- plot + geom_scattermore(position = position,
                                      pointsize = pt.size, pixels = raster.dpi)
    }else {
      plot <- plot + ggplot2::geom_point(mapping = ggplot2::aes_string(size = size.by),
                                         position = position, color = "black", pch = 21)
    }
  }
  if (!is.null(x = cols)) {
    cols.scale <- if (length(x = cols) == 1 && cols %in%
                      rownames(x = brewer.pal.info)) {
      ggplot2::scale_fill_brewer(palette = cols)
    }else {
      ggplot2::scale_fill_manual(values = cols, na.value = na.value)
    }

    cols.scale2 <- if (length(x = cols) == 1 && cols %in%
                       rownames(x = brewer.pal.info)) {
      ggplot2::scale_color_brewer(palette = cols)
    }else {
      ggplot2::scale_color_manual(values = cols, na.value = na.value)
    }
    plot <- plot + cols.scale + cols.scale2
    if (!is.null(x = rows.highlight)) {
      plot <- plot + ggplot2::guides(color = FALSE)
    }
    plot <- plot + ggplot2::guides(color = "none")
  }
  plot <- plot + cowplot::theme_cowplot() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  if (!is.null(x = span)) {
    plot <- plot + ggplot2::geom_smooth(mapping = ggplot2::aes_string(x = names.plot[1],
                                                                      y = names.plot[2]),
                                        method = "loess", span = span)
  }

  return(plot)
}


#' Dim plot with SuperCell
#'
#' DimPlot Seurat function using metacell sizes
#'
#' @param object A Seurat metacell object
#' @inheritParams Seurat::DimPlot
#' @return A FeatureScatter sctter plot for metacells
#' @import Seurat
#' @export

DimPlot.SuperCell <- function (object, dims = c(1, 2), cells = NULL, cols = NULL,
                               pt.size = NULL, reduction = NULL, group.by = NULL, split.by = NULL, size.by = "size",
                               shape.by = NULL, order = NULL, shuffle = FALSE, seed = 1,
                               label = FALSE, label.size = 4, label.color = "black", label.box = FALSE,
                               repel = FALSE, alpha = 1, cells.highlight = NULL, cols.highlight = "#DE2D26",
                               sizes.highlight = 1, na.value = "grey50", ncol = NULL, combine = TRUE,
                               raster = NULL, raster.dpi = c(512, 512))
{
  if (!rlang::is_integerish(x = dims, n = 2L, finite = TRUE) || !all(dims >
                                                                     0L)) {
    rlang::abort(message = "'dims' must be a two-length integer vector")
  }
  reduction <- reduction %||% DefaultDimReduc(object = object)
  cells <- cells %||% Seurat::Cells(x = object, assay = DefaultAssay(object = object[[reduction]]))
  dims <- paste0(Key(object = object[[reduction]]), dims)
  orig.groups <- group.by
  group.by <- group.by %||% "ident"
  data <- Seurat::FetchData(object = object, vars = c(dims, group.by,size.by),
                            cells = cells, clean = "project")

  #group.by <- colnames(x = data)[3:ncol(x = data)]
  for (group in group.by) {
    if (!is.factor(x = data[, group])) {
      data[, group] <- factor(x = data[, group])
    }
  }

  if(!is.null(x = size.by)) {
    data[,size.by] <- as.numeric(data[,size.by])
  }

  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  if (!is.null(x = split.by)) {
    split <- Seurat::FetchData(object = object, vars = split.by,
                               clean = TRUE)[split.by]
    data <- data[rownames(split), ]
    data[, split.by] <- split
  }
  if (isTRUE(x = shuffle)) {
    set.seed(seed = seed)
    data <- data[sample(x = 1:nrow(x = data)), ]
  }
  plots <- lapply(X = group.by, FUN = function(x) {
    plot <- SingleDimPlot.SuperCell(data = data[, c(dims, x, split.by,
                                                    shape.by,size.by)], dims = dims, col.by = x, cols = cols,
                                    pt.size = pt.size, shape.by = shape.by, order = order,
                                    alpha = alpha, label = FALSE, cells.highlight = cells.highlight,
                                    cols.highlight = cols.highlight, sizes.highlight = sizes.highlight,
                                    na.value = na.value, raster = raster, raster.dpi = raster.dpi)
    if (label) {
      plot <- Seurat::LabelClusters(plot = plot, id = x, repel = repel,
                                    size = label.size, split.by = split.by, box = label.box,
                                    color = label.color)
    }
    if (!is.null(x = split.by)) {
      plot <- plot + Seurat:::FacetTheme() + ggplot2::facet_wrap(facets = dplyr::vars(!!rlang::sym(x = split.by)),
                                                                 ncol = if (length(x = group.by) > 1 || is.null(x = ncol)) {
                                                                   length(x = unique(x = data[, split.by]))
                                                                 }
                                                                 else {
                                                                   ncol
                                                                 })
    }
    plot <- if (is.null(x = orig.groups)) {
      plot + ggplot2::labs(title = NULL)
    }
    else {
      plot + Seurat::CenterTitle()
    }
  })
  if (!is.null(x = split.by)) {
    ncol <- 1
  }
  if (combine) {
    plots <- patchwork::wrap_plots(plots, ncol = orig.groups %iff% ncol)
  }
  return(plots)
}



SingleDimPlot.SuperCell <- function (data, dims, col.by = NULL, cols = NULL, pt.size = NULL,size.by = "size",
                                     shape.by = NULL, alpha = 1, alpha.by = NULL, order = NULL,
                                     label = FALSE, repel = FALSE, label.size = 4, cells.highlight = NULL,
                                     cols.highlight = "#DE2D26", sizes.highlight = 1, na.value = "grey50",
                                     raster = NULL, raster.dpi = NULL)
{
  if ((nrow(x = data) > 1e+05) & is.null(x = raster)) {
    message("Rasterizing points since number of points exceeds 100,000.",
            "\nTo disable this behavior set `raster=FALSE`")
  }
  raster <- raster %||% (nrow(x = data) > 1e+05)
  pt.size <- pt.size %||% Seurat::AutoPointSize(data = data, raster = raster)
  if (!is.null(x = cells.highlight) && pt.size != Seurat:AutoPointSize(data = data,
                                                                       raster = raster) && sizes.highlight != pt.size && isTRUE(x = raster)) {
    warning("When `raster = TRUE` highlighted and non-highlighted cells must be the same size. Plot will use the value provided to 'sizes.highlight'.")
  }
  if (!is.null(x = raster.dpi)) {
    if (!is.numeric(x = raster.dpi) || length(x = raster.dpi) !=
        2)
      stop("'raster.dpi' must be a two-length numeric vector")
  }
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  if (!is.data.frame(x = data)) {
    data <- as.data.frame(x = data)
  }
  if (is.character(x = dims) && !all(dims %in% colnames(x = data))) {
    stop("Cannot find dimensions to plot in data")
  }
  else if (is.numeric(x = dims)) {
    dims <- colnames(x = data)[dims]
  }
  if (!is.null(x = cells.highlight)) {
    if (inherits(x = cells.highlight, what = "data.frame")) {
      stop("cells.highlight cannot be a dataframe. ", "Please supply a vector or list")
    }
    highlight.info <- Seurat:::SetHighlight(cells.highlight = cells.highlight,
                                            cells.all = rownames(x = data), sizes.highlight = sizes.highlight %||%
                                              pt.size, cols.highlight = cols.highlight, col.base = cols[1] %||%
                                              "#C3C3C3", pt.size = pt.size, raster = raster)
    order <- highlight.info$plot.order
    data$highlight <- highlight.info$highlight
    col.by <- "highlight"
    pt.size <- highlight.info$size
    cols <- highlight.info$color
  }
  if (!is.null(x = order) && !is.null(x = col.by)) {
    if (typeof(x = order) == "logical") {
      if (order) {
        data <- data[order(!is.na(x = data[, col.by]),
                           data[, col.by]), ]
      }
    }
    else {
      order <- rev(x = c(order, setdiff(x = unique(x = data[,
                                                            col.by]), y = order)))
      data[, col.by] <- factor(x = data[, col.by], levels = order)
      new.order <- order(x = data[, col.by])
      data <- data[new.order, ]
      if (length(x = pt.size) == length(x = new.order)) {
        pt.size <- pt.size[new.order]
      }
    }
  }
  if (!is.null(x = size.by) && !col.by %in% colnames(x = data)) {
    warning("Cannot find ", size.by, " in plotting data, not scaling point size")
    size.by <- NULL
  }

  if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
    warning("Cannot find ", col.by, " in plotting data, not coloring plot")
    col.by <- NULL
  }
  else {
    col.index <- match(x = col.by, table = colnames(x = data))
    if (grepl(pattern = "^\\d", x = col.by)) {
      col.by <- paste0("x", col.by)
    }
    else if (grepl(pattern = "-", x = col.by)) {
      col.by <- gsub(pattern = "-", replacement = ".",
                     x = col.by)
    }
    colnames(x = data)[col.index] <- col.by
  }
  if (!is.null(x = shape.by) && !shape.by %in% colnames(x = data)) {
    warning("Cannot find ", shape.by, " in plotting data, not shaping plot")
  }
  if (!is.null(x = alpha.by) && !alpha.by %in% colnames(x = data)) {
    warning("Cannot find alpha variable ", alpha.by, " in data, setting to NULL",
            call. = FALSE, immediate. = TRUE)
    alpha.by <- NULL
  }
  if(!is.null(x = size.by)) {
    data[,size.by] <- as.numeric(data[,size.by])
  }
  plot <- ggplot2::ggplot(data = data)
  plot <- if (isTRUE(x = raster)) {
    plot + scattermore::geom_scattermore(mapping = aes_string(x = dims[1],
                                                              y = dims[2], fill = paste0("`", col.by, "`"), shape = shape.by,
                                                              alpha = alpha.by), pointsize = pt.size, alpha = alpha,
                                         pixels = raster.dpi)
  }
  else {
    plot + ggplot2::geom_point(mapping =  ggplot2::aes_string(x = dims[1], y = dims[2], size = size.by,
                                                              fill = paste0("`", col.by, "`"),
                                                              shape = shape.by,
                                                              alpha = alpha.by),
                               color = "black",
                               pch = 21,
                               alpha = alpha)
  }
  plot <- plot + ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 3,
                                                                                   alpha = 1))) + ggplot2::labs(fill = NULL, title = col.by) + Seurat::CenterTitle()
  if (label && !is.null(x = col.by)) {
    plot <- Seurat::LabelClusters(plot = plot, id = col.by, repel = repel,
                                  size = label.size)
  }
  if (!is.null(x = cols)) {
    if (length(x = cols) == 1 && (is.numeric(x = cols) ||
                                  cols %in% rownames(x = brewer.pal.info))) {
      scale <- ggplot2::scale_fill_brewer(palette = cols, na.value = na.value)
    }
    else if (length(x = cols) == 1 && (cols %in% c("alphabet",
                                                   "alphabet2", "glasbey", "polychrome", "stepped"))) {
      colors <- Seurat::DiscretePalette(length(unique(data[[col.by]])),
                                        palette = cols)
      scale <-  ggplot2::scale_fill_manual(values = colors, na.value = na.value)
    }
    else {
      scale <-  ggplot2::scale_fill_manual(values = cols, na.value = na.value)
    }
    plot <- plot + scale
  }
  plot <- plot +  cowplot::theme_cowplot()
  return(plot)
}

#' Expand metacells to single cells
#'
#' Given size of metacells create a single cell seurat object
#' in which each metacell of size s is expanded (same data and metadata)
#' to s single cells.
#'
#'
#' @param object Seurat metacell object
#' @return A seurat single cell object of expanded metacells
#' @import Seurat
#' @export

ExpandMetacellSeuratAssay5 <- function(object,
                                       features,
                                       assay = "RNA",
                                       #slot = "data",
                                       meta.data.vars = NULL) {
  membership <- rep(1:ncol(object), object$size)
  DefaultAssay(object) <- assay
  feature.data <- GetAssayData(
    object = object, assay = assay, layer = "data"
  )

  meta.data <- FetchData(object,vars = c(c("orig.ident","ident"),meta.data.vars))
  expanded.meta.data <- meta.data[membership,]
  feature.data <- feature.data[features,]
  expanded.data <-  as(feature.data[,membership],"CsparseMatrix")
  colnames(expanded.data) <- rownames(expanded.meta.data)

  if(!sum(dim(object[[assay]]$counts)) == 0) {
    feature.counts <- GetAssayData(
      object = object, assay = assay, layer = "counts"
    )
    feature.counts <- feature.counts[features,]
    expanded.counts <-  as(feature.counts[,membership],"CsparseMatrix")
    colnames(expanded.counts) <- rownames(expanded.meta.data)
  } else {
    expanded.counts <- expanded.data
  }


  #colnames(expanded.data) <- colnames(expanded.data)
  expanded.sobj <- CreateSeuratObject(counts = expanded.counts,
                                      meta.data = expanded.meta.data,
                                      assay = assay)
  expanded.sobj[[assay]]$data <- expanded.data
  Idents(expanded.sobj) <- "ident"
  return(expanded.sobj)
}

#' Dot plot of metacells
#' Take into account metacell sizes using ExpandMetacellSeuratAssay5
#' Same parameter as the original DotPlot Seurat function
#'
#' @param object A Seurat metacell object
#' @inheritParams Seurat::DotPlot
#' @return As Seurat DotPLot, a ggplot object
#' @import Seurat
#' @export

DotPlot.SuperCell <- function (object, features, assay = NULL,  group.by = NULL, split.by = NULL,...)
{
  assay <- assay %||% DefaultAssay(object = object)
  features.plot <- features

  if (length(features) == 1) { # take another feature
    features <- c(c(features),rownames(object[[assay]])[rownames(object[[assay]])!= features][1])
  }

  object <- ExpandMetacellSeuratAssay5(object,features = features,
                                       assay = assay,
                                       meta.data.vars = c("ident",c(split.by,group.by)[!is.null(c(split.by,group.by))]))

  plot <- DotPlot(object, features.plot, assay = assay,  group.by = group.by, split.by = split.by,...)
  remove(object)
  gc()



  return(plot)
}

#' Violin plot of metacells
#' Take into account metacell sizes using ExpandMetacellSeuratAssay5
#' Same parameter as the original VlnPlot Seurat function
#'
#' @param object A Seurat metacell object
#' @inheritParams Seurat::VlnPlot
#' @return As Seurat VlnPLot, a patchworked ggplot object if combine = TRUE; otherwise, a list of ggplot objects
#' @import Seurat
#' @export

VlnPlot.SuperCell <- function(object, features, assay = NULL,  group.by = NULL, split.by = NULL,...) {
  assay <- assay %||% DefaultAssay(object = object)
  features.plot <- features
  if (length(features) == 1) { # take another feature
    features <- c(c(features),rownames(object[[assay]])[rownames(object[[assay]])!= features][1])
  }

  object <- ExpandMetacellSeuratAssay5(object,features = features,assay = assay,
                                       meta.data.vars = c("ident","orig.ident",c(split.by,group.by)[!is.null(c(split.by,group.by))]))

  plot <- VlnPlot(object, features.plot, assay = assay,  group.by = group.by, split.by = split.by,pt.size = 0,...)
  remove(object)
  gc()
  return(plot)
}
