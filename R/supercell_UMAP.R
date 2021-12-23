#' Compute UMAP of super-cells
#'
#' Computes UMAP of super-cells
#' @param SC super-cell structure (output of \link{SCimplify}) with a field \code{PCA_name} containig PCA result
#' @param PCA_name name of \code{SC} field containing result of \link{supercell_prcomp}
#' @param n.comp number of vector of principal components to use for computing UMAP
#' @param ... other parameters of \link[umap]{umap}
#'
#' @return \link[umap]{umap} result
#'
#' @export
supercell_UMAP <- function(
  SC,
  PCA_name = "SC_PCA",
  n.comp = NULL,
  n_neighbors = 15,
  ...
){

  if(!(PCA_name %in% names(SC))){
    stop(paste("PCA_name", PCA_name, "field is not found in SC. Please, first compute pca with `SC[['SC_PCA']] <- supercell_prcomp(SC)` or provide correct name of PCA field!"))
  }

  SC_pca <- SC[[PCA_name]]

  if(is.null(n.comp)){
    n.comp <- 1:ncol(SC_pca$x)
  } else if(length(n.comp) == 1){
    n.comp <- 1:n.comp
  }

  SC_pca <- SC_pca$x[,n.comp]

  custom.config <- umap::umap.defaults
  custom.config$n_neighbors <- n_neighbors

  SC_umap <- umap::umap(SC_pca, config = custom.config, ...)

  return(SC_umap)
}

#' Plot super-cell UMAP
#' Plots super-cell UMAP (result of \link{supercell_UMAP})
#'
#' @param
#'
#' @return \link[ggplot2]{ggplot}
#' @export

supercell_plot_UMAP <- function(
  SC,
  groups,
  UMAP_name = "SC_UMAP",
  color.use = NULL,
  asp = 1,
  alpha = 0.7,
  title = NULL,
  ...
){

  N.SC <- SC$N.SC

  if(!(UMAP_name %in% names(SC))){
    stop(paste("UMAP_name", UMAP_name, "field is not found in SC. Please, first compute umap with `SC[['SC_UMAP']] <- supercell_UMAP(SC)` or provide correct name of UMAP field!"))
  }

  if(length(groups) == 1){ #groups stands for SC field name
    if(!(groups %in% names(SC))){
      stop(paste("groups", groups, "is not found in SC. Please provide correct  group name field"))
    }
    real_groups <- SC[[groups]]
  } else {
    real_groups <- groups
  }

  if(length(real_groups) != N.SC){
    stop(paste0("Length of groups (n = ", length(real_groups), ") != real number of super-cells (n = ", N.SC, ")"))
  }

  lay.df <- data.frame(SC[[UMAP_name]]$layout)
  lay.df$groups <- real_groups
  lay.df$supercell_size <- SC$supercell_size


  g <- ggplot2::ggplot(
    lay.df,
    ggplot2::aes(x = X1, y = X2, color = groups, fill = groups, size = sqrt(supercell_size))
  ) +
    ggplot2::scale_size_continuous(range = c(0.5, 0.5*max(log1p((SC$supercell_size))))) +
    ggplot2::labs(x = "UMAP-1", y = "UMAP-2",  title = title) +
    ggplot2::geom_point(alpha = alpha) +
    ggplot2::theme_classic() +
    ggplot2::theme(asp = asp)

  if(!is.null(color.use)){
    g <- g + ggplot2::scale_color_manual(values = color.use)
    g <- g + ggplot2::scale_fill_manual(values = color.use)
  }

 plot(g)

 return(invisible(g))
}

