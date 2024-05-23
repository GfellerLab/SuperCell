## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "figures/",
  fig.width = 6, fig.height = 6,
  eval = FALSE
)

## ----library, warning=FALSE, eval=FALSE---------------------------------------
#  if (!requireNamespace("remotes")) install.packages("remotes")
#  remotes::install_github("GfellerLab/SuperCell")

## ----load libraries, warning=FALSE--------------------------------------------
#  library(SuperCell)
#  library(velocyto.R)

## ----load data----------------------------------------------------------------
#  data("pancreas")

## ----compute metacells for RNA velocity---------------------------------------
#  SC_pancreas <- SCimplify_for_velocity(
#    pancreas$emat,
#    pancreas$nmat,
#    gamma = 10)

## ----compute RNA velocity-----------------------------------------------------
#  Vel <- supercell_estimate_velocity(SC_pancreas$emat, SC_pancreas$nmat)

## ----plot RNA velocity--------------------------------------------------------
#  # Assign clusters to metacells based on the given single clustering
#  clusters <- supercell_assign(pancreas$meta$clusters, SC_pancreas$membership)

## ----set up the color scheme--------------------------------------------------
#  # Set up the color scheme
#  N.clusters <- length(unique(clusters))
#  pal <- setNames(colorRampPalette(RColorBrewer::brewer.pal(8, name= "Set1"))(N.clusters),
#                  as.character(unique(clusters)))
#  color <- setNames(pal[as.character(clusters)], names(clusters))
#  
#  map_cluster_to_cluster_name <- c('Ngn3 low EP', 'Alpha', 'Delta', 'Beta', 'Pre-endocrine', 'Ngn3 high EP', 'Ductal', 'Epsilon')

## ----Plot tSNE, fig.asp=1, fig.height=5, fig.width=5, fig.ext='pdf'-----------
#  size <- sqrt(1+table(SC_pancreas$membership))/sqrt(2)/2
#  set.seed(3)
#  tsne10 <- tSNE.velocity.plot(Vel,
#                               cell.colors = ac(color, alpha = 0.6),
#                               scale = "log", do.par = T,
#                               delta.norm = FALSE, nPcs = 15, norm.nPcs = 15 * 10, perplexity = 50,
#                               show.grid.flow = FALSE, grid.n = 20, grid.sd = NULL, min.grid.cell.mass = 1,
#                               pcount = 1, verbose = TRUE, min.arrow.median.ratio = 1/10,
#                               max.arrow.quantile = 0.9, arrow.scale = 1, arrow.lwd = 1,
#                               xlab = "tSNE-x", ylab = "tSNE-y", size.norm = FALSE, cex = size
#  )
#  legend("topleft", legend=map_cluster_to_cluster_name,
#         fill=pal,  cex=0.8)
#  

