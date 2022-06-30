## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "figures/",
  fig.width = 6, fig.height = 6
)

## ----package and data loading-------------------------------------------------
library(SuperCell)
data(cell_lines)

GE <- cell_lines$GE
cell.meta <- cell_lines$meta


## ----parameters---------------------------------------------------------------
gamma <- 20 # graining level
n.pc  <- 10 # number of PCs

## ----two samples--------------------------------------------------------------
cell.idx.HCC827 <- which(cell.meta == "HCC827")
cell.idx.H838   <- which(cell.meta == "H838")

## ----combined scimplify-------------------------------------------------------
SC.HCC827.H838 <- SCimplify(
  GE[,c(cell.idx.HCC827, cell.idx.H838)],  # log-normalized gene expression matrix
  gamma = gamma, # graining level
  cell.split.condition = cell.meta[c(cell.idx.HCC827, cell.idx.H838)], # metacell do not mix cells from different cell lines
  n.pc = n.pc) # number of proncipal components to use

genes.use <- SC.HCC827.H838$genes.use

SC.HCC827.H838$cell.line <- supercell_assign(cell.meta[c(cell.idx.HCC827, cell.idx.H838)], 
                                             supercell_membership = SC.HCC827.H838$membership)

SC.GE.HCC827.H838 <- supercell_GE(GE[,c(cell.idx.HCC827, cell.idx.H838)], groups = SC.HCC827.H838$membership)

SC.HCC827.H838$SC_PCA <- supercell_prcomp(
  Matrix::t(SC.GE.HCC827.H838),
  supercell_size = SC.HCC827.H838$supercell_size, 
  genes.use = genes.use)

SC.HCC827.H838$SC_UMAP <- supercell_UMAP(
  SC.HCC827.H838, 
  n_neighbors = 10)

supercell_plot_UMAP(
    SC.HCC827.H838,
    group = "cell.line",
    title = paste0("Combined construction of HCC827 and H838 metacells")
  )


## ----independent scimplify----------------------------------------------------
SC.HCC827 <- SCimplify(GE[,cell.idx.HCC827],  # log-normalized gene expression matrix
                gamma = gamma, # graining level
                n.pc = n.pc, # number of proncipal components to use
                genes.use = genes.use) # using the same set of genes as for the combined analysis

SC.HCC827$cell.line <- supercell_assign(cell.meta[cell.idx.HCC827], supercell_membership = SC.HCC827$membership)

SC.H838 <- SCimplify(GE[,cell.idx.H838],  # log-normalized gene expression matrix
                gamma = gamma, # graining level
                n.pc = n.pc, # number of proncipal components to use
                genes.use = genes.use) # using the same set of genes as for the combined analysis
SC.H838$cell.line <- supercell_assign(cell.meta[cell.idx.H838], supercell_membership = SC.H838$membership)

SC.merged <- supercell_merge(list(SC.HCC827, SC.H838), fields = c("cell.line"))

# compute metacell gene expression for SC.HCC827
SC.GE.HCC827 <- supercell_GE(GE[, cell.idx.HCC827], groups = SC.HCC827$membership)
# compute metacell gene expression for SC.H838
SC.GE.H838 <- supercell_GE(GE[, cell.idx.H838], groups = SC.H838$membership)
# merge GE matricies
SC.GE.merged <- supercell_mergeGE(list(SC.GE.HCC827, SC.GE.H838))

SC.merged$SC_PCA <- supercell_prcomp(
  Matrix::t(SC.GE.merged),
  supercell_size = SC.merged$supercell_size, 
  genes.use = genes.use)

SC.merged$SC_UMAP <- supercell_UMAP(
  SC.merged, 
  n_neighbors = 10)

g <- supercell_plot_UMAP(
    SC.merged,
    group = "cell.line",
    title = paste0("Independent construction of HCC827 and H838 metacells")
  )


## ----heatmap metacell membership----------------------------------------------
heatmap(as.matrix(table(SC.merged$membership, SC.HCC827.H838$membership)), scale = "none")

## ----size distribution--------------------------------------------------------
summary(SC.merged$supercell_size)
summary(SC.HCC827.H838$supercell_size)

## ----combined vs independent analysis-----------------------------------------
## Combined analysis
# actual graining level for H838 cell line
length(cell.idx.H838)/sum(SC.HCC827.H838$cell.line == "H838")
# actual graining level for H838 cell line
length(cell.idx.HCC827)/sum(SC.HCC827.H838$cell.line == "HCC827")

## Independent analysis 
# actual graining level for H838 cell line
length(cell.idx.H838)/sum(SC.merged$cell.line == "H838")
# actual graining level for HCC827 cell line
length(cell.idx.HCC827)/sum(SC.merged$cell.line == "HCC827")
# actual overall graining level in the combined analysis
length(SC.merged$membership)/max(SC.merged$membership)


