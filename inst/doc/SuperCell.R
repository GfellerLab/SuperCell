## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "figures/",
  fig.width = 6, fig.height = 6
)

## ----library, warning=FALSE---------------------------------------------------
if (!requireNamespace("remotes")) install.packages("remotes") 
remotes::install_github("GfellerLab/SuperCell")

library(SuperCell)

## ----load data----------------------------------------------------------------
data(cell_lines) # list with GE - gene expression matrix (logcounts), meta - cell meta data
GE <- cell_lines$GE
dim(GE) # genes as rows and cells as columns
cell.meta <- cell_lines$meta

## ----Simplification, warning=FALSE, paged.print=FALSE-------------------------
gamma <- 20 # graining level
k.knn <- 5

SC <- SCimplify(GE,  # gene expression matrix 
                k.knn = k.knn, # number of nearest neighbors to build kNN network
                gamma = gamma, # graining level
                n.var.genes = 1000 # number of the top variable genes to use for dimentionality reduction 
)


# plot network of metacells

supercell_plot(SC$graph.supercells, # network
               color.use = "gray", # color of the nodes
               main = paste("Metacell network, gamma =", gamma), 
               seed = 1) 

# plot single-cell network
supercell_plot(SC$graph.singlecell, # network
               group = cell.meta, # colored by cell line assignment
               do.frames = F, # not drawing frames around each node 
               main = paste("Single-cell network, N =", dim(GE)[2]), 
               lay.method = "components") # method to compute the network 2D embedding 

## ----average gene expression--------------------------------------------------
SC.GE <- supercell_GE(GE, SC$membership)
dim(SC.GE) 

## ----assign metacells to cell line infromation--------------------------------
SC$cell_line <- supercell_assign(clusters = cell.meta, # single-cell assigment to cell lines (clusters)
                                 supercell_membership = SC$membership, # single-cell assignment to metacells
                                 method = "jaccard")


seed <- 1 # seed for network plotting 

# plot network of metacells colored by cell line assignment 
supercell_plot(SC$graph.supercells, 
               group = SC$cell_line, 
               seed = seed, 
               main = "Metacells colored by cell line assignment")

## ----purity of supercell in terms of cell line composition--------------------
# compute purity of metacells in terms of cell line composition
purity <- supercell_purity(clusters = cell.meta, 
                           supercell_membership = SC$membership, method = 'entropy')
hist(purity, main = "Purity of metacells \nin terms of cell line composition")

## ----plotting options---------------------------------------------------------
## rotate network to be more consistent with the single-cell one
supercell_plot(SC$graph.supercells, 
               group = SC$cell_line, 
               seed = seed, 
               alpha = -pi/2,
               main  = "Metacells colored by cell line assignment (rotated)")

## alternatively, any layout can be provided as 2xN numerical matrix, where N is number of nodes (cells)

## Let's plot metacell network using the layout of the single-cell network:
## 1) get single-cell network layout 
my.lay.sc <- igraph::layout_components(SC$graph.singlecell) 

## 2) compute metacell network layout averaging coordinates withing metacells
my.lay.SC <- Matrix::t(supercell_GE(ge = t(my.lay.sc), groups = SC$membership))

## 3) provide layout with the parameter $lay$
supercell_plot(SC$graph.supercells, 
               group = SC$cell_line, 
               lay = my.lay.SC,
               main  = "Metacells colored by cell line assignment (averaged coordinates)")

## ----clustering---------------------------------------------------------------
#dimensionality reduction 
SC.PCA         <- supercell_prcomp(Matrix::t(SC.GE), # metacell gene exptression matrix
                                   genes.use = SC$genes.use, # genes used for the coarse-graining, but any set can be provided
                                   supercell_size = SC$supercell_size, # sample-weighted pca
                                   k = 20) 
## compute distance
D              <- dist(SC.PCA$x)

## cluster metacells
SC.clusters    <- supercell_cluster(D = D, k = 5, supercell_size = SC$supercell_size) 
SC$clustering  <- SC.clusters$clustering

## ----assign metacell clustering results to cell line information--------------
## mapping metacell cluster to cell line 
map.cluster.to.cell.line    <- supercell_assign(supercell_membership = SC$clustering, clusters  = SC$cell_line)
## clustering as cell line
SC$clustering_reordered     <- map.cluster.to.cell.line[SC$clustering]

supercell_plot(SC$graph.supercells, 
               group = SC$clustering_reordered, 
               seed = seed,
               alpha = -pi/2,
               main = "Metacells colored by cluster")

## ----differential expression analysis-----------------------------------------
markers.all.positive <- supercell_FindAllMarkers(ge = SC.GE, # metacell gene expression matrix
                                                 supercell_size = SC$supercell_size, # size of metacell for sample-weighted method
                                                 clusters = SC$clustering_reordered, # clustering
                                                 logfc.threshold = 1, # mininum log fold-change
                                                 only.pos = T) # keep only upregulated genes
markers.all.positive$H1975[1:20,]

## ----Violin plots-------------------------------------------------------------
genes.to.plot <- c("DHRS2", "MT1P1", "TFF1", "G6PD", "CCL2", "C1S")

supercell_VlnPlot(ge = SC.GE, 
                  supercell_size = SC$supercell_size, 
                  clusters = SC$clustering_reordered,
                  features = genes.to.plot,
                  idents = c("H1975", "H2228", "A549"), 
                  ncol = 3)

supercell_GeneGenePlot(ge = SC.GE, 
                       gene_x = genes.to.plot[1:3],
                       gene_y = genes.to.plot[4:6],
                       supercell_size = SC$supercell_size, 
                       clusters = SC$clustering_reordered,)

## -----------------------------------------------------------------------------
SC10 <- supercell_rescale(SC, gamma = 10)

SC10$cell_line <- supercell_assign(clusters = cell.meta, # single-cell assigment to cell lines (clusters)
                                 supercell_membership = SC10$membership, # single-cell assignment to metacells
                                 method = "jaccard")

supercell_plot(SC10$graph.supercells, 
               group = SC10$cell_line, 
               seed = 1,
               main  = "Metacells at gamma = 10 colored by cell line assignment")

### don't forget to recompute metacell gene expression matrix for a new grainig level with 
# GE10 <- supercell_GE(GE, SC10$membership)
### if you are going to perform downstream analyses at the new graining level

## ----Seurat-------------------------------------------------------------------
#install.packages("Seurat")
library(Seurat)

m.seurat <- supercell_2_Seurat(SC.GE = SC.GE, SC = SC, fields = c("cell_line", "clustering", "clustering_reordered"))

## ----PCAplot------------------------------------------------------------------
PCAPlot(m.seurat)

### cluster SuperCell network (unweighted clustering)
m.seurat <- FindClusters(m.seurat, graph.name = "RNA_nn") # now RNA_nn is metacell network

m.seurat <- FindNeighbors(m.seurat, verbose = FALSE)  # RNA_nn has been replaced with kNN graph of metacell (unweigted)
m.seurat <- FindClusters(m.seurat, graph.name = "RNA_nn") 

