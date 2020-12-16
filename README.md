Installation
============

``` r
#devtools::install_github("GfellerLab/SuperCell")
library(SuperCell)
```

Analysis
========

Load data of 5 cancer cell lines sequenced with 10x from [Tian et al., 2019](https://doi.org/10.1038/s41592-019-0425-8)
-----------------------------------------------------------------------------------------------------------------------

``` r
GE <- readRDS("./data/5cancer_cell_lines_10x_GE.Rds") # load gene expression matrix 
dim(GE) # genes as rows and cells as columns
## [1] 11786  3918
cell.meta <- readRDS("./data/5cancer_cell_lines_10x_cell_line_assignment.Rds") # load cell assignment to a cancer cell line
```

Simplify single-cell data reducing numer of super-cell *g**a**m**m**a* = 20 times comparing to the number of single cells
-------------------------------------------------------------------------------------------------------------------------

``` r
gamma <- 20
SC    <- SCimplify(GE, k.knn = 5, gamma = gamma, n.var.genes = 1000, use.nn2 = FALSE)
set.seed(12345)
p.SC <- supercell_plot(SC$graph.supercells, main = paste("Super-cell network, gamma =", gamma), seed = 1)
```

![](figures/Simplification-1.png)

``` r
p.sc <- supercell_plot(SC$graph.singlecell, group = cell.meta, do.frames = F, main = paste("Single-cell network, N =", dim(GE)[2]), lay.method = "components")
```

![](figures/Simplification-2.png)

Compute gene expression for simplified data
-------------------------------------------

``` r
SC.GE <- supercell_GE(GE, SC$membership)
dim(SC.GE)
## [1] 11786   196
```

\#Map each super-cell to a particular cell line

``` r
SC2cellline  <- supercell_assign(cell.meta, SC$membership)
SC$cell_line <- SC2cellline

p <- supercell_plot(SC$graph.supercells, group = SC$cell_line, seed = 1, main = "Super-cell colored by cell line assignment")
```

![](figures/unnamed-chunk-5-1.png)

Cluster super-cell data

``` r
SC.PCA         <- supercell_prcomp(t(SC.GE), genes.use = SC$genes.use, supercell_size = SC$supercell_size, k = 20)
D              <- dist(SC.PCA$x)

SC.clusters    <- supercell_cluster(D = D, k = 5, supercell_size = SC$supercell_size)
SC$clustering  <- SC.clusters$clustering
```

``` r
map.cluster.to.cell.line    <- supercell_assign(supercell_membership = SC$clustering, clusters  = SC$cell_line)
SC$clustering_reordered     <- map.cluster.to.cell.line[SC$clustering]
p <- supercell_plot(SC$graph.supercells, group = SC$clustering_reordered, seed = 1, main = "Super-cell colored by cell line cluster")
```

![](figures/unnamed-chunk-7-1.png)

Differential expression analysis of clustered super-cell data

``` r
markers.all.positive <- supercell_FindAllMarkers(ge = SC.GE, 
                                                 supercell_size = SC$supercell_size,
                                                 clusters = SC$clustering_reordered,
                                                 logfc.threshold = 1,
                                                 only.pos = T)
markers.all.positive$H1975[1:20,]
##         p.value adj.p.value     pct.1     pct.2    logFC  w.mean.1   w.mean.2
## DHRS2         0           0 1.0000000 0.8224648 4.079618 3.6237875 0.10054192
## MT1E          0           0 1.0000000 0.9339270 3.702128 4.6560576 0.69254917
## PEG10         0           0 1.0000000 0.9388107 2.561314 3.2359195 0.89548751
## LGALS1        0           0 1.0000000 1.0000000 2.394566 6.5798440 4.00039835
## S100A2        0           0 1.0000000 0.9956909 2.238731 3.1041134 1.24595811
## ZNF880        0           0 1.0000000 0.6058604 1.819325 1.6675772 0.08666412
## MT2A          0           0 1.0000000 1.0000000 1.718104 6.3843539 4.12274024
## CT45A2        0           0 1.0000000 0.5808676 1.691487 1.7640705 0.27318252
## XAGE1B        0           0 1.0000000 1.0000000 1.667468 6.5436587 3.66859782
## HSPB1         0           0 1.0000000 1.0000000 1.645820 5.5983707 3.87890578
## PRDX2         0           0 1.0000000 0.9686872 1.631436 3.4980150 1.33455473
## IFI27         0           0 1.0000000 0.9488653 1.620494 3.5275622 1.32877161
## MT1P1         0           0 1.0000000 0.9951163 1.548610 3.2663495 1.50666087
## TMEM134       0           0 1.0000000 1.0000000 1.492498 2.4398008 1.01156248
## TNNT1         0           0 1.0000000 0.9824763 1.490033 3.2983849 1.63217280
## DMKN          0           0 1.0000000 0.9859236 1.459521 2.6011521 0.83344367
## XAGE1A        0           0 1.0000000 0.9801781 1.458303 2.8883699 1.06459003
## MT2P1         0           0 1.0000000 0.9948291 1.451145 2.7769310 1.18613463
## CAV1          0           0 1.0000000 1.0000000 1.390797 4.1486580 2.37528254
## MGP           0           0 0.9244851 0.2074117 1.389683 0.6735631 0.01058892
```
