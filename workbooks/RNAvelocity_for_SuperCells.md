RNA velocity combining with super-cells
---------------------------------------

``` r
#if (!requireNamespace("remotes")) install.packages("remotes")
#remotes::install_github("GfellerLab/SuperCell")

library(SuperCell)
library(velocyto.R)
```

    ## Loading required package: Matrix

We show an example of how SuperCell can be used prior to velocyto.R to
compute RNA velocity.

### Load data and compute Super-cells’ spliced and un-spliced counts

We use a built-in *pancreas* dataset from [Bastidas-Ponce et al.
(2019)](https://journals.biologists.com/dev/article/146/12/dev173849/19483/Comprehensive-single-cell-mRNA-profiling-reveals-a)
For RNA velocity, we need spliced and un-spliced count matrices.

``` r
data("pancreas")
```

``` r
SC_pancreas <- SCimplify_for_velocity(pancreas$emat, 
                                        pancreas$nmat, gamma = 10)
```

### Compute RNA velocity

``` r
Vel <- supercell_estimate_velocity(SC_pancreas$emat, SC_pancreas$nmat)
```

    ## Warning in supercell_estimate_velocity(SC_pancreas$emat, SC_pancreas$nmat): supercell_size was replaced with 1

    ## calculating cell knn ... done
    ## calculating convolved matrices ... done
    ## fitting gamma coefficients ... done. succesfful fit for 14807 genes
    ## filtered out 3406 out of 14807 genes due to low nmat-emat correlation
    ## filtered out 1043 out of 11401 genes due to low nmat-emat slope
    ## calculating RNA velocity shift ... done
    ## calculating extrapolated cell state ... done

### Plot RNA velocity on a tSNE’s coordinates

``` r
# Assign clusters to super-cells based on the given single clustering
clusters <- supercell_assign(pancreas$meta$clusters, SC_pancreas$membership)
```

#### Set up the color scheme

``` r
# Set up the color scheme
N.clusters <- length(unique(clusters))
pal <- setNames(colorRampPalette(RColorBrewer::brewer.pal(8, name= "Set1"))(N.clusters), 
                as.character(unique(clusters))) 
color <- setNames(pal[as.character(clusters)], names(clusters))

map_cluster_to_cluster_name <- c('Ngn3 low EP', 'Alpha', 'Delta',
                                          'Beta', 'Pre-endocrine', 'Ngn3 high EP',
                                          'Ductal', 'Epsilon')
```

#### Plot tSNE

``` r
size <- sqrt(1+table(SC_pancreas$membership))/sqrt(2)/2
set.seed(3)
tsne10 <- tSNE.velocity.plot(Vel, 
                             cell.colors = ac(color, alpha = 0.6), 
                             scale = "log", do.par = T,
                             delta.norm = FALSE, nPcs = 15, norm.nPcs = 15 * 10, perplexity = 50,
                             show.grid.flow = FALSE, grid.n = 20, grid.sd = NULL, min.grid.cell.mass = 1,
                             pcount = 1, verbose = TRUE, min.arrow.median.ratio = 1/10,
                             max.arrow.quantile = 0.9, arrow.scale = 1, arrow.lwd = 1,
                             xlab = "tSNE-x", ylab = "tSNE-y", size.norm = FALSE, cex = size
)
```

    ## log ... pca ... tSNE ...Performing PCA
    ## Read the 740 x 15 data matrix successfully!
    ## Using no_dims = 2, perplexity = 50.000000, and theta = 0.500000
    ## Computing input similarities...
    ## Building tree...
    ## Done in 0.10 seconds (sparsity = 0.249310)!
    ## Learning embedding...
    ## Iteration 50: error is 60.800991 (50 iterations in 0.11 seconds)
    ## Iteration 100: error is 50.997427 (50 iterations in 0.10 seconds)
    ## Iteration 150: error is 50.493617 (50 iterations in 0.10 seconds)
    ## Iteration 200: error is 50.302629 (50 iterations in 0.11 seconds)
    ## Iteration 250: error is 50.263116 (50 iterations in 0.10 seconds)
    ## Iteration 300: error is 0.470060 (50 iterations in 0.10 seconds)
    ## Iteration 350: error is 0.376875 (50 iterations in 0.11 seconds)
    ## Iteration 400: error is 0.350275 (50 iterations in 0.11 seconds)
    ## Iteration 450: error is 0.342156 (50 iterations in 0.11 seconds)
    ## Iteration 500: error is 0.337179 (50 iterations in 0.11 seconds)
    ## Iteration 550: error is 0.331093 (50 iterations in 0.11 seconds)
    ## Iteration 600: error is 0.325433 (50 iterations in 0.11 seconds)
    ## Iteration 650: error is 0.323196 (50 iterations in 0.10 seconds)
    ## Iteration 700: error is 0.320226 (50 iterations in 0.11 seconds)
    ## Iteration 750: error is 0.318516 (50 iterations in 0.11 seconds)
    ## Iteration 800: error is 0.316694 (50 iterations in 0.11 seconds)
    ## Iteration 850: error is 0.315215 (50 iterations in 0.11 seconds)
    ## Iteration 900: error is 0.314090 (50 iterations in 0.11 seconds)
    ## Iteration 950: error is 0.313261 (50 iterations in 0.11 seconds)
    ## Iteration 1000: error is 0.312096 (50 iterations in 0.11 seconds)
    ## Fitting performed in 2.13 seconds.
    ## delta norm ... done

``` r
legend("topleft", legend=map_cluster_to_cluster_name,
       fill=pal,  cex=0.8)
```

![](RNAvelocity_for_SuperCells_files/figure-markdown_github/Plot%20tSNE-1.png)
