[![R-CMD-check](https://github.com/GfellerLab/SuperCell/workflows/R-CMD-check/badge.svg)](https://github.com/GfellerLab/SuperCell/actions)
[1E![DOI](https://img.shields.io/badge/DOI%3A-0.1186/s12859--022--04861--1-brightgreen)](https://doi.org/10.1186/s12859-022-04861-1)

# Coarse-graining of large single-cell RNA-seq data into metacells

SuperCell is an R package for coarse-graining large single-cell RNA-seq
data into metacells and performing downstream analysis at the metacell
level.

The exponential scaling of scRNA-seq data represents an important hurdle
for downstream analyses. One of the solutions to facilitate the analysis
of large-scale and noisy scRNA-seq data is to merge transcriptionally
highly similar cells into *metacells*. This concept was first introduced
by [*Baran et al., 2019*](https://doi.org/10.1186/s13059-019-1812-2)
(MetaCell) and by [*Iacono et al.,
2018*](https://genome.cshlp.org/content/28/6/878) (bigSCale). More recent
methods to build *metacells* have been described in [*Ben-Kiki et
al. 2022*](https://doi.org/10.1186/s13059-022-02667-1) (MetaCell2),
[*Bilous et al., 2022*](https://doi.org/10.1186/s12859-022-04861-1)
(SuperCell) and [*Persad et al.,
2022*](https://doi.org/10.1038/s41587-023-01716-9) (SEACells). Despite
some differences in the implementation, all the methods are
network-based and can be summarized as follows:

**1.** A single-cell network is computed based on cell-to-cell
similarity (in transcriptomic space)

**2.** Highly similar cells are identified as those forming dense
regions in the single-cell network and merged together into metacells
(coarse-graining)

**3.** Transcriptomic information within each metacell is combined
(average or sum).

**4.** Metacell data are used for the downstream analyses instead of
large-scale single-cell data


Unlike clustering, the aim of metacells is not to identify large groups
of cells that comprehensively capture biological concepts, like cell
types, but to merge cells that share highly similar profiles, and may
carry repetitive information. **Therefore metacells represent a
compromise structure that optimally remove redundant information in
scRNA-seq data while preserving the biologically relevant
heterogeneity.**

An important concept when building metacells is the **graining level**
(*γ*), which we define as the ratio between the number of single cells
in the initial data and the number of metacells. We suggest applying *γ*
between 10 and 50, which significantly reduces the computational
resources needed to perform the downstream analyses while preserving
most of the result of the initial (i.e., single-cell) analyses.

## Installation

SuperCell requires
[igraph](https://CRAN.R-project.org/package=igraph),
[RANN](https://CRAN.R-project.org/package=RANN),
[WeightedCluster](https://CRAN.R-project.org/package=WeightedCluster),
[corpcor](https://CRAN.R-project.org/package=corpcor),
[weights](https://CRAN.R-project.org/package=weights),
[Hmisc](https://CRAN.R-project.org/package=Hmisc),
[Matrix](https://CRAN.R-project.org/package=Matrix),
[matrixStats](https://CRAN.R-project.org/package=matrixStats),
[plyr](https://CRAN.R-project.org/package=plyr),
[irlba](https://CRAN.R-project.org/package=irlba),
[grDevices](https://stat.ethz.ch/R-manual/R-devel/library/grDevices/html/00Index.html),
[patchwork](https://CRAN.R-project.org/package=patchwork),
[ggplot2](https://CRAN.R-project.org/package=ggplot2).
SuperCell uses [velocyto.R](https://github.com/velocyto-team/velocyto.R)
for RNA velocity.

``` r
install.packages("igraph")
install.packages("RANN")
install.packages("WeightedCluster")
install.packages("corpcor")
install.packages("weights")
install.packages("Hmisc")
install.packages("Matrix")
install.packages("patchwork")
install.packages("plyr")
install.packages("irlba")
```

Installing SuperCell package from gitHub

``` r
if (!requireNamespace("remotes")) install.packages("remotes")
remotes::install_github("GfellerLab/SuperCell")

library(SuperCell)
```

## Examples

1.  [Building and analyzing metacells with
    SuperCell](./vignettes/a_SuperCell.Rmd)
2.  [Building metacells with SuperCell and alayzing them with a standard
    Seurat
    pipeline](https://github.com/GfellerLab/SIB_workshop/blob/main/workbooks/Workbook_1__cancer_cell_lines.md)
3.  [Data integration of metacells built with
    SuperCell](https://github.com/GfellerLab/SIB_workshop/blob/main/workbooks/Workbook_2__COVID19_integration.md)

## [License]

SuperCell is developed by the group of David Gfeller at University of
Lausanne.

SuperCell is available under GPL-3 License.

For scientific questions, please contact Mariia Bilous
(<mariia.bilous@unil.ch>) or David Gfeller (<David.Gfeller@unil.ch>).

## How to cite

If you use SuperCell in a publication, please cite: [Bilous et
al. Metacells untangle large and complex single-cell transcriptome
networks, BMC Bioinformatics
(2022).](https://doi.org/10.1186/s12859-022-04861-1)
