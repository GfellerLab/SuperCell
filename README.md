# SuperCell2.0 enables semi-supervised construction of multimodal metacell atlases

SuperCell2.0 now handles single-cell multimodal data such as CITE-seq
(joint measurements of RNA and epitope in single cells) and 10X multiome
(joint measurements of RNA and ATAC in single nuclei).

Leveraging the Weighted Nearest Neighbor (WNN) framework of Seurat,
SuperCell2.0 performs **multimodal metacell identfication**.
SuperCell2.0 also proposes a **semi-supervised** workflow in which
partial cell annotation can be used to help metacell identification. See
our [tutorials](#tutorials) for examples.

## Coarse-graining of large single-cell data into metacells

SuperCell2.0 is an R package for coarse-graining large single-cell
(multi)omics data into metacells and performing downstream analysis at
the metacell level.

The exponential scaling of single cell data represents an important
hurdle for downstream analyses. One of the solutions to facilitate the
analysis of large-scale and noisy single-cell data is to merge highly
similar cells into *metacells*. This concept was first introduced by
[*Baran et al., 2019*](https://doi.org/10.1186/s13059-019-1812-2)
(MetaCell) and by [*Iacono et al.,
2018*](https://genome.cshlp.org/content/28/6/878) (bigSCale). More
recent methods to build *metacells* have been described in [*Ben-Kiki et
al. 2022*](https://doi.org/10.1186/s13059-022-02667-1) (MetaCell2),
[*Bilous et al., 2022*](https://doi.org/10.1186/s12859-022-04861-1)
(SuperCell) and [*Persad et al.,
2022*](https://doi.org/10.1038/s41587-023-01716-9) (SEACells). Despite
some differences in the implementation, all the methods are
network-based and can be summarized as follows:

**1.** A single-cell network is computed based on cell-to-cell
similarity.

**2.** Highly similar cells are identified as those forming dense
regions in the single-cell network and merged together into metacells
(coarse-graining)

**3.** information within each metacell is combined (average or sum).

**4.** Metacell data are used for the downstream analyses instead of
large-scale single-cell data

Unlike clustering, the aim of metacells is not to identify large groups
of cells that comprehensively capture biological concepts, like cell
types, but to merge cells that share highly similar profiles, and may
carry repetitive information. **Therefore metacells represent a
compromise structure that optimally remove redundant information in
single-cell data while preserving the biologically relevant
heterogeneity.**

An important concept when building metacells is the **graining level**
(*γ*), which we define as the ratio between the number of single cells
in the initial data and the number of metacells. We suggest applying *γ*
between 10 and 75, which significantly reduces the computational
resources needed to perform the downstream analyses while preserving
most of the result of the initial (i.e., single-cell) analyses.

## Installation

SuperCell2.0 requires [Seurat](https://github.com/satijalab/seurat) for
standard single-cell assays, such as RNA and Protein assays, and
[Signac](https://github.com/stuart-lab/signac) for chromatin assays as
weel as other R packages (the full list of dependencies is available
[here](DESCRIPTION)). 

To facilitate dependencies installation, we recommend to use the conda
[environment](tutorials/supercell_tuto_env.yaml) we provide for tutorials. Then
you can install everything using conda like this:

``` bash
conda env create -n supercell_tuto_env -f tutorials/supercell_tuto_env.yaml
conda activate supercell_tuto_env
Rscript -e "remotes::install_github('GfellerLab/SuperCell@develop',upgrade = 'never');library(SuperCell)"
```

Otherwise, you can install install SuperCell2.0 in R like this:

``` r
if (!requireNamespace("remotes")) install.packages("remotes")
remotes::install_github("GfellerLab/SuperCell")

library(SuperCell)
```

## Tutorials

1.  [Building and analyzing metacells in Bone Marrow CITE-seq data with
    SuperCell2.0](https://htmlpreview.github.io/?https://github.com/GfellerLab/SuperCell/blob/develop/docs/tutorials/SuperCell2.0_BM_CITE_seq.html)
2.  [Building and analyzing metacells in PBMC 10X multiome data with
    SuperCell2.0](https://htmlpreview.github.io/?https://github.com/GfellerLab/SuperCell/blob/develop/docs/tutorials/SuperCell2.0_PBMC_10x_multiome.html)

## [License]

SuperCell2.0 is developed by the group of David Gfeller at University of
Lausanne.

SuperCell2.0 is available under GPL-3 License.

For scientific questions, please contact Léonard Hérault
([leonard.herault\@gustaveroussy.fr](mailto:leonard.herault@gustaveroussy.fr))
or David Gfeller
([David.Gfeller\@unil.ch](mailto:David.Gfeller@unil.ch)).

## How to cite

If you use SuperCell2.0 in a publication, please cite:

-   [Bilous et al. Metacells untangle large and complex single-cell
    transcriptome networks, BMC Bioinformatics
    (2022).](https://doi.org/10.1186/s12859-022-04861-1)

-   [Bilous et al. Building and analyzing metacells in single-cell
    genomics data, Mol Syst Bio
    (2024).](https://doi.org/10.1038/s44320-024-00045-6)
