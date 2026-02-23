# SuperCell2.0 enables semi-supervised construction of multimodal metacell atlases

SuperCell2.0 now handles single-cell multimodal data such as CITE-seq (joint measurements of RNA and epitope in single cells) and 10X multiome (joint measurements of RNA and ATAC in single nuclei).

Leveraging the Weighted Nearest Neighbor (WNN) framework of [Seurat](https://satijalab.org/seurat/), SuperCell2.0 performs **multimodal metacell identfication**. SuperCell2.0 also proposes a **semi-supervised** workflow in which partial cell annotation can be used to help metacell identification. See our [tutorials](#tutorials) for examples.

<p align="center">

<img src="docs/figure_1_workflow.png" width="750"/>

</p>

## Installation

SuperCell2.0 requires [Seurat](https://github.com/satijalab/seurat) for standard single-cell assays, such as RNA and Protein assays, and [Signac](https://github.com/stuart-lab/signac) for chromatin assays as weel as other R packages (the full list of dependencies is available [here](DESCRIPTION)).

To facilitate dependencies installation, we recommend to use the conda [environment](tutorials/supercell_tuto_env.yaml) we provide for tutorials. Then you can install everything using conda like this:

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

1.  [Building and analyzing metacells in Bone Marrow CITE-seq data with SuperCell2.0](https://htmlpreview.github.io/?https://github.com/GfellerLab/SuperCell/blob/supercell-2.0/docs/tutorials/SuperCell2.0_BM_CITE_seq.html)
2.  [Building and analyzing metacells in PBMC 10X multiome data with SuperCell2.0](https://htmlpreview.github.io/?https://github.com/GfellerLab/SuperCell/blob/supercell-2.0/docs/tutorials/SuperCell2.0_PBMC_10x_multiome.html)
3.  [Building and analyzing a PBMC CITE-seq atlas with SuperCell2.0 and STACAS](https://htmlpreview.github.io/?https://github.com/GfellerLab/SuperCell/blob/supercell-2.0/docs/tutorials/SuperCell2.0_PBMC_CITE_seq_atlas.html)

## [License]

SuperCell2.0 is developed by the group of David Gfeller at University of Lausanne.

SuperCell2.0 is available under GPL-3 License.

For scientific questions, please contact Léonard Hérault ([leonard.herault\@gustaveroussy.fr](mailto:leonard.herault@gustaveroussy.fr)) or David Gfeller ([David.Gfeller\@unil.ch](mailto:David.Gfeller@unil.ch)).

## How to cite

If you use SuperCell2.0 in a publication, please cite:

-   [Hérault et al. SuperCell2.0 enables semi-supervised construction of multimodal metacell atlases](https://doi.org/10.64898/2026.02.19.706848)

-   [Bilous et al. Metacells untangle large and complex single-cell transcriptome networks, BMC Bioinformatics (2022).](https://doi.org/10.1186/s12859-022-04861-1)

-   [Bilous et al. Building and analyzing metacells in single-cell genomics data, Mol Syst Bio (2024).](https://doi.org/10.1038/s44320-024-00045-6)
