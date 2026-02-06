#' Cancer cell lines dataset
#'
#' ScRNA-seq data of 5 cancer cell lines from [Tian et al., 2019](https://doi.org/10.1038/s41592-019-0425-8).
#'
#' Data available at authors' [GitHub](https://github.com/LuyiTian/sc_mixology/blob/master/data/) under file name *sincell_with_class_5cl.Rdata*.
#'
#' @format A list with gene expression (i.e., log-normalized counts) (GE), and metadata data (meta):
#' \describe{
#'   \item{GE}{gene expression (log-normalized counts) matrix}
#'   \item{meta}{cells metadata (cell line annotation)}
#' }
#' @source \doi{https://doi.org/10.1038/s41592-019-0425-8}

"cell_lines"


#' Subset of the bone marrow CITE-seq dataset
#'
#' CITE-seq data of one tube (HumanHTO7) of BM cells from [Stuart T, Butler A, et al. 2019](https://doi.org/10.1016/j.cell.2019.05.031).
#'
#' Data available as a seurat object in the SeuratData package' [GitHub](https://github.com/satijalab/seurat-data).
#' Data are derived from GEO: GSE128639
#'
#' @format Seurat object
#' @source \doi{https://doi.org/10.1016/j.cell.2019.05.031

"bmcite_HumanHTO7"