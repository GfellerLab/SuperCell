#' Pancreatic cell dataset
#'
#' Spliced and un-spliced scRNA-seq counts of 3696 pancreatic cells from  \href{https://dev.biologists.org/content/146/12/dev173849.abstract}{Bastidas-Ponce et al. (2018)}.
#'
#'
#' @format A list with spliced count matrix (emat), un-spliced count matrix (nmat) and metadata data frame (meta):
#' \describe{
#'   \item{emat}{spliced (exonic) count matrix}
#'   \item{nmat}{un-spliced (intronic) count matrix}
#' }
#' @source \url{https://scvelo.readthedocs.io/Pancreas.html}

"pancreas"


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
#' @source \url{https://doi.org/10.1038/s41592-019-0425-8}

"cell_lines"
