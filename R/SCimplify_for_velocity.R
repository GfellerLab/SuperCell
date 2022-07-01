#' Construct super-cells from spliced and un-spliced matrices
#'
#' @param emat spliced (exonic) count matrix
#' @param nmat unspliced (nascent) count matrix
#' @param gamma graining level of data (proportion of number of single cells in the initial dataset to the number of super-cells in the final dataset)
#' @param membership metacell membership vector (if provided, will be used for \code{emat}, \code{nmat} metacell matrices averaging)
#' @param ... other parameters from \link{SCimplify}
#'
#' @return list containing vector of membership, spliced count and un-spliced count matrices
#' @export



SCimplify_for_velocity <- function(emat,
                                   nmat,
                                   gamma = NULL,
                                   membership = NULL,
                                   ...){
        SCim <- NULL
        if (is.null(membership)){
                if (is.null(gamma)){
                        stop("Please specify supercell membership or graining level gamma.")
                } else {

                        #print("before: counts")
                        mat <- methods::as(Matrix::as.matrix(emat) + Matrix::as.matrix(nmat), 'sparseMatrix')

                        #print("before: logp1")
                        ge <- scater::normalizeCounts(mat)
                        #ge <- log1p(sweep(mat, 2, Matrix::colSums(mat), FUN = "/") * 1e4)

                        #print("before: SCimplify")

                        SCim <- SCimplify(ge, k.knn = 5, gamma = gamma, ...)
                        membership <- SCim$membership
                }


        }


        #print("before: Create spliced and unspliced counts")
        ## Create spliced and unspliced counts
        SC.emat <- methods::as(supercell_GE(emat, membership), "sparseMatrix")
        SC.nmat <- methods::as(supercell_GE(nmat, membership), "sparseMatrix")
        colnames(SC.emat) <- paste0(1:ncol(SC.emat))
        colnames(SC.nmat) <- paste0(1:ncol(SC.nmat))
        # output
        return(list(membership = membership, emat = SC.emat, nmat = SC.nmat, SC = SCim))
}
