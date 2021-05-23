#' compute PCA for super-cell data (sample-weighted data)
#'
#' @param X super-cell transposed gene expression matrix (! where rows represent super-cells and cols represent genes)
#' @param genes.use genes to use for dimensionality reduction
#' @param genes.exclude genes to exclude from dimensionaloty reduction
#' @param supercell_size a vector with supercell sizes (ordered the same way as in X)
#' @param k number of components to compute
#' @param do.scale scale data before PCA
#' @param do.center center data before PCA
#' @param fast.pca whether to run fast PCA (works for datasets with |super-cells| > 50)
#' @param seed a seed to use for \code{set.seed}
#'
#' @return the same object as \link[stats]{prcomp} result
#' @export
#'
#'


supercell_prcomp <- function(X,
                             genes.use = NULL,
                             genes.exclude = NULL,
                             supercell_size = NULL,
                             k = 20,
                             do.scale = TRUE,
                             do.center = TRUE,
                             double.centering = FALSE, # depricated
                             fast.pca = TRUE, #
                             seed = 12345) {
  if(is.null(supercell_size)){
    supercell_size = rep(1, nrow(X))
  }

  if(is.null(genes.use)) genes.use <- colnames(X)
  if(!is.null(genes.exclude)) genes.use <- setdiff(genes.use, genes.exclude)

  X <- X[, genes.use]

  if(do.scale | do.center){
    X <- corpcor::wt.scale(X, w = supercell_size, center = do.center, scale = do.scale)
  }

  W         <- diag(supercell_size)

  # X for sample-weighted PCA
  X.for.PCA <- (1/sum(supercell_size)) * t(X) %*% W %*% X

  set.seed(seed)
  if(!fast.pca){
    pca                   <- prcomp(x = X.for.PCA, center = FALSE, scale. = FALSE, rank. = k)
    pca$center            <- do.center
    pca$scale             <- do.scale
  } else {
    pca                   <- irlba::irlba(X.for.PCA, nv = k) # is equvivalent to prcomp when there is no centering and scaling (then prcomp is the same as svd)
    pca$rotation          <- pca$v
    pca$v                 <- NULL
    pca$center            <- do.center
    pca$scale             <- do.scale
    pca$sdev              <- sqrt(pca$d)
  }


  pca$x <- X %*% pca$rotation
  return(pca)
}


## irlba::irlba is a good approximation of svd and svd is equvivalent to prcomp, when center = F and scale = F,
## thus, in all my cases, irlba is an approximartion of prcomp, because I usualy do scaling in advance
## irlba has less memory allocation and is faster

test_irlba_and_prcomp <- function(){
  if(0){
    set.seed(12345)
    X <- readRDS("./data/5cancer_cell_lines_10x_GE.Rds")

    var.genes <- apply(X, 1, var)
    var.genes <- order(var.genes, decreasing = T)[1:500]

    Xt <- t(X[var.genes,])
    #Xt <- scale(Xt, center = T, scale = T)


    t.pr <- bench::mark(pr <- prcomp(Xt, rank. = 10, center = F, scale. = F))
    d.pr <- dist(pr$x)
    plot(x = pr$x[,1], y = pr$x[,2])

    t.sv <- bench::mark("svd" = {sv <- svd(Xt, nu = 10)
    sv$x <- sv$u %*% diag(sv$d[1:10])})
    d.sv <- dist(sv$x)
    plot(x = sv$x[,1], y = sv$x[,2])

    summary(abs(d.sv - d.pr))

    t.ir.pr <- bench::mark(ir.pr <- irlba::prcomp_irlba(Xt, n = 10))
    t.ir.svd <- bench::mark(ir.svd <- irlba::irlba(Xt, nv = 10))

    ir <- irlba::irlba(Xt, nv = 10)
    ir$x <- ir$u %*% diag(ir$d)
    d.ir <- dist(ir$x)

    plot(x = ir$x[,1], y = ir$x[,2])

    summary(abs(d.ir - d.pr))

    pca <- supercell_prcomp(t(X[1:100,]), supercell_size = sample(1:10, ncol(X), replace = T), double.centering = F)
    pca.dc <- supercell_prcomp(t(X[1:1000,]), supercell_size = sample(1:10, ncol(X), replace = T))
  }
}
