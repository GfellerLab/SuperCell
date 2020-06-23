#' compute PCA for super-cell data (sample-weighted data)
#'
#' @param X an averaged gene expression matrix (! where rows represent super-cells and cols represent genes)
#' @param supercell_size a vector with supercell size (ordered the same way as in X)
#' @param k number of components to compute
#' @param seed a seed to use for \code{set.seed}
#'
#' @return the same object as \link[stats]{prcomp} result
#' @export
#'
#'



supercell_prcomp <- function(X, genes.use = NULL, supercell_size = NULL, k = 20, do.scale = TRUE, do.center = TRUE, seed = 12345){
  if(is.null(supercell_size)){
    supercell_size = rep(1, nrow(X))
  }

  if(is.null(genes.use)) genes.use <- colnames(X)

  X <- X[, genes.use]

  if(do.scale | do.center){
    X <- corpcor::wt.scale(X, w = supercell_size, center = do.center, scale = do.scale)
  }

  n.i  <- nrow(X)
  n.j  <- ncol(X)
  m    <- rep(0, n.j)
  for(i in 1:n.i){
    m <- m + supercell_size[i]*X[i,]
  }
  m         <- m/sum(supercell_size)
  M         <- matrix(m, byrow = T, nrow = n.i, ncol = n.j) # matrix where rows are m's

  X.cent    <- X - M

  W         <- diag(supercell_size)

  X.for.PCA <- (1/sum(supercell_size)) * t(X.cent) %*% W %*% X.cent

  #print(X.for.PCA)
  set.seed(seed)
  pca       <- prcomp(x = X.for.PCA, center = FALSE, scale. = FALSE, rank. = k)

  # n.sdev    <- length(pca$sdev)
  # scale     <- sqrt(1/(sum(supercell_size) - 1))
  # for(i in 1:n.sdev){
  #   w.mean       <- weighted.mean(pca$x[,i], supercell_size)
  #   pca$sdev[i]  <- scale * sqrt(sum((supercell_size*pca$x[,i] - w.mean)^2))
  # }

  pca$x <- X.cent %*% pca$rotation

  return(pca)
}
