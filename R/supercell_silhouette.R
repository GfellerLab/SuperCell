#' Compute Silhouette index accounting for samlpe size (super cells size) ###
#'
#' @param x -- clustering
#' @param dist - distance among super-cells
#' @param supercell_size -- super-cell size
#'
#' @return silhouette result
#' @export


supercell_silhouette <- function(x, dist, supercell_size = NULL){
  a <- c()
  b <- c()
  s <- c()
  dist <- as.matrix(dist)

  n.cl   <- length(unique(x))
  n.elem <- length(x)
  clusters <- unique(x)

  if(is.null(supercell_size)){
    supercell_size <- rep(1, n.elem)
  }

  for(i in 1:n.elem){
    Ci      <- which(x == x[i])
    nCi     <- which(x != x[i])
    size.Ci <- sum(supercell_size[Ci])
    # print(paste0("size.Ci[",i, "]:"))
    # print(size.Ci)

    a.cur   <- 0
    for(j in Ci){
      a.cur <- a.cur + dist[i, j]* supercell_size[j]
    }
    if(size.Ci > supercell_size[i]){
      a.cur <- a.cur/(size.Ci-supercell_size[i])
    } else {
      a.cur <- 0
    }
    a <- c(a, a.cur)
  #  print("a.cur:")
  #  print(a.cur)

    b.cur   <- c() ## compute distance from a cell i to centers of all other clusters

    for(cl in clusters){
      b.curj  <- Inf
      if(cl != x[i]){  ## if not a current cluster
        b.curj   <- 0
        Cj       <- which(x == cl)
        size.Cj  <- sum(supercell_size[Cj])
        for(j in Cj){
          b.curj <- b.curj + dist[i, j] * supercell_size[j]
        }
        b.curj <- b.curj/size.Cj
      }
      b.cur <- c(b.cur, b.curj)
    }
    b.cur <- min(b.cur, na.rm = TRUE)
    b <- c(b, b.cur)

  #  print("b.cur:")
  #  print(b.cur)


    if(size.Ci > supercell_size[i]){
      s <- c(s, (b.cur - a.cur)/(max(a.cur, b.cur)))
    } else {
      s <- c(s, 0)
    }

  }

  s.matrix <- cbind(x, s)
  colnames(s.matrix) <- c("cluster", "silhouette width")

  #### cluster avg width
  clus.avg.widths        <- rep(-1, n.cl)
  names(clus.avg.widths) <- unique(x)
  for(cl in unique(x)){
    Ci      <- which(x == cl)
    size.Ci <- sum(supercell_size[Ci])

    clus.avg.widths[cl] <- 0
    for(i in Ci){
      clus.avg.widths[cl] <- clus.avg.widths[cl] + supercell_size[i] * s[i]
    }
    clus.avg.widths[cl] <- clus.avg.widths[cl]/size.Ci
  }


  result <- list(s = s.matrix,
                 clus.avg.widths = clus.avg.widths,
                 avg.width = weighted.mean(s, supercell_size))
  return(result)

}



