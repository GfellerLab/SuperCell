ComputeUnimodalKnn <- function(seurat, 
                               k.knn = 30, 
                               kith = NULL, 
                               kernel = T, 
                               graph.name = "nn", 
                               assay = c("RNA"), 
                               reduction = list("pca"), 
                               dims = list(c(1:30)), 
                               label = NULL, 
                               subsetLabel = NULL,
                               verbose = FALSE) 
{
  if (!is.null(subsetLabel) & !is.null(label)) {
    seurat <- seurat[, seurat[[label]][, 1] == subsetLabel]
    # print(seurat)
    if (dim(seurat)[2] <= 2) {
      k.knn = dim(seurat)[2]
    }
    if (k.knn > dim(seurat)[2]) {
      k.knn <- dim(seurat)[2] - 1
    }
  }
  # print(k.knn)
  if (k.knn == 1) {
    graph <- igraph::make_empty_graph(n = 1, directed = F)
  }
  else {
    seurat <- Seurat::FindNeighbors(seurat, 
                                    reduction = reduction[[1]], 
                                    dims = dims[[1]], 
                                    k.param = k.knn, 
                                    verbose = verbose,
                                    return.neighbor = T)
    if (kernel) {
      if (verbose) {message("computing kernel")}
      if (is.null(kith)) {
        kith = k.knn%/%2
      }
      if (kith == 1) {
        kith = kith + 1
      }
      if (verbose) {
        message("Using assay:")
        message(assay[[1]])
      }
      graph.name = paste0(assay[[1]], ".", graph.name)
      # print(graph.name)
      # print(k.knn)
      # print(is(seurat))
      j <- as.numeric(x = t(x = seurat@neighbors[[graph.name]]@nn.idx))
      i <- ((1:length(x = j)) - 1)%/%k.knn + 1
      x <- as.numeric(x = t(x = seurat@neighbors[[graph.name]]@nn.dist))
      graph.adj <- as(object = Matrix::sparseMatrix(i = i, 
                                                    j = j, x = x, dims = c(ncol(x = seurat), ncol(x = seurat))), 
                      Class = "Graph")
      rownames(x = graph.adj) <- colnames(x = seurat)
      colnames(x = graph.adj) <- colnames(x = seurat)
      graph <- igraph::graph_from_adjacency_matrix(graph.adj, 
                                                   diag = F, mode = "directed", weighted = T)
      graph <- igraph::as.undirected(graph, edge.attr.comb = "mean")
      graph.adj <- igraph::as_adjacency_matrix(graph, attr = "weight")
      sigmas <- seurat@neighbors[[graph.name]]@nn.dist[, 
                                                       kith]
      graph.adj = expm1(-(graph.adj %*% Matrix::Diagonal(x = 1/sigmas))^(2)) + 
        igraph::as_adjacency_matrix(graph, sparse = T)
      graph.adj <- graph.adj + Matrix::t(graph.adj)
      graph <- igraph::graph_from_adjacency_matrix(graph.adj, 
                                                   diag = F, mode = "undirected", weighted = T)
    }
    else {
      Seurat::DefaultAssay(seurat) <- assay[[1]]
      seurat <- Seurat::FindNeighbors(seurat, reduction = reduction[[1]], 
                                      dims = dims[[1]], k.param = k.knn)
      graph.name = paste0(assay[[1]], "_", graph.name)
      # print(graph.name)
      # print("making graph symmetric")
      graph.adj <- seurat@graphs[[graph.name]] + Matrix::t(seurat@graphs[[graph.name]])
      graph <- igraph::graph_from_adjacency_matrix(graph.adj, 
                                                   diag = F, mode = "undirected", weighted = T)
      if (!grepl("snn", graph.name)) {
        igraph::E(graph)$weight <- 1
      }
    }
  }
  igraph::V(graph)$name <- colnames(seurat)
  return(graph)
}


ComputeUnimodalKnn_v5 <- ComputeUnimodalKnn



ComputeMultimodalKnn <- function(seurat, 
                                 k.knn = 30, 
                                 kith = NULL, 
                                 kernel = T, 
                                 graph.name = "knn", 
                                 assay = c("RNA","ADT"), 
                                 reduction = list("pca","apca"), 
                                 dims = list(c(1:30),c(1:30)),
                                 label = NULL,
                                 subsetLabel = NULL,
                                 verbose = FALSE) {
  kernelOri <- kernel
  if (!is.null(subsetLabel) & !is.null(label)) {
    seurat <- seurat[,seurat[[label]][,1] == subsetLabel]
    if (dim(seurat)[2] <= 2) {
      k.knn = dim(seurat)[2]
    }
    if (k.knn> dim(seurat)[2]) {
      k.knn <- dim(seurat)[2] -1
    }
  }
  #print(k.knn)
  if (k.knn == 1) {
    graph <- igraph::make_empty_graph(n=1,directed = F)
  } else {
    # adapt knn.range (approximate nerghbors) computed if needed, default value in FindMultimodalNeighbors is 200)
    #
    if (dim(seurat)[2] <  200) {
      knn.range = dim(seurat)[2]
    } else {
      knn.range = 200
    }
    # print(knn.range)
    searchingMn <- T
    # sometimes knn.range need to be decrease by more than the number of cells 
    while (searchingMn) {
      
      # print("k.knn:")
      # print(k.knn)
      # print("knn.range:")
      # print(knn.range)
      searchingMn <- F
      tryCatch( { seurat <- Seurat::FindMultiModalNeighbors(seurat,
                                                            reduction = reduction,
                                                            dims = dims,
                                                            k.nn = k.knn,
                                                            knn.range = knn.range,
                                                            verbose = verbose) }
                , error = function(e) {searchingMn <<- T})
      knn.range <- knn.range - 1
      if(k.knn >= knn.range & searchingMn) {
        k.knn <- knn.range - 1
      }
      
      #if we cannot compute multimodal neighbors with seurat we make a complete graph
      if (k.knn < 1) {
        # print(subsetLabel)
        complete <- TRUE
        kernel <- FALSE
        kernelOri <- TRUE
        seurat@graphs[[paste0("w", graph.name)]] <- as.Graph(matrix(data = 1,
                                                                    nrow = ncol(seurat),
                                                                    ncol = ncol(seurat),
                                                                    dimnames = list(colnames(seurat),colnames(seurat))))
        break
      }
      
    }
    
    if (verbose) {message("multimodal neighbors found")}
    
    # FindMultiModalNeighbors Seurat function does not consider the cell itself as the first neighbor (contrary to FindNeighbors)
    
    
    if (kernel) {
      if (verbose) {message("computing kernel")}
      if (is.null(kith)) { # regarding previous comment on exact number of neighbors, this will differ a little bit from unimodal mode
        kith = k.knn%/%2
      }
      if(kith == 1) { #This is not needed regarding previous comment on exact number of neighbors contrary to the unimodal mode  
        kith = kith+1
      }
      if (verbose) {
        message("Using assay:")
        message(cat(assay))
      }
      
      
      # print(k.knn)
      graph.name <- "weighted.nn"
      j <- as.numeric(x = t(x = seurat@neighbors[[graph.name]]@nn.idx))
      i <- ((1:length(x = j)) - 1) %/% k.knn + 1
      ## seurat@neighbors$weighted.nn@nn.dist = sqrt(x = MinMax(data = (1 - multimodalAffinity) / 2, min = 0, max = 1)) see https://github.com/satijalab/seurat/blob/master/R/clustering.R
      x <- as.numeric(x = t(x = 1 - 2 * (seurat@neighbors$weighted.nn@nn.dist^2))) 
      graph.adj <- as(object = Matrix::sparseMatrix(i = i, 
                                                    j = j, 
                                                    x = x, 
                                                    dims = c(ncol(x = seurat), ncol(x = seurat))), 
                      Class = "Graph")
      rownames(x = graph.adj) <- colnames(x = seurat)
      colnames(x = graph.adj) <- colnames(x = seurat)
      
      # Using symetrization here and then undirected mode to create the graph is preferred to have undirected graph as input for walktrap
      graph.adj <- graph.adj + Matrix::t(graph.adj)  
      
      
      graph <- igraph::graph_from_adjacency_matrix(graph.adj, 
                                                   diag = F, mode = "undirected", weighted = T)
      
      # walktrap cannot take negative weights that can be obtained in rare cases when we retrieve multimodal kernel affinites above
      # This is when multimodal aff is zero (but the edge still retained in the knn) we retrieve a multimodal affinities just below zero
      # because of numerical precision
      # We set this weight to an epsilon value to keep the link for the walktrap
      if(length(which(igraph::E(graph)$weight < 0)) > 0) {
        
        min.weight <- 1e-16
        igraph::E(graph)$weight <- Seurat::MinMax(igraph::E(graph)$weight,min = min.weight, max = max(igraph::E(graph)$weight)) 
        
      }
    }
    else {
      graph.name = paste0("w", graph.name)
      # print(graph.name)
      graph.adj <- seurat@graphs[[graph.name]] + Matrix::t(seurat@graphs[[graph.name]])
      # Should not be needed as seurat multimodal knn looks symmetrical
      graph <- igraph::graph_from_adjacency_matrix(graph.adj, 
                                                   diag = F, mode = "undirected", weighted = T)
      
      # Remove weight from graph after symmetrization (keep it if the graph used is snn)
      # Should not be needed as seurat multimodal knn looks symmetrical
      if (!grepl("snn", graph.name)) {
        igraph::E(graph)$weight <- 1
      }
    }
    #In case a complete graph has to be computed and kernel asked, set weights to -1 and then assign them the median weight in the merging procedure
    if (kernelOri & !kernel) {
      igraph::E(graph)$weight <- -1
    }
  }
  igraph::V(graph)$name <- colnames(seurat)
  
  return(graph)
}


#' DimPlot of metacells in single cell space  
#'
#' Plots metacell in single cell space (pca, umap,..) taking the average coordinates of single cells in each metacell
#'
#'
#' @param seurat the original seurat object at the single cell level with the computed dimension reduction 
#' @param seurat.mc a seurat object at the metacell level
#' @param metacell.col name of metacell col in Seurat object (no need to use with default metacell obtained with \link{SCimplify_for_Seurat})
#' @param sc.col name of single-cell col in Seurat object 
#' @param dims numerical vector of 2 dimension components to plot
#' @param reduction computed dimension reduction in the seurat object to use
#' @param mc.color colors for metacell idents
#' @param sc.color colors for single-cell idents
#' @param alpha transparency value for the single-cell points
#' @param pt_size size the single-cell points
#' @param continuous_metric boolean indicating if the metric variable is continuous or not. If TRUE a continuous color scale will be used
#' 
#'
#'
#'@return a plot of metacells in single cell space
#'
#' @export

DimPlotSC <- function (seurat, 
                       seurat.mc, 
                       metacell.col = NULL, 
                       sc.col = NULL, 
                       dims = c(1, 2), 
                       reduction = "wnn.umap", 
                       mc.color = NULL, 
                       sc.color = NULL,
                       alpha = 1, 
                       pt_size = 0, 
                       metric = "size", 
                       continuous_metric = F) 
{
  if (!is.null(seurat.mc@misc$membership)) {
    membership <- seurat.mc@misc$membership
  }  else  {
    if (!is.null(seurat.mc@misc$gamma)) {
      gamma <- seurat.mc@misc$gamma
    }
    else {
      gamma <- floor(ncol(seurat)/ncol(seurat.mc))
    }
    membership <- igraph::cut_at(seurat.mc@misc$metacells_hierarchy, 
                                 no = floor(ncol(seurat)/gamma))
  }
  seurat$Metacell <- membership
  seuratCoord <- Embeddings(seurat[[reduction]])
  seuratCoordMetacell <- cbind(seuratCoord, membership)
  
  centroids <- stats::aggregate(seuratCoord ~ membership, seuratCoord, mean)
  rownames(centroids) <- centroids[, 1]
  #get rid of potential outliers discarded with MetaCell2
  matches <- unlist(regmatches(colnames(seurat.mc), gregexpr("[[:digit:]]+", colnames(seurat.mc))))
  
  centroids <- centroids[matches, ]
  
  centroids[[metric]] <- seurat.mc[[metric]][, 1]
  if (is.null(metacell.col)) {
    metacell.col <- "SuperCell_MC"
    centroids[[metacell.col]] <- rep("red", length(seurat.mc$size))
  }
  else {
    centroids[[metacell.col]] <- seurat.mc[[metacell.col]][, 1]
    # print(head(centroids))
  }
  if (!is.null(sc.col)) {
    seuratCoord <- data.frame(seuratCoord)
    seuratCoord[[sc.col]] <- seurat[[sc.col]][, 1]
    p <- ggplot2::ggplot(seuratCoord, aes_string(colnames(seuratCoord)[dims[1]], 
                                                 colnames(seuratCoord)[dims[2]], color = sc.col)) + 
      ggplot2::geom_point(size = pt_size, alpha = alpha)
  }
  else {
    p <- ggplot2::ggplot(data.frame(seuratCoord), aes_string(colnames(seuratCoord)[dims[1]], 
                                                             colnames(seuratCoord)[dims[2]])) + 
      ggplot2::geom_point(size = pt_size, color = "grey", alpha = 1)
    # if (!is.null(sc.color)) {
    #   print("coloring sc...")
    #   p <- p + ggplot2::scale_color_manual(values = sc.color)
    # }
    
  }
  if (!continuous_metric) {
    p <- p + ggplot2::geom_point(data = centroids, aes_string(colnames(centroids)[1 + dims[1]], colnames(centroids)[1 + dims[2]], fill = metacell.col, 
                                                              size = metric), colour = "black", pch = 21)
  }
  else {
    p <- p + ggplot2::geom_point(data = centroids, aes_string(colnames(centroids)[1 + dims[1]], colnames(centroids)[1 + dims[2]], fill = metric), 
                                 colour = "black", pch = 21, size = 2)
  }
  if (!is.null(metacell.col) & !is.null(sc.col) & !is.null(mc.color)) {
    if (metacell.col == sc.col & !is.null(mc.color)) {
      sc.color = mc.color
    }
  }
  if (!is.null(mc.color) & !continuous_metric) {
    p <- p + ggplot2::scale_fill_manual(values = mc.color) 
  }
  if (!is.null(sc.color)) {
    p <- p + ggplot2::scale_color_manual(values = sc.color)
  }
  
  p <- p +  ggplot2::theme_classic() + guides(fill = guide_legend(override.aes = list(size = 3.5)))
  return(p)
}

#' Feature-feature correlation plot with a Seurat object 
#'
#' Plots feature-feature expression and computes their correlation
#'
#' @param seurat.mc a seurat object
#' @param feature_x feature or vector of features (if vector, has to be the same lenght as feature_y)
#' @param feature_y feature or vector of features (if vector, has to be the same lenght as feature_x)
#' @param method a character string indicating which correlation coefficient is to be used for the test. One of "pearson", "kendall", or "spearman", can be abbreviated.
#' @param assays Assays to use in a vector (c(Assay_for_feature_x,Assay_for_feature_y))
#' @param clusters seurat metadata column with clustering information 
#' @param is.normalised if FALSE ADT assay will be normalized with "CLR" and RNA assay with a "LogNormalize" 
#' @param plot if FALSE do not plot and only return the computed statistics
#' @param color.use colors for idents
#' @param use.size whether to compute weighted correlation by size (required that the seurat object has a size column in its metadata)
#'
#'
#'@return a data.frame with correlation coefficient for each feature couples
#'
#' @export

supercell_FeatureFeaturePlot_Seurat <- function(seurat.mc,
                                                feature_x,
                                                feature_y,
                                                method = c("pearson", "kendall", "spearman"),
                                                assays = c("RNA","ADT"),
                                                cluster = "celltype.l2",
                                                is.normalized = F,
                                                plot = T,
                                                color.use = NULL,
                                                use.size = T) {
  method <- match.arg(arg = method)
  Seurat::DefaultAssay(seurat.mc) <- assays[1]
  if(assays[1] == "RNA") {
    if (!is.normalized) {
      seurat.mc <- Seurat::NormalizeData(seurat.mc, 
                                         normalization.method = "LogNormalize", 
                                         margin = 1)
    }
  } else {
    if (!is.normalized) {
      seurat.mc <- Seurat::NormalizeData(seurat.mc, 
                                         normalization.method = "CLR", 
                                         margin = 2)
    }
  }
  
  fe1 <- Seurat::GetAssayData(seurat.mc,slot = "data",assay = assays[1])[feature_x,]
  
  feature_x <- paste0(tolower(assays[1]),"_",feature_x)
  
  rownames(fe1) <- feature_x
  
  Seurat::DefaultAssay(seurat.mc) <- assays[2]
  if(assays[2] == "ADT") {
    if (!is.normalized) {
      seurat.mc <- Seurat::NormalizeData(seurat.mc, 
                                         normalization.method = "LogNormalize", 
                                         margin = 1)
    }
  } else {
    if (!is.normalized) {
      seurat.mc <- Seurat::NormalizeData(seurat.mc, 
                                         normalization.method = "CLR", 
                                         margin = 2)
    }
  }
  
  fe2 <- Seurat::GetAssayData(seurat.mc,slot = "data",assay = assays[2])[feature_y,]
  
  rownames(fe2) <- feature_y
  
  fe <- rbind(fe1,fe2)
  
  if(use.size) {
    sizes <- as.numeric(seurat.mc$size)
  }   else {
    sizes <-rep(1,length(ncol(seurat.mc)))
  }
  
  
  
  res <- supercell_FeatureFeaturePlot(fe,
                                      feature_x = feature_x,
                                      feature_y = feature_y,
                                      method = method,
                                      supercell_size = sizes,
                                      cluster = seurat.mc[[cluster]][,1],
                                      color.use = color.use,
                                      combine = F)
  
  
  #w.cor <- as.numeric(res$w.co)
  w.cor <- data.frame(features = names(res$w.cor),w.cor = as.numeric(res$w.cor))
  if (plot) {
    for (i in names(res$p)) {
      
      plot(res$p[[i]])
    }
  }
  
  return(w.cor)
  
}


ExpandMetacellSeurat <- function(metacell.sobj,
                                 features,
                                 assay = "RNA",
                                 #slot = "data",
                                 meta.data.vars = c('celltype',"input","method")) {
  membership <- rep(1:ncol(metacell.sobj), metacell.sobj$size)
  DefaultAssay(metacell.sobj) <- assay
  feature.data <- GetAssayData(
    object = metacell.sobj, assay = assay, slot = "data"
  )
  feature.counts <- GetAssayData(
    object = metacell.sobj, assay = assay, slot = "counts"
  )
  meta.data <- FetchData(metacell.sobj,vars = meta.data.vars)
  expanded.meta.data <- meta.data[membership,]
  
  feature.data <- feature.data[features,]
  expanded.data <-  as(feature.data[,membership],"CsparseMatrix")
  colnames(expanded.data) <- rownames(expanded.meta.data)
  
  feature.counts <- feature.counts[features,]
  expanded.counts <-  as(feature.counts[,membership],"CsparseMatrix")
  colnames(expanded.counts) <- rownames(expanded.meta.data)
  
  #colnames(expanded.data) <- colnames(expanded.data)
  expanded.sobj <- CreateSeuratObject(counts = expanded.data,meta.data = expanded.meta.data,assay = assay)
  expanded.sobj[[assay]]@counts <- expanded.counts
  Idents(expanded.sobj) <- Idents(metacell.sobj)
  return(expanded.sobj)
} 