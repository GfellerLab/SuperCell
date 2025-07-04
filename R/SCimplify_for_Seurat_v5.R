#' SCimplify_for_Seurat
#'
#' \code{SCimplify_for_Seurat} 
#' Build metacells from a Seurat single-cell object. 
#' @param seurat A Seurat single-cell object. It has to be preprocessed (eg. latent space computed) for the assay(s) used to identify metacells
#' @param sobj.mc A metacell seurat object that will be rescaled (optional) it requires a metacell_hierarchy object in the slot misc.
#' @param gamma graining level.
#' @param assay a list of one or two assays to use to build the knn graph on which metacell are identified. 
#' @param reduction a list of corresponding reduction name in the seurat single cell object.
#' @param dims a list of corresponding dimensions to use.
#' @param membership a vector of metacell membership as in SuperCell v1 to use to directly aggregate the data (optionnal).
#' @return  A Seurat metacell object with metacell_hierarchy and memebrships in the slot misc.
#' @examples
#' sobj.mc <- SCimplify_for_Seurat(seurat = pbmc,
#'                          gamma = 30)
#' @import Seurat
#' @export

SCimplify_for_Seurat <- function(seurat, 
                                 seurat.mc = NULL,
                                 k.knn = 30, 
                                 kith = NULL, 
                                 kernel = T, 
                                 gamma = 20,
                                 graph.name = NULL, 
                                 assay = c("RNA"), 
                                 reduction = list("pca"), 
                                 dims = list(c(1:30)), 
                                 membership = NULL, 
                                 metacellNormalization = F, 
                                 avg.in.data = F,
                                 fragmentFiles = NULL, 
                                 tmpPath = NULL, 
                                 outputDirMcFragment = NULL, 
                                 bgzip_path = NULL,
                                 tabix_path = NULL,
                                 prefixMC = "",
                                 peakSep = c("-", "-"), 
                                 label = NULL, 
                                 return.seurat = T, 
                                 nb_cl = NULL) 
{
  if (!is.null(label)) {
    seurat[[paste0(label,"_with_unknown")]] <- seurat[[label]]
    seurat[[paste0(label,"_with_unknown")]][is.na(seurat[[paste0(label,"_with_unknown")]])] <- "unknown"
  }
  if (is.null(seurat.mc) & is.null(membership)) {
    if (length(assay) == 1) {
      if (is.null(graph.name)) {graph.name = "nn"}
      if (is.null(label)) {
        graph <- ComputeUnimodalKnn(seurat = seurat, 
                                       k.knn = k.knn, kith = kith, kernel = kernel, 
                                       graph.name = graph.name, assay = assay, reduction = reduction, 
                                       dims = dims)
      }
      else {
        if (length(which(is.na(seurat[[label]][, 1]))) > 
            0) {
          print("using partial annotation")
          unknowns <- colnames(seurat)[which(is.na(seurat[[label]][, 
                                                                   1]))]
          graphUnknown <- ComputeUnimodalKnn(seurat = seurat, 
                                                k.knn = k.knn, kith = kith, kernel = kernel, 
                                                graph.name = graph.name, assay = assay, reduction = reduction, 
                                                dims = dims)
          igraph::V(graphUnknown)$label <- seurat[[label]][, 
                                                           1]
          unknownAndNeighbors <- unique(c(unknowns, unlist(sapply(which(is.na(seurat[[label]][, 
                                                                                              1])), FUN = function(X) {
                                                                                                names(igraph::neighbors(graphUnknown, X))
                                                                                              }))))
          graphUnknown <- igraph::subgraph(graphUnknown, 
                                           unknownAndNeighbors)
          edgeList <- igraph::as_edgelist(graphUnknown)
          igraph::E(graphUnknown)$kept <- is.na(igraph::V(graphUnknown)[edgeList[, 
                                                                                 1]]$label) | is.na(igraph::V(graphUnknown)[edgeList[, 
                                                                                                                                     2]]$label)
          graphUnknown = igraph::subgraph.edges(graphUnknown, 
                                                which(igraph::E(graphUnknown)$kept), delete.vertices = FALSE)
        }
        graphList <- lapply(X = as.vector(na.exclude(unique(seurat[[label]][, 
                                                                            1]))), FUN = function(X) {
                                                                              ComputeUnimodalKnn(seurat = seurat, k.knn = k.knn, 
                                                                                                    kith = kith, kernel = kernel, graph.name = graph.name, 
                                                                                                    assay = assay, reduction = reduction, dims = dims, 
                                                                                                    label = label, subsetLabel = X)
                                                                            })
        if (length(which(is.na(seurat[[label]][, 1]))) > 
            0) {
          graphList[[length(graphList) + 1]] <- graphUnknown
        }
        if (length(graphList) > 1) {
          graph <- do.call(igraph::union, graphList)
          allWeigths <- lapply(X = 1:length(graphList), 
                               FUN = function(X) {
                                 igraph::get.edge.attribute(graph, name = paste0("weight_", 
                                                                                 X))
                               })
          dfWeights <- do.call("cbind", allWeigths)
          igraph::E(graph)$weight <- rowSums(dfWeights, 
                                             na.rm = T)
          graph = igraph::permute(graph, match(igraph::V(graph)$name, 
                                               colnames(seurat)))
        } else {
          graph <- graphList[[1]]
        }
      }
    }
    else {
      if (is.null(graph.name)) {graph.name = "knn"}
      if (is.null(label)) {
        graph <- ComputeMultimodalKnn(seurat = seurat, 
                                      k.knn = k.knn, kith = kith, kernel = kernel, 
                                      graph.name = graph.name, assay = assay, reduction = reduction, 
                                      dims = dims)
      }
      else {
        if (length(which(is.na(seurat[[label]][, 1]))) > 
            0) {
          print("using partial annotation")
          unknowns <- colnames(seurat)[which(is.na(seurat[[label]][, 
                                                                   1]))]
          graphUnknown <- ComputeMultimodalKnn(seurat = seurat, 
                                               k.knn = k.knn, kith = kith, kernel = kernel, 
                                               graph.name = graph.name, assay = assay, reduction = reduction, 
                                               dims = dims)
          igraph::V(graphUnknown)$label <- seurat[[label]][, 
                                                           1]
          unknownAndNeighbors <- unique(c(unknowns, unlist(sapply(which(is.na(seurat[[label]][, 
                                                                                              1])), FUN = function(X) {
                                                                                                names(igraph::neighbors(graphUnknown, X))
                                                                                              }))))
          graphUnknown <- igraph::subgraph(graphUnknown, 
                                           unknownAndNeighbors)
          edgeList <- igraph::as_edgelist(graphUnknown)
          igraph::E(graphUnknown)$kept <- is.na(igraph::V(graphUnknown)[edgeList[, 
                                                                                 1]]$label) | is.na(igraph::V(graphUnknown)[edgeList[, 
                                                                                                                                     2]]$label)
          graphUnknown = igraph::subgraph.edges(graphUnknown, 
                                                which(igraph::E(graphUnknown)$kept), delete.vertices = FALSE)
        }
        graphList <- lapply(X = na.exclude(unique(seurat[[label]][, 
                                                                  1])), FUN = function(X) {
                                                                    ComputeMultimodalKnn(seurat = seurat, k.knn = k.knn, 
                                                                                         kith = kith, kernel = kernel, graph.name = graph.name, 
                                                                                         assay = assay, reduction = reduction, dims = dims, 
                                                                                         label = label, subsetLabel = X)
                                                                  })
        if (length(which(is.na(seurat[[label]][, 1]))) > 
            0) {
          graphList[[length(graphList) + 1]] <- graphUnknown
        }
        if (length(graphList) > 1) {
          graph <- do.call(igraph::union, graphList)
          allWeigths <- lapply(X = 1:length(graphList), 
                               FUN = function(X) {
                                 igraph::edge_attr(graph, name = paste0("weight_", 
                                                                        X))
                               })
          dfWeights <- do.call("cbind", allWeigths)
          if (any(dfWeights == -1, na.rm = T)) {
            medWeights <- median(na.exclude(dfWeights[dfWeights > 
                                                        0]))
            dfWeights[which(dfWeights == -1)] <- medWeights
          }
          igraph::E(graph)$weight <- rowSums(dfWeights, 
                                             na.rm = T)
          graph = igraph::permute(graph, match(igraph::V(graph)$name, 
                                               colnames(seurat)))
        }
        else {
          graph <- graphList[[1]]
        }
      }
    }
    walktrap <- igraph::cluster_walktrap(graph)
    seurat[[paste0("walktrap_clusters_", assay[[1]])]] <- walktrap$membership
    membership <- igraph::cut_at(walktrap, no = floor(ncol(seurat)/gamma))
    names(membership) <- colnames(seurat)
    print("metacells identified")
  }
  else {
    if (is.null(membership) & !is.null(seurat.mc)) {
      walktrap <- seurat.mc@misc$metacells_hierarchy
      membership <- igraph::cut_at(walktrap, no = floor(ncol(seurat)/gamma))
      names(membership) <- colnames(seurat)
    }
    else {
      if (!is.null(membership)) {
        walktrap <- list(membership = membership)
        gamma = floor(length(membership)/length(unique(membership)))
      }
    }
  }
  seurat[[paste0("metacell_g", gamma)]] <- membership
  if (return.seurat) {
    assaysToAgg <- Assays(seurat)[sapply(X = Assays(seurat), 
                                         FUN = function(X) {
                                           !is.null(colnames(seurat[[X]]$counts))
                                         })]
    
    isChromAssay <- sapply(X = assaysToAgg, FUN = function(X) {
      is(GetAssay(seurat,assay = X))[1] == "ChromatinAssay"
    })
    
    chrom.assay.list <- list()
    
    for (chromAssay in assaysToAgg[isChromAssay]) {
      if (avg.in.data) {
        chrom.assay.list[[chromAssay]] <- CreateChromatinAssay(counts =  MetacellExpression(seurat, 
                                                                                            assays = chromAssay, 
                                                                                            group.by = paste0("metacell_g", gamma), 
                                                                                            #layer = "counts", 
                                                                                            return.seurat = F)[[chromAssay]],
                                                               genome = genome(seurat[[chromAssay]]),
                                                               ranges = Signac::StringToGRanges(rownames(seurat[[chromAssay]]), 
                                                                                                sep = peakSep), 
                                                               annotation = Signac::Annotation(seurat[[chromAssay]]))
        chrom.assay.list[[chromAssay]]$data <- MetacellExpression(seurat, 
                                                                  assays = chromAssay, 
                                                                  pb.method = "average",
                                                                  group.by = paste0("metacell_g", gamma), 
                                                                  layer = "data", 
                                                                  return.seurat = F)[[chromAssay]]
      } else {
        chrom.assay.list[[chromAssay]] <- CreateChromatinAssay(counts =  MetacellExpression(seurat, 
                                                                                            assays = chromAssay, 
                                                                                            group.by = paste0("metacell_g", gamma), 
                                                                                            #layer = "counts", 
                                                                                            return.seurat = F)[[chromAssay]],
                                                               genome = genome(seurat[[chromAssay]]),
                                                               ranges = Signac::StringToGRanges(rownames(seurat[[chromAssay]]), 
                                                                                                sep = peakSep), 
                                                               annotation = Signac::Annotation(seurat[[chromAssay]]))
      }
      
      
      if (!is.null(fragmentFiles[[chromAssay]])) {
        if (is.null(tmpPath)) {
          tmpPath <- "./tmp/"
        }
        mcfragmentFileName <- transform_fragment_file_parallel(input_file = fragmentFiles[[chromAssay]], 
                                                               prefixMC = prefixMC, 
                                                               tmp_path = tmpPath, 
                                                               output_path = outputDirMcFragment, 
                                                               membership = membership, 
                                                               returnOutputFileName = T, 
                                                               bgzip_path = bgzip_path,
                                                               tabix_path = tabix_path,
                                                               nb_cl = nb_cl)
        print("Fragment file aggregated")
        mcFragments <- CreateFragmentObject(mcfragmentFileName,
                                            cells = colnames(chrom.assay.list[[chromAssay]]))
        
        Fragments(chrom.assay.list[[chromAssay]]) <- mcFragments
        # metacell.name <- paste0("Metacell_",c(1:ncol(data.return[[i]])))
        
      }
      
    }
    
    if (length(which(!isChromAssay)) >0 ) {
      if (avg.in.data) {
        std.assay.list <- list()
        
        for (assay in assaysToAgg[!isChromAssay]) {
          std.assay.list[[assay]] <- CreateAssay5Object(counts = MetacellExpression(seurat, 
                                                                                    assays = assay, 
                                                                                    group.by = paste0("metacell_g", gamma), 
                                                                                    #layer = "counts", 
                                                                                    return.seurat = F)[[assay]],
                                                        data =  MetacellExpression(seurat, 
                                                                                   assays = assay, pb.method = "average",
                                                                                   group.by = paste0("metacell_g", gamma), 
                                                                                   layer = "data", 
                                                                                   return.seurat = F)[[assay]]
          )
        }
        
        seurat.mc <- CreateSeuratObject(std.assay.list[[1]],assay = names(std.assay.list)[1])
        for (std.a in names(std.assay.list)[-1]) {
          seurat.mc[[std.a]] <- std.assay.list[[std.a]]
        }
        for (a in names(chrom.assay.list)) {
          seurat.mc[[a]] <- chrom.assay.list[[a]]
        }
      } else {
        seurat.mc <- MetacellExpression(seurat, 
                                        assays = assaysToAgg[!isChromAssay], 
                                        group.by = paste0("metacell_g", gamma), 
                                        #layer = "counts", 
                                        return.seurat = T)
      }
      for (a in names(chrom.assay.list)) {
        seurat.mc[[a]] <- chrom.assay.list[[a]]
      }
      
      
    } else {
      
      seurat.mc <- CreateSeuratObject(chrom.assay.list[[1]],assay = names(chrom.assay.list)[1])
      for (a in names(chrom.assay.list)[-1]) {
        seurat.mc[[a]] <- chrom.assay.list[[a]]
      }
      
    }
    if (metacellNormalization) {
      normalizations <- names(seurat@commands)[startsWith(names(seurat@commands), 
                                                          prefix = "NormalizeData")]
      for (n in normalizations) {
        DefaultAssay(seurat.mc) <- seurat@commands[[n]]$assay
        seurat.mc <- NormalizeData(object = seurat.mc, 
                                   normalization.method = seurat@commands[[n]]$normalization.method, 
                                   scale.factor = seurat@commands[[n]]$scale.factor, 
                                   margin = seurat@commands[[n]]$margin, verbose = TRUE)
      }
    }
  }
  else {
    seurat.mc <- list(membership = membership, supercell_size = as.numeric(table(membership)), 
                      h_membership = walktrap)
  }
  fields <- sapply(X = colnames(seurat@meta.data), FUN = function(X) {
    is.character(seurat[[X]][, 1]) | is.factor(seurat[[X]][, 
                                                           1])
  })
  print("metadata assignement")
  for (f in colnames(seurat@meta.data)[fields]) {
    assign_res <- supercell_assign(clusters = seurat[[f]][,1], 
                                   supercell_membership =  paste0("Metacell_",membership), method = "absolute")
    #names(assign_res) <- paste0("g", names(assign_res))
    seurat.mc[[f]] <- assign_res
    purity_res <- supercell_purity(clusters = seurat[[f]][, 1], supercell_membership = paste0("Metacell_",membership))
    #names(purity_res) <- paste0("g", names(purity_res))
    seurat.mc[[paste0(f, "_purity")]] <- purity_res
  }
  if (!is.null(label)) {
    # annotate metacell containing only unknown cell as unknown
    seurat.mc[[label]][seurat.mc[[paste0(label,"_purity")]]==0] <- "unknown"
  }
  if (return.seurat) {
    seurat.mc$size <- as.numeric(table(membership))
    seurat.mc@misc$metacells_hierarchy <- walktrap
    seurat.mc@misc$walktrap_clusters <- walktrap$membership
    seurat.mc@misc$gamma <- gamma
    seurat.mc@misc$membership <- membership
  }
  return(seurat.mc)
}


#' SCimplify_for_Seurat_v5
#'
#' \code{SCimplify_for_Seurat_v5} 
#' Copy of SCimplify_for_Seurat for nmanuscript workflow compatibility
#' Build metacells from a Seurat single-cell object. 
#' @param seurat A Seurat single-cell object. It has to be preprocessed (eg. latent space computed) for the assay(s) used to identify metacells
#' @param sobj.mc A metacell seurat object that will be rescaled (optional) it requires a metacell_hierarchy object in the slot misc.
#' @param gamma graining level.
#' @param assay a list of one or two assays to use to build the knn graph on which metacell are identified. 
#' @param reduction a list of corresponding reduction name in the seurat single cell object.
#' @param dims a list of corresponding dimensions to use.
#' @param membership a vector of metacell membership as in SuperCell v1 to use to directly aggregate the data (optionnal).
#' @return  A Seurat metacell object with metacell_hierarchy and memebrships in the slot misc.
#' @examples
#' sobj.mc <- SCimplify_for_Seurat_v5(seurat = pbmc,
#'                          gamma = 30)
#' @import Seurat
#' @export
SCimplify_for_Seurat_v5 <- SCimplify_for_Seurat


