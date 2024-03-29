# Test data load

data("cell_lines")

true_dims <- c(2000, 3918) # true dimension of GE matrix (i.e., 2000 genes and 3918 cells)
true_cell_ids <- c("Lib90_00000", "Lib90_00001", "Lib90_00002", "Lib90_00003", "Lib90_00004",
                   "Lib90_00005", "Lib90_00006", "Lib90_00007", "Lib90_00008","Lib90_00009")
true_gene_names <- c("PIGG", "KAT7", "CGB7", "CDS1", "POLR2M", "FH", "STAT5B", "COMMD10", "SERINC3", "LAMA3")

true_meta <- c("HCC827", "HCC827", "H838",   "HCC827", "HCC827", "HCC827", "H1975",  "H1975",  "HCC827", "HCC827")

test_that("Demo data (cell_lines) loads correctly", {
  expect_equal(names(cell_lines), c("GE", "meta"))
  expect_equal(colnames(cell_lines$GE)[1:10], true_cell_ids)
  expect_equal(rownames(cell_lines$GE)[1:10], true_gene_names, info = "Inconsistent gene names")
  expect_equal(cell_lines$meta[1:10], true_meta)
})


# Test metacell construction with `SCimplify()`
true_gamma <- 20
true_k.knn <- 5

load(file = file.path(paste0("../output/SC_cell_line_gamma_", true_gamma,".Rda"))) # (true) SC
load(file = file.path(paste0("../output/SC_GE_first_10_cell_line_gamma_", true_gamma, ".Rda")))  # (true) SC.GE.first.10
load(file = file.path(paste0("../output/SC_PCA_cell_line_gamma_", true_gamma, ".Rda")))  # (true) SC.PCA
load(file = file.path(paste0("../output/SC_markers_cell_line_gamma_", true_gamma, ".Rda"))) # (true) markers.all.positive

SC_test <- SCimplify(
  X = cell_lines$GE,
  gamma = true_gamma,
  k.knn = true_k.knn,
  n.var.genes = 1000
)

test_that("SCimplify of demo data (cell_lines) work correctly (as before)", {
  expect_equal(names(SC_test)[1:10], names(SC)[1:10], info = "SC names have been changed")
  expect_equal(SC_test$gamma, SC$gamma, info = "Inconsistent graning level")
  expect_equal(SC_test$genes.use, SC$genes.use, info = "Inconsistent set of features used to build metacells")
  expect_equal(SC_test$do.approx, SC$do.approx, info = "Approximate or exact metacell construction is not preserved")
  expect_equal(SC_test$n.pc, SC$n.pc, info = "Inconsistent number of PC dimension used for metacell construction")
  expect_equal(SC_test$N.SC, SC$N.SC, info = "Inconsistent number of metacells")
  expect_equal(SC_test$membership, SC$membership, info = "Metacell membweship vector is not preserved")
})

# Test metacell construction with `SCimplify_from_embedding()`
# first get embedding used in `SCimplify()`
X <- cell_lines$GE
colnames(X) <- colnames(cell_lines$GE)
rownames(X) <- rownames(cell_lines$GE)
gamma <- true_gamma
k.knn <- true_k.knn
n.var.genes <- 1000
genes.use <- genes.exclude <- NULL
fast.pca <- TRUE
n.pc <- 10
do.scale <- TRUE
seed <- 12345

N.c <- ncol(X)


if(is.null(rownames(X))){
  if(!(is.null(genes.use) | is.null(genes.exclude))){
    stop("rownames(X) is Null \nGene expression matrix X is expected to have genes as rownames")
  } else {
    warning("colnames(X) is Null, \nGene expression matrix X is expected to have genes as rownames! \ngenes will be created automatically in a form 'gene_i' ")
    rownames(X) <- paste("gene", 1:nrow(X), sep = "_")
  }
}

if(is.null(colnames(X))){
  warning("colnames(X) is Null, \nGene expression matrix X is expected to have cellIDs as colnames! \nCellIDs will be created automatically in a form 'cell_i' ")
  colnames(X) <- paste("cell", 1:N.c, sep = "_")
}

cell.ids <- colnames(X)

genes.use <- SC_test$genes.use
X <- X[genes.use,]


presampled.cell.ids <- cell.ids
rest.cell.ids       <- c()


X.for.pca            <- Matrix::t(X[genes.use, presampled.cell.ids])
if(do.scale){ X.for.pca            <- scale(X.for.pca) }
X.for.pca[is.na(X.for.pca)] <- 0

if(is.null(n.pc[1]) | min(n.pc) < 1){stop("Please, provide a range or a number of components to use: n.pc")}
if(length(n.pc)==1) n.pc <- 1:n.pc

if(fast.pca & (N.c < 1000)){
  warning("Normal pca is computed because number of cell is low for irlba::irlba()")
  fast.pca <- FALSE
}

if(!fast.pca){
  PCA.presampled          <- stats::prcomp(X.for.pca, rank. = max(n.pc), scale. = F, center = F)
} else {
  set.seed(seed)
  PCA.presampled          <- irlba::irlba(X.for.pca, nv = max(n.pc, 25))
  PCA.presampled$x        <- PCA.presampled$u %*% diag(PCA.presampled$d)
}

rownames(PCA.presampled$x) <- cell.ids

# build metacells using the same PCA as one used in `SCimplify`
SC_test_from_emb <- SCimplify_from_embedding(
  X = PCA.presampled$x[,n.pc],
  gamma = true_gamma,
  k.knn = true_k.knn
)

test_that("SCimplify_from_embedding of demo data (cell_lines) work correctly (as SCimplify)", {
  expect_equal(names(SC_test_from_emb)[1:10], names(SC)[1:10], info = "SC names have been changed")
  expect_equal(SC_test_from_emb$gamma, SC$gamma, info = "Inconsistent graning level")
  expect_equal(SC_test_from_emb$do.approx, SC$do.approx, info = "Approximate or exact metacell construction is not preserved")
  expect_equal(SC_test_from_emb$n.pc, SC$n.pc, info = "Inconsistent number of PC dimension used for metacell construction")
  expect_equal(SC_test_from_emb$N.SC, SC$N.SC, info = "Inconsistent number of metacells")
  expect_equal(SC_test_from_emb$membership, SC$membership, info = "Metacell membweship vector is not preserved")
})

# Test averaging gene expression

SC_test.GE <- supercell_GE(cell_lines$GE, groups = SC_test$membership, mode = "average")

test_that("Metacell gene expression of gemo data (cell_lines) works as before", {
  expect_equal(SC_test.GE[1:10,], SC.GE.first10)
})

# Test metacell PCA

SC_test.PCA <- supercell_prcomp(
  Matrix::t(SC_test.GE),
  genes.use = SC_test$genes.use,
  supercell_size = SC_test$supercell_size,
  k = 20
)

test_that("Metacell PCA of demo data (cell_lines) works as before", {
  expect_equal(abs(SC_test.PCA$x), abs(SC.PCA$x)) # on Windows, the sign of some PCs is opposite, so also test distance next
  expect_equal(dist(SC_test.PCA$x[1:10,]), dist(SC.PCA$x[1:10,]), ignore_attr = TRUE)
})

# Test metacell assignment
SC_test$cell_line <- supercell_assign(supercell_membership = SC_test$membership, clusters = cell_lines$meta, method = "jaccard")

test_that("Metacell assignment works as before", {
  expect_equal(SC_test$cell_line, SC$cell_line)
})

# Test metacell clustering
D                  <- dist(SC_test.PCA$x)
SC_test$clustering <- supercell_cluster(D = D, k = 5, supercell_size = SC_test$supercell_size)$clustering

test_that("Clustering of demo data (cell_lines) works as before", {
  expect_equal(SC_test$clustering, SC$clustering)
})

# Test DEA of metacells
## mapping metacell cluster to cell line
map.cluster.to.cell.line    <- supercell_assign(supercell_membership = SC_test$clustering, clusters  = SC_test$cell_line)
## clustering as cell line
SC_test$clustering_reordered     <- map.cluster.to.cell.line[SC_test$clustering]
markers_test.all.positive        <- supercell_FindAllMarkers(
  ge = SC_test.GE,
  supercell_size = SC_test$supercell_size,
  clusters = SC_test$clustering_reordered,
  logfc.threshold = 1,
  only.pos = T)

test_that("DEA results of demo data (cell_lines) are preserved", {
  expect_equal(names(markers_test.all.positive), names(markers.all.positive))
  expect_equal(markers_test.all.positive[["A549"]], markers.all.positive[["A549"]])
  expect_equal(markers_test.all.positive[["H1975"]], markers.all.positive[["H1975"]])
  expect_equal(markers_test.all.positive[["H2228"]], markers.all.positive[["H2228"]])
  expect_equal(markers_test.all.positive[["H838"]], markers.all.positive[["H838"]])
  expect_equal(markers_test.all.positive[["HCC827"]], markers.all.positive[["HCC827"]])
})

