# DiffWeightedLM <- function(data.use, 
#                            cells.1, 
#                            cells.2,
#                            weights.1,
#                            weights.2,
#                            use.sapply = FALSE,
#                            return.all = TRUE,
#                            verbose = TRUE) 
# {
#   group <- c(rep(1, length(cells.1)), rep(2, length(cells.2)))
#   weights <- c(weights.1, weights.2)
#   
#   if (use.sapply) {
#     my.sapply <- ifelse(test = verbose, 
#                         yes = pbapply::pbsapply, no = sapply)
#     
#     p_val <- unlist(x = my.sapply(X = 1:nrow(data.use), FUN = function(x) {
#       y <- as.numeric(data.use[x, c(cells.1, cells.2)])
#       model <- lm(y ~ factor(group), weights = weights)
#       summary(model)$coefficients[2, "Pr(>|t|)"]
#     }))
#     
#     to.return <- data.frame(p_val = p_val, row.names = rownames(x = data.use))
#     
#   } else {
#     lm.test <- apply(data.use, 1, function(x) {
#       y <- as.numeric(x[c(cells.1, cells.2)])
#       model <- lm(y ~ factor(group), weights = weights)
#       summary(model)$coefficients[2, ]
#     })
#     
#     if (return.all) { 
#       to.return <- data.frame(
#         p_val = sapply(lm.test, function(res) res["Pr(>|t|)"]),
#         t_value = sapply(lm.test, function(res) res["t value"]),
#         df = sapply(lm.test, function(res) res["df"]),
#         row.names = rownames(data.use)
#       )
#     } else {
#       to.return <- data.frame(p_val = sapply(lm.test, function(res) res["Pr(>|t|)"]), 
#                               row.names = rownames(data.use))
#     }
#   }
#   return(to.return)
# }

DiffWeightedLM <- function(data.use, 
                           cells.1, 
                           cells.2,
                           weights.1,
                           weights.2,
                           use.sapply = T,
                           return.all = F,
                           verbose = TRUE) 
{
  group <- c(rep(1, length(cells.1)), rep(2, length(cells.2)))
  weights <- c(weights.1, weights.2)
  
  if (use.sapply) {
    my.sapply <- ifelse(test = verbose, 
                        yes = pbapply::pbsapply, no = sapply)
    
    lm.results <- my.sapply(X = 1:nrow(data.use), FUN = function(x) {
      y <- as.numeric(data.use[x, c(cells.1, cells.2)])
      model <- lm(y ~ factor(group), weights = weights)
      summary(model)$coefficients[2, ]  # Extract group coefficient row
    })
    
    p_val <- lm.results["Pr(>|t|)",]
    
    if (return.all) {
      to.return <- data.frame(
        p_val = p_val,
        t_value = lm.results["t value",],
        df = rep(length(cells.1) + length(cells.2) - 2, length(p_val)), # Manually set Df
        row.names = rownames(data.use)
      )
    } else {
      to.return <- data.frame(
        p_val = p_val, 
        row.names = rownames(data.use)
      )
    }
    
  } else {
    lm.results <- apply(data.use, 1, function(x) {
      y <- as.numeric(x[c(cells.1, cells.2)])
      model <- lm(y ~ factor(group), weights = weights)
      summary(model)$coefficients[2, ]  # Extract group coefficient row
    })
    
    if (return.all) {
      to.return <- data.frame(
        p_val = lm.results["Pr(>|t|)",],
        t_value = lm.results["t value",],
        df = rep(length(cells.1) + length(cells.2) - 2, length(lm.results)), # Manually set Df
        row.names = rownames(data.use)
      )
    } else {
      to.return <- data.frame(
        p_val = lm.results["Pr(>|t|)",], 
        row.names = rownames(data.use)
      )
    }
  }
  
  return(to.return)
}

DiffSurveyWeightedTTest <- function(data.use, 
                                    cells.1, 
                                    cells.2,
                                    weights.1,
                                    weights.2,
                                    use.sapply = FALSE,
                                    return.all = TRUE,
                                    verbose = TRUE) 
{
  library(survey)
  
  group <- c(rep(1, length(cells.1)), rep(2, length(cells.2)))
  weights <- c(weights.1, weights.2)
  
  if (use.sapply) {
    my.sapply <- ifelse(test = verbose, 
                        yes = pbapply::pbsapply, no = sapply)
    
    p_val <- unlist(x = my.sapply(X = 1:nrow(data.use), FUN = function(x) {
      y <- as.numeric(data.use[x, c(cells.1, cells.2)])
      design <- svydesign(ids = ~1, weights = ~weights, data = data.frame(y, group))
      res <- svyttest(y ~ factor(group), design)
      res$p.value
    }))
    
    to.return <- data.frame(p_val = p_val, row.names = rownames(x = data.use))
    
  } else {
    survey.test <- apply(data.use, 1, function(x) {
      y <- as.numeric(x[c(cells.1, cells.2)])
      design <- svydesign(ids = ~1, weights = ~weights, data = data.frame(y, group))
      res <- svyttest(y ~ factor(group), design)
      res
    })
    if (length((survey.test)) > 1) {
      to.return <- data.frame(matrix(unlist(survey.test),ncol=10,
                                     byrow = T))
      colnames(to.return) <- names(unlist(survey.test[[1]]))
      rownames(to.return) <- rownames(x = data.use)
      
      
      # t.test(x = data.use[x, cells.1], y = data.use[x, cells.2])$p.value
    }
    
    if (return.all) { 
      to.return <- data.frame(
        p_val = as.numeric(to.return$`p.value.factor(group)2`),
        t_value =  as.numeric(to.return$statistic.t),
        df =  as.numeric(to.return$parameter.df),
        row.names = rownames(data.use)
      )
    } else {
      to.return <- data.frame(p_val = res$p.value, 
                              row.names = rownames(data.use))
    }
  }
  return(to.return)
}

wtd.t.test.ess <- function(x, y = 0, weight = NULL, weighty = NULL, samedata = TRUE, 
                           alternative = "two.tailed", mean1 = FALSE, bootse = FALSE, 
                           bootp = FALSE, bootn = 1000, drops = "pairwise") 
{
  if (is.null(weight)) {
    weight <- rep(1, length(x))
  }
  if (bootse == FALSE & bootp == TRUE) 
    warning("bootp can only be used with bootstrapped standard errors")
  if (length(y) != length(x) & length(y) > 1) {
    if (samedata == TRUE) 
      warning("Treating data for x and y separately because they are of different lengths")
    samedata <- FALSE
  }
  if (length(y) == 1) 
    samedata <- FALSE
  if (samedata == TRUE & drops == "pairwise") {
    use <- !is.na(x) & !is.na(y) & !is.na(weight)
    x <- x[use]
    if (length(y) > 1) 
      y <- y[use]
    weight <- weight[use]
  }
  if (is.null(weighty) & samedata == TRUE) {
    weighty <- weight
  }
  if (is.null(weighty) & samedata == FALSE & length(y) > 1) {
    warning("y has no weights, weights for y are assumed to be 1")
    weighty <- rep(1, length(y))
  }
  if (mean1 == TRUE) {
    weight <- weight/mean(weight, na.rm = TRUE)
    if (length(y) > 1) 
      weighty <- weighty/mean(weighty, na.rm = TRUE)
  }
  n <- sum(weight[!is.na(x)])^2 / sum(weight[!is.na(x)]^2)
  # print(n)
  # n <- sum(weight[!is.na(x)], na.rm = TRUE)
  mx <- wtd.mean(x, weight, na.rm = TRUE)
  vx <- wtd.var(x, weight, na.rm = TRUE)
  if (length(y) == 1) {
    dif <- mx - y
    sx <- sqrt(vx)
    se <- sx/sqrt(n)
    if (bootse == TRUE) {
      samps <- lapply(1:bootn, function(g) sample(1:length(x), 
                                                  round(sum(weight, na.rm = TRUE), 0), replace = TRUE, 
                                                  prob = weight))
      sepests <- sapply(samps, function(q) mean(x[q], na.rm = TRUE)) - 
        y
      se <- sqrt(var(sepests))
    }
    t <- (mx - y)/se
    df <- n - 1
    p.value <- (1 - pt(abs(t), df)) * 2
    if (alternative == "greater") 
      p.value <- pt(t, df, lower.tail = FALSE)
    if (alternative == "less") 
      p.value <- pt(t, df, lower.tail = TRUE)
    if (bootp == TRUE & bootse == TRUE) 
      p.value <- 2 * min(c(sum(sepests > y & !is.na(sepests))/sum(!is.na(sepests)), 
                           sum(sepests < y & !is.na(sepests))/sum(!is.na(sepests))))
    if (bootp == TRUE & bootse == TRUE & alternative == "greater") 
      p.value <- sum(sepests > y & !is.na(sepests))/sum(!is.na(sepests))
    if (bootp == TRUE & bootse == TRUE & alternative == "less") 
      p.value <- sum(sepests < y & !is.na(sepests))/sum(!is.na(sepests))
    coef <- c(t, df, p.value)
    out2 <- c(dif, mx, y, se)
    names(coef) <- c("t.value", "df", "p.value")
    names(out2) <- c("Difference", "Mean", "Alternative", 
                     "Std. Err")
    out <- list("One Sample Weighted T-Test", coef, out2)
    names(out) <- c("test", "coefficients", "additional")
  }
  if (length(y) > 1) {
    n2 <- sum(weighty[!is.na(y)])^2 / sum(weighty[!is.na(y)]^2)
    # print(n2)
    # n2 <- sum(weighty[!is.na(y)], na.rm = TRUE)
    my <- wtd.mean(y, weighty, na.rm = TRUE)
    vy <- wtd.var(y, weighty, na.rm = TRUE)
    dif <- mx - my
    sxy <- sqrt((vx/n) + (vy/n2))
    if (bootse == TRUE) {
      samps1 <- lapply(1:bootn, function(g) sample(1:length(x), 
                                                   round(sum(weight, na.rm = TRUE), 0), replace = TRUE, 
                                                   prob = weight))
      samps2 <- lapply(1:bootn, function(g) sample(1:length(y), 
                                                   round(sum(weighty, na.rm = TRUE), 0), replace = TRUE, 
                                                   prob = weighty))
      sepests1 <- sapply(samps1, function(q) mean(x[q], 
                                                  na.rm = TRUE))
      sepests2 <- sapply(samps2, function(q) mean(y[q], 
                                                  na.rm = TRUE))
      sxy <- sqrt(var(sepests1 - sepests2, na.rm = TRUE))
    }
    df <- (((vx/n) + (vy/n2))^2)/((((vx/n)^2)/(n - 1)) + 
                                    ((vy/n2)^2/(n2 - 1)))
    t <- (mx - my)/sxy
    p.value <- (1 - pt(abs(t), df)) * 2
    if (alternative == "greater") 
      p.value <- pt(t, df, lower.tail = FALSE)
    if (alternative == "less") 
      p.value <- pt(t, df, lower.tail = TRUE)
    if (bootp == TRUE & bootse == TRUE) 
      p.value <- 2 * min(c(sum(sepests1 > sepests2 & !is.na(sepests1))/sum(!is.na(sepests1)), 
                           sum(sepests1 < sepests2 & !is.na(sepests1))/sum(!is.na(sepests1))))
    if (bootp == TRUE & bootse == TRUE & alternative == "greater") 
      p.value <- sum(sepests1 > sepests2 & !is.na(sepests1))/sum(!is.na(sepests1))
    if (bootp == TRUE & bootse == TRUE & alternative == "less") 
      p.value <- sum(sepests1 < sepests2 & !is.na(sepests1))/sum(!is.na(sepests1))
    coef <- c(t, df, p.value)
    out2 <- c(dif, mx, my, sxy)
    names(coef) <- c("t.value", "df", "p.value")
    names(out2) <- c("Difference", "Mean.x", "Mean.y", "Std. Err")
    out <- list("Two Sample Weighted T-Test (Welch)", coef, 
                out2)
    names(out) <- c("test", "coefficients", "additional")
  }
  out
}

#' Find Markers of metacells
#'
#' Same as Seurat FindMarkers function taking into account metacell size.
#' Only weighted t test using weights::wtd.t.test funnction is currently supported
#' 
#'
#' @param object Seurat metacell object
#' @param test.use onky "weighted_t" is supported.
#' @inheritParams Seurat::FindMarkers
#' @return data.frame with a ranked list of putative markers as rows, and associated statistics as columns (p-values, t statistic, degree of freedom (df) of the weighted t-test). 
#' avg_logFC, pct.1, pc.2 are computed using Seurat v5 approach (pseudocounts added at the group level) taking into account metacell sizes.
#' @import Seurat
#' @rdname FindMarkers.SuperCell
#' @export

FindMarkers.SuperCell <- function(object, ...) 
{
  UseMethod(generic = "FindMarkers.SuperCell", object = object)
}



#' Find Markers of metacells
#'
#' Same as Seurat FindMarkers function taking into account metacell size.
#' Only weighted t test using weights::wtd.t.test funnction is currently supported
#' 
#'
#' @param object Seurat metacell object
#' @param test.use onky "weighted_t" is supported.
#' @inheritParams Seurat::FindMarkers
#' @return data.frame with a ranked list of putative markers as rows, and associated statistics as columns (p-values, t statistic, degree of freedom (df) of the weighted t-test). 
#' avg_logFC, pct.1, pc.2 are computed using Seurat v5 approach (pseudocounts added at the group level) taking into account metacell sizes.
#' @import Seurat
#' @rdname FindMarkers.SuperCell
#' @concept differential_expression
#' @export
#' @method FindMarkers.SuperCell default

FindMarkers.SuperCell.default <- function(object, slot = "data", cells.1 = NULL, cells.2 = NULL, 
                                          weights.1 = NULL, weights.2 = NULL,
                                          features = NULL, logfc.threshold = 0.1, test.use = "survey_weighted_t", 
                                          min.pct = 0.01, min.diff.pct = -Inf, verbose = TRUE, only.pos = FALSE, 
                                          max.cells.per.ident = Inf, random.seed = 1, latent.vars = NULL, 
                                          min.cells.feature = 3, min.cells.group = 3, fc.results = NULL, 
                                          densify = FALSE, ...) 
{
  Seurat:::ValidateCellGroups(object = object, cells.1 = cells.1, cells.2 = cells.2, 
                              min.cells.group = min.cells.group)
  features <- features %||% rownames(x = object)
  if (test.use %in% Seurat:::DEmethods_noprefilter()) {
    features <- rownames(x = object)
    min.diff.pct <- -Inf
    logfc.threshold <- 0
  }
  alpha.min <- pmax(fc.results$pct.1, fc.results$pct.2)
  names(x = alpha.min) <- rownames(x = fc.results)
  features <- names(x = which(x = alpha.min >= min.pct))
  if (length(x = features) == 0) {
    warning("No features pass min.pct threshold; returning empty data.frame")
    return(fc.results[features, ])
  }
  alpha.diff <- alpha.min - pmin(fc.results$pct.1, fc.results$pct.2)
  features <- names(x = which(x = alpha.min >= min.pct & alpha.diff >= 
                                min.diff.pct))
  if (length(x = features) == 0) {
    warning("No features pass min.diff.pct threshold; returning empty data.frame")
    return(fc.results[features, ])
  }
  if (slot != "scale.data") {
    total.diff <- fc.results[, 1]
    names(total.diff) <- rownames(fc.results)
    features.diff <- if (only.pos) {
      names(x = which(x = total.diff >= logfc.threshold))
    }
    else {
      names(x = which(x = abs(x = total.diff) >= logfc.threshold))
    }
    features <- intersect(x = features, y = features.diff)
    if (length(x = features) == 0) {
      warning("No features pass logfc.threshold threshold; returning empty data.frame")
      return(fc.results[features, ])
    }
  }
  if (max.cells.per.ident < Inf) {
    set.seed(seed = random.seed)
    if (length(x = cells.1) > max.cells.per.ident) {
      cells.1 <- sample(x = cells.1, size = max.cells.per.ident)
    }
    if (length(x = cells.2) > max.cells.per.ident) {
      cells.2 <- sample(x = cells.2, size = max.cells.per.ident)
    }
    if (!is.null(x = latent.vars)) {
      latent.vars <- latent.vars[c(cells.1, cells.2), , 
                                 drop = FALSE]
    }
  }
  if (inherits(x = object, what = "IterableMatrix")) {
    if (test.use != "wilcox") {
      stop("Differential expression with BPCells currently only supports the 'wilcox' method.", 
           " Please rerun with test.use = 'wilcox'")
    }
    data.use <- object[features, c(cells.1, cells.2), drop = FALSE]
    groups <- c(rep("foreground", length(cells.1)), rep("background", 
                                                        length(cells.2)))
    de.results <- suppressMessages(BPCells::marker_features(data.use, 
                                                            group = groups, method = "wilcoxon"))
    de.results <- subset(de.results, de.results$foreground == 
                           "foreground")
    de.results <- data.frame(feature = de.results$feature, 
                             p_val = de.results$p_val_raw)
    rownames(de.results) <- de.results$feature
    de.results$feature <- NULL
  }
  else {
    de.results <- PerformDE.SuperCell(object = object, cells.1 = cells.1, 
                                      cells.2 = cells.2, weights.1 = weights.1, weights.2 = weights.2,
                                      features = features, test.use = test.use, 
                                      verbose = verbose, min.cells.feature = min.cells.feature, 
                                      latent.vars = latent.vars, densify = densify, ...)
  }
  if (ncol(de.results) > 1){
    de.results <- cbind(data.frame("p_val"=de.results[,1]), fc.results[rownames(x = de.results), 
                                                                       , drop = FALSE],
                        de.results[,2:ncol(de.results)])
  } else {
    de.results <- cbind(de.results, fc.results[rownames(x = de.results), 
                                                                       , drop = FALSE])                
  }
  if (only.pos) {
    de.results <- de.results[de.results[, 2] > 0, , drop = FALSE]
  }
  if (test.use %in% Seurat:::DEmethods_nocorrect()) {
    de.results <- de.results[order(-de.results$power, -de.results[, 
                                                                  1]), ]
  }
  else {
    de.results <- de.results[order(de.results$p_val, -abs(de.results$pct.1 - 
                                                            de.results$pct.2)), ]
    de.results$p_val_adj = p.adjust(p = de.results$p_val, 
                                    method = "bonferroni", n = nrow(x = object))
  }
  return(de.results)
}



PerformDE.SuperCell <- function(object, cells.1, cells.2, 
                                weights.1,weights.2,
                                features, test.use, verbose, 
                                min.cells.feature, latent.vars, densify, ...) 
{
  if (!(test.use %in% Seurat:::DEmethods_latent()) && !is.null(x = latent.vars)) {
    warning("'latent.vars' is only used for the following tests: ", 
            paste(Seurat:::DEmethods_latent(), collapse = ", "), call. = FALSE, 
            immediate. = TRUE)
  }
  if (!test.use %in% Seurat:::DEmethods_checkdots()) {
    CheckDots(...)
  }
  data.use <- object[features, c(cells.1, cells.2), drop = FALSE]
  if (densify) {
    data.use <- as.matrix(x = data.use)
  }
  
  de.results <- switch(EXPR = test.use, weighted_t = DiffWeightedTTest(data.use = data.use,
                                                                       cells.1 = cells.1, cells.2 = cells.2, 
                                                                       weights.1 = weights.1, weights.2 = weights.2,
                                                                       verbose = verbose),
                       survey_weighted_t = DiffSurveyWeightedTTest(data.use = data.use,
                                                      cells.1 = cells.1, cells.2 = cells.2, 
                                                      weights.1 = weights.1, weights.2 = weights.2,
                                                      verbose = verbose),
                       weighted_lm = DiffWeightedLM(data.use = data.use,
                                                    cells.1 = cells.1, cells.2 = cells.2, 
                                                    weights.1 = weights.1, weights.2 = weights.2,
                                                    verbose = verbose),
                       wilcox = Seurat:::WilcoxDETest(data.use = data.use, 
                                                      cells.1 = cells.1, cells.2 = cells.2, verbose = verbose, 
                                                      ...), 
                       wilcox_limma = Seurat:::WilcoxDETest(data.use = data.use, 
                                                            cells.1 = cells.1, cells.2 = cells.2, verbose = verbose, 
                                                            limma = TRUE, ...), 
                       bimod = Seurat:::DiffExpTest(data.use = data.use, 
                                                    cells.1 = cells.1, cells.2 = cells.2, verbose = verbose), 
                       roc = Seurat:::MarkerTest(data.use = data.use, cells.1 = cells.1, 
                                                 cells.2 = cells.2, verbose = verbose), 
                       t = Seurat:::DiffTTest(data.use = data.use, 
                                              cells.1 = cells.1, cells.2 = cells.2, verbose = verbose), 
                       negbinom = Seurat:::GLMDETest(data.use = data.use, cells.1 = cells.1, 
                                                     cells.2 = cells.2, min.cells = min.cells.feature, 
                                                     latent.vars = latent.vars, test.use = test.use, verbose = verbose), 
                       poisson = Seurat:::GLMDETest(data.use = data.use, cells.1 = cells.1, 
                                                    cells.2 = cells.2, min.cells = min.cells.feature, 
                                                    latent.vars = latent.vars, test.use = test.use, verbose = verbose), 
                       MAST = MASTDETest(data.use = data.use, cells.1 = cells.1, 
                                         cells.2 = cells.2, latent.vars = latent.vars, verbose = verbose, 
                                         ...), 
                       DESeq2 = Seurat:::DESeq2DETest(data.use = data.use, 
                                                      cells.1 = cells.1, cells.2 = cells.2, verbose = verbose, 
                                                      ...), 
                       LR = Seurat:::LRDETest(data.use = data.use, cells.1 = cells.1, 
                                              cells.2 = cells.2, latent.vars = latent.vars, verbose = verbose), 
                       stop("Unknown test: ", test.use))
  return(de.results)
}

# Differential expression testing using  Student's t-test
#
# Identify differentially expressed genes between two groups of cells using
#  Student's t-tests
#
# @return Returns a p-value ranked matrix of putative differentially expressed
# genes.
#
#' @importFrom pbapply pbsapply
#' @importFrom future.apply future_sapply
#' @importFrom future nbrOfWorkers
#' @export
#
DiffTTest <- function(data.use, cells.1, cells.2, verbose = TRUE) 
{
  my.sapply <- ifelse(test = verbose && future:::nbrOfWorkers() == 1, 
                      yes = pbapply::pbsapply, no = future.apply::future_sapply)
  p_val <- unlist(x = my.sapply(X = 1:nrow(data.use), FUN = function(x) {
    t.test(x = data.use[x, cells.1], y = data.use[x, cells.2])$p.value
  }))
  to.return <- data.frame(p_val = p_val, pval_rep =p_val, row.names = rownames(x = data.use))
  return(to.return)
}


#' # Differential expression testing using Wilcoxon Rank Sum test
#' #
#' # Identify differentially expressed genes between two groups of cells using
#' # Wilcoxon Rank Sum test
#' #
#' # @return Returns a p-value ranked matrix of putative differentially expressed
#' # genes.
#' #
#' #' @importFrom pbapply pbsapply
#' #' @importFrom future.apply future_sapply
#' #' @importFrom future nbrOfWorkers
#' #' @export
#' #' 
#' WilcoxDETest <- function (data.use, cells.1, cells.2, verbose = TRUE, limma = FALSE, 
#'           ...) 
#' {
#'   data.use <- data.use[, c(cells.1, cells.2), drop = FALSE]
#'   j <- seq_len(length.out = length(x = cells.1))
#'   my.sapply <- ifelse(test = verbose && nbrOfWorkers() == 1, 
#'                       yes = pbsapply, no = future_sapply)
#'   overflow.check <- ifelse(test = is.na(x = suppressWarnings(length(x = data.use[1, 
#'   ]) * length(x = data.use[1, ]))), yes = FALSE, no = TRUE)
#'   presto.check <- PackageCheck("presto", error = FALSE)
#'   limma.check <- PackageCheck("limma", error = FALSE)
#'   group.info <- data.frame(row.names = c(cells.1, cells.2))
#'   group.info[cells.1, "group"] <- "Group1"
#'   group.info[cells.2, "group"] <- "Group2"
#'   group.info[, "group"] <- factor(x = group.info[, "group"])
#'   if (presto.check[1] && (!limma)) {
#'     data.use <- data.use[, rownames(group.info), drop = FALSE]
#'     res <- presto::wilcoxauc(X = data.use, y = group.info[, 
#'                                                           "group"])
#'     res <- res[1:(nrow(x = res)/2), ]
#'     p_val <- res$pval
#'   }
#'   else {
#'     if (getOption("Seurat.presto.wilcox.msg", TRUE) && (!limma)) {
#'       message("For a (much!) faster implementation of the Wilcoxon Rank Sum Test,", 
#'               "\n(default method for FindMarkers) please install the presto package", 
#'               "\n--------------------------------------------", 
#'               "\ninstall.packages('devtools')", "\ndevtools::install_github('immunogenomics/presto')", 
#'               "\n--------------------------------------------", 
#'               "\nAfter installation of presto, Seurat will automatically use the more ", 
#'               "\nefficient implementation (no further action necessary).", 
#'               "\nThis message will be shown once per session")
#'       options(Seurat.presto.wilcox.msg = FALSE)
#'     }
#'     if (limma.check[1] && overflow.check) {
#'       p_val <- my.sapply(X = 1:nrow(x = data.use), FUN = function(x) {
#'         return(min(2 * min(limma::rankSumTestWithCorrelation(index = j, 
#'                                                              statistics = data.use[x, ])), 1))
#'       })
#'     }
#'     else {
#'       if (limma && overflow.check) {
#'         stop("To use the limma implementation of the Wilcoxon Rank Sum Test,\n        please install the limma package:\n        --------------------------------------------\n        install.packages('BiocManager')\n        BiocManager::install('limma')\n        --------------------------------------------")
#'       }
#'       else {
#'         data.use <- data.use[, rownames(x = group.info), 
#'                              drop = FALSE]
#'         p_val <- my.sapply(X = 1:nrow(x = data.use), 
#'                            FUN = function(x) {
#'                              return(wilcox.test(data.use[x, ] ~ group.info[, 
#'                                                                            "group"], ...)$p.value)
#'                            })
#'       }
#'     }
#'   }
#'   return(data.frame(p_val, row.names = rownames(x = data.use)))
#' }



# Differential expression testing using weigted Student's t-test
#
# Identify differentially expressed genes between two groups of cells using
# weighted Student's t-tests
#
# @return Returns a p-value ranked matrix of putative differentially expressed
# genes.
#
#' @importFrom pbapply pbsapply
#' @importFrom future.apply future_sapply
#' @importFrom future nbrOfWorkers
#' @export

DiffWeightedTTest <- function(data.use, 
                              cells.1, 
                              cells.2,
                              weights.1,
                              weights.2, 
                              use.sapply = FALSE,
                              return.all = TRUE,
                              do.bootstrapping = FALSE,
                              verbose = TRUE) 
{
  
  
  if (do.bootstrapping) {
    bootse <- TRUE
    bootp <- TRUE
  }
  else {
    bootse <- FALSE
    bootp <- FALSE
  }
  
  
  if (use.sapply) {
    my.sapply <- ifelse(test = verbose && future:::nbrOfWorkers() == 1, 
                        yes = pbapply::pbsapply, no = future.apply::future_sapply)
    
    p_val <- unlist(x = my.sapply(X = 1:nrow(data.use), FUN = function(x) {
      res <- weights::wtd.t.test(x = data.use[x,cells.1], y = data.use[x,cells.2],
                                 weight = weights.1, weighty = weights.2,
                                 mean1 = FALSE, samedata = FALSE, bootse = bootse,
                                 bootp = bootp)
      return(res$coefficients[['p.value']])
    }))
    
    to.return <- data.frame(p_val = p_val, row.names = rownames(x = data.use))
    
  } else {
    w.t.test <- apply(data.use, 1, function(x) {
      weights::wtd.t.test(x = x[cells.1], y = x[cells.2],
                          weight = weights.1, weighty = weights.2,
                          mean1 = TRUE, samedata = FALSE, bootse = bootse,
                          bootp = bootp)
    })
    if (length((w.t.test)) > 1) {
      to.return <- data.frame(matrix(unlist(w.t.test), ncol = 8, 
                                     byrow = T))
      colnames(to.return) <-c("test", "t.value", "df", "p.value", 
                              "difference", "mean_x", "mean_y", "std_err")
      rownames(to.return) <- names(w.t.test)
      
      
      # t.test(x = data.use[x, cells.1], y = data.use[x, cells.2])$p.value
    }
    if (return.all) { 
      to.return   <-data.frame(p_val = as.numeric(to.return$p.value), 
                               t_value = as.numeric(to.return$t.value),
                               df = as.numeric(to.return$df),
                               row.names = rownames(to.return))
      
    } else {
      to.return <- data.frame(p_val = as.numeric(to.return$p.value), row.names = rownames(to.return))
    }
  }
  
  return(to.return)
}

#' @inheritParams Seurat::FindMarkers
#' @rdname FindMarkers.SuperCell
#' @concept differential_expression
#' @export
#' @method FindMarkers.SuperCell Assay


FindMarkers.SuperCell.Assay <- function(object, slot = "data", cells.1 = NULL, cells.2 = NULL, 
                                        weights.1 = NULL, weights.2 = NULL,
                                        features = NULL, test.use = "survey_weighted_t", fc.slot = "data", pseudocount.use = 1, 
                                        norm.method = NULL, mean.fxn = NULL, fc.name = NULL, base = 2, 
                                        ...) 
{
  data.slot <- ifelse(test = test.use %in% Seurat:::DEmethods_counts(), 
                      yes = "counts", no = slot)
  if (length(x = Layers(object = object, search = slot)) > 
      1) {
    stop(slot, " layers are not joined. Please run JoinLayers")
  }
  data.use <- GetAssayData(object = object, slot = data.slot)
  fc.results <- FoldChange.SuperCell(object = object, slot = fc.slot, 
                                     cells.1 = cells.1, cells.2 = cells.2, 
                                     weights.1 = weights.1, weights.2 = weights.2,
                                     features = features, 
                                     pseudocount.use = pseudocount.use, mean.fxn = mean.fxn, 
                                     fc.name = fc.name, base = base, norm.method = norm.method)
  de.results <- FindMarkers.SuperCell(object = data.use, cells.1 = cells.1, 
                                      weights.1 = weights.1, weights.2 = weights.2,
                                      cells.2 = cells.2, features = features, test.use = test.use, 
                                      fc.results = fc.results, ...)
  return(de.results)
}

#' @inheritParams Seurat::FindMarkers
#' @rdname FindMarkers.SuperCell
#' @concept differential_expression
#' @export
#' @method FindMarkers.SuperCell Seurat

FindMarkers.SuperCell.Seurat <- function (object, ident.1 = NULL, ident.2 = NULL, latent.vars = NULL, 
                                          group.by = NULL, subset.ident = NULL, assay = NULL, reduction = NULL, 
                                          ...) 
{
  if (!is.null(x = group.by)) {
    if (!is.null(x = subset.ident)) {
      object <- subset(x = object, idents = subset.ident)
    }
    Idents(object = object) <- group.by
  }
  if (!is.null(x = assay) && !is.null(x = reduction)) {
    stop("Please only specify either assay or reduction.")
  }
  if (length(x = ident.1) == 0) {
    stop("At least 1 ident must be specified in `ident.1`")
  }
  if (is.null(x = reduction)) {
    assay <- assay %||% DefaultAssay(object = object)
    data.use <- object[[assay]]
    cellnames.use <- colnames(x = data.use)
  }
  else {
    data.use <- object[[reduction]]
    cellnames.use <- rownames(x = data.use)
  }
  cells <- Seurat:::IdentsToCells(object = object, ident.1 = ident.1, 
                                  ident.2 = ident.2, cellnames.use = cellnames.use)
  cells <- sapply(X = cells, FUN = intersect, y = cellnames.use, 
                  simplify = FALSE, USE.NAMES = TRUE)
  
  weights.1 <- Seurat::FetchData(object,vars = "size", cells = cells$cells.1)$size
  weights.2 <- Seurat::FetchData(object,vars = "size", cells = cells$cells.2)$size
  
  
  if (!all(vapply(X = cells, FUN = length, FUN.VALUE = integer(length = 1L)))) {
    abort(message = "Cells in one or both identity groups are not present in the data requested")
  }
  if (!is.null(x = latent.vars)) {
    latent.vars <- Seurat::FetchData(object = object, vars = latent.vars, 
                                     cells = c(cells$cells.1, cells$cells.2))
  }
  norm.command <- paste0("NormalizeData.", assay)
  norm.method <- if (norm.command %in% Command(object = object) && 
                     is.null(x = reduction)) {
    Command(object = object, command = norm.command, value = "normalization.method")
  }
  else if (length(x = intersect(x = c("FindIntegrationAnchors", 
                                      "FindTransferAnchors"), y = Command(object = object)))) {
    command <- intersect(x = c("FindIntegrationAnchors", 
                               "FindTransferAnchors"), y = Command(object = object))[1]
    Command(object = object, command = command, value = "normalization.method")
  }
  else {
    NULL
  }
  de.results <- FindMarkers.SuperCell(object = data.use, latent.vars = latent.vars, 
                                      cells.1 = cells$cells.1, cells.2 = cells$cells.2, 
                                      weights.1 = weights.1, weights.2 = weights.2,
                                      norm.method = norm.method, 
                                      ...)
  return(de.results)
}

#' Fold Change
#'
#' Calculate log fold change and percentage of cells expressing each feature
#' for different identity classes taking into account metacell sizes.
#'
#' If the slot is \code{scale.data} or a reduction is specified, average difference
#' is returned instead of log fold change and the column is named "avg_diff".
#' Otherwise, log2 fold change is returned with column named "avg_log2_FC".
#'
#' @examples
#' \dontrun{
#' data("pbmc_small")
#' FoldChange.SuperCell(pbmc_small, ident.1 = 1)
#' }
#' 
#' @param object A Seurat metacell object with a size metadata column.
#' @param ... Arguments passed to other methods
#' @rdname FoldChange.SuperCell
#' @export FoldChange.SuperCell
#' @return Returns a data.frame
#' @seealso \code{FindMarkers.SuperCell}
FoldChange.SuperCell <- function(object, ...) 
{
  UseMethod(generic = "FoldChange.SuperCell", object = object)
}

#' Weighted mean 
#' @rdname weighted.mean.fxn
#' @param x vector of observation (eg. metacell counts for a gene)*
#' @param weights vector of weights, same length as x (eg. metacell sizes)
#' @export
weighted.mean.fxn <-  function(x, weights) {Matrix::rowSums(x %*% Matrix::Diagonal(x = weights)) / sum(weights)}


#' @rdname FoldChange.SuperCell
#' @concept differential_expression
#' @export
#' @method FoldChange.SuperCell Assay

FoldChange.SuperCell.Assay <- function(object, 
                                       cells.1, cells.2,
                                       weights.1, weights.2,
                                       features = NULL, slot = "data", 
                                       pseudocount.use = 1, fc.name = NULL, mean.fxn = NULL, base = 2, 
                                       norm.method = NULL, ...) 
{
  data <- GetAssayData(object = object, slot = slot)
  
  
  log1pdata.mean.fxn <- function(x, weights) {
    return(log1p(x = (((Matrix::rowSums(expm1(x) %*% Matrix::Diagonal(x = weights))+pseudocount.use) / sum(weights)) -1))/log(base))
  }
  
  # log1pdata.mean.fxn <- function(x) {
  #   return(log(x = (rowSums(x = expm1(x = x)) + pseudocount.use)/NCOL(x),
  #              base = base))
  # }
  
  scaledata.mean.fxn <-  function(x, weights) {Matrix::rowSums(x %*% Matrix::Diagonal(x = weights)) / sum(weights)}
  # scaledata.mean.fxn <- rowMeans
  
  counts.mean.fxn <- function(x, weights) {
    return(log1p(x = (((Matrix::rowSums(x %*% Matrix::Diagonal(x = weights))+pseudocount.use) / sum(weights)) -1))/log(base))
  }
  
  # counts.mean.fxn <- function(x) {
  #   return(log(x = (rowSums(x = x) + pseudocount.use)/NCOL(x), 
  #              base = base))
  # }
  
  
  if (!is.null(x = norm.method)) {
    if (norm.method != "LogNormalize") {
      new.mean.fxn <- counts.mean.fxn
    }
    else {
      new.mean.fxn <- counts.mean.fxn
      if (slot == "data") {
        new.mean.fxn <- log1pdata.mean.fxn
      }
      else if (slot == "scale.data") {
        new.mean.fxn <- scaledata.mean.fxn
      }
    }
  }
  else {
    new.mean.fxn <- switch(EXPR = slot, data = log1pdata.mean.fxn, 
                           scale.data = scaledata.mean.fxn, counts = counts.mean.fxn, 
                           log1pdata.mean.fxn)
  }
  mean.fxn <- mean.fxn %||% new.mean.fxn
  base.text <- ifelse(test = base == exp(1), yes = "", no = base)
  fc.name <- fc.name %||% ifelse(test = slot == "scale.data", 
                                 yes = "avg_diff", no = paste0("avg_log", base.text, "FC"))
  FoldChange.SuperCell(object = data, cells.1 = cells.1, cells.2 = cells.2, 
                       weights.1 = weights.1, weights.2 = weights.2,
                       features = features, mean.fxn = mean.fxn, fc.name = fc.name)
}


#' @rdname FoldChange.SuperCell
#' @inheritParams Seurat::FoldChange
#' @concept differential_expression
#' @export
#' @method FoldChange.SuperCell default

FoldChange.SuperCell.default <- function(object, 
                                         cells.1, cells.2,
                                         weights.1, weights.2,
                                         mean.fxn, fc.name, features = NULL, 
                                         ...) 
{
  
  features <- features %||% rownames(x = object)
  thresh.min <- 0
  pct.1 <- rowSums((object[features, cells.1] > thresh.min) %*% Matrix::Diagonal(x = weights.1)) / sum(weights.1)
  pct.2 <- rowSums((object[features, cells.2]  > thresh.min) %*% Matrix::Diagonal(x = weights.2)) / sum(weights.2)
  
  data.1 <- mean.fxn(object[features, cells.1,drop = FALSE],weights.1)
  data.2 <- mean.fxn(object[features, cells.2,drop = FALSE],weights.2)
  fc <- (data.1 - data.2)
  fc.results <- as.data.frame(x = cbind(fc, pct.1, pct.2))
  colnames(fc.results) <- c(fc.name, "pct.1", "pct.2")
  return(fc.results)
}


#' @inheritParams Seurat::FoldChange.Seurat
#' @param object Metacell seurat object with a size meta.Data column
#' @rdname FoldChange.SuperCell
#' @concept differential_expression
#' @export
#' @method FoldChange.SuperCell Seurat
FoldChange.SuperCell.Seurat <- function(object, ident.1 = NULL, ident.2 = NULL, group.by = NULL, 
                                        subset.ident = NULL, assay = NULL, slot = "data", reduction = NULL, 
                                        features = NULL, pseudocount.use = 1, mean.fxn = NULL, base = 2, 
                                        fc.name = NULL, ...) 
{
  if (!is.null(x = group.by)) {
    if (!is.null(x = subset.ident)) {
      object <- subset(x = object, idents = subset.ident)
    }
    Idents(object = object) <- group.by
  }
  if (!is.null(x = assay) && !is.null(x = reduction)) {
    stop("Please only specify either assay or reduction.")
  }
  if (is.null(x = reduction)) {
    assay <- assay %||% DefaultAssay(object = object)
    data.use <- object[[assay]]
    cellnames.use <- colnames(x = data.use)
  }
  else {
    data.use <- object[[reduction]]
    cellnames.use <- rownames(data.use)
  }
  cells <- IdentsToCells(object = object, ident.1 = ident.1, 
                         ident.2 = ident.2, cellnames.use = cellnames.use)
  
  weights.1 <- Seurat::FetchData(object,vars = "size", cells = cells$cells.1)$size
  weights.2 <- Seurat::FetchData(object,vars = "size", cells = cells$cells.2)$size
  
  norm.command <- paste0("NormalizeData.", assay)
  norm.method <- if (norm.command %in% Command(object = object) && 
                     is.null(x = reduction)) {
    Command(object = object, command = norm.command, value = "normalization.method")
  }
  else if (length(x = intersect(x = c("FindIntegrationAnchors", 
                                      "FindTransferAnchors"), y = Command(object = object)))) {
    command <- intersect(x = c("FindIntegrationAnchors", 
                               "FindTransferAnchors"), y = Command(object = object))[1]
    Command(object = object, command = command, value = "normalization.method")
  }
  else {
    NULL
  }
  fc.results <- FoldChange(object = data.use, cells.1 = cells$cells.1, 
                           cells.2 = cells$cells.2, weights.1 = weights.1, weights.2 = weights.2,
                           features = features, slot = slot, 
                           pseudocount.use = pseudocount.use, mean.fxn = mean.fxn, 
                           base = base, fc.name = fc.name, norm.method = norm.method)
  return(fc.results)
}

#' Gene expression markers for all identity classes
#'
#' Same as Seurat FindAllMarkers function taking into account metacell size.
#' Only weighted t test using weights::wtd.t.test funnction is currently supported
#' 
#'
#' @param object Seurat metacell object with a metadata size column
#' @param test.use onky "weighted_t" is supported.
#' @inheritParams Seurat::FindAllMarkers
#' @return data.frame with a ranked list of putative markers as rows, and associated statistics as columns (p-values, t statistic, degree of freedom (df) of the weighted t-test). 
#' avg_logFC, pct.1, pc.2 are computed using Seurat v5 approach (pseudocounts added at the group level) taking into account metacell sizes.
#' @import Seurat
#' @rdname FindAllMarkers.SuperCell
#' @export

FindAllMarkers.SuperCell <- function(object, assay = NULL, features = NULL, logfc.threshold = 0.1, 
                                     test.use = "survey_weighted_t", slot = "data", min.pct = 0.01, min.diff.pct = -Inf, 
                                     node = NULL, verbose = TRUE, only.pos = FALSE, max.cells.per.ident = Inf, 
                                     random.seed = 1, latent.vars = NULL, min.cells.feature = 3, 
                                     min.cells.group = 3, mean.fxn = NULL, fc.name = NULL, base = 2, 
                                     return.thresh = 0.01, densify = FALSE, ...) 
{
  MapVals <- function(vec, from, to) {
    vec2 <- setNames(object = to, nm = from)[as.character(x = vec)]
    vec2[is.na(x = vec2)] <- vec[is.na(x = vec2)]
    return(unname(obj = vec2))
  }
  if ((test.use == "roc") && (return.thresh == 0.01)) {
    return.thresh <- 0.7
  }
  if (is.null(x = node)) {
    idents.all <- sort(x = unique(x = Idents(object = object)))
  }
  else {
    if (!PackageCheck("ape", error = FALSE)) {
      stop(cluster.ape, call. = FALSE)
    }
    tree <- Tool(object = object, slot = "BuildClusterTree")
    if (is.null(x = tree)) {
      stop("Please run 'BuildClusterTree' before finding markers on nodes")
    }
    descendants <- DFT(tree = tree, node = node, include.children = TRUE)
    all.children <- sort(x = tree$edge[, 2][!tree$edge[, 
                                                       2] %in% tree$edge[, 1]])
    descendants <- MapVals(vec = descendants, from = all.children, 
                           to = tree$tip.label)
    drop.children <- setdiff(x = tree$tip.label, y = descendants)
    keep.children <- setdiff(x = tree$tip.label, y = drop.children)
    orig.nodes <- c(node, as.numeric(x = setdiff(x = descendants, 
                                                 y = keep.children)))
    tree <- ape::drop.tip(phy = tree, tip = drop.children)
    new.nodes <- unique(x = tree$edge[, 1, drop = TRUE])
    idents.all <- (tree$Nnode + 2):max(tree$edge)
  }
  genes.de <- list()
  messages <- list()
  for (i in 1:length(x = idents.all)) {
    if (verbose) {
      message("Calculating cluster ", idents.all[i])
    }
    genes.de[[i]] <- tryCatch(expr = {
      FindMarkers.SuperCell(object = object, assay = assay, ident.1 = if (is.null(x = node)) {
        idents.all[i]
      }
      else {
        tree
      }, ident.2 = if (is.null(x = node)) {
        NULL
      }
      else {
        idents.all[i]
      }, features = features, logfc.threshold = logfc.threshold, 
      test.use = test.use, slot = slot, min.pct = min.pct, 
      min.diff.pct = min.diff.pct, verbose = verbose, 
      only.pos = only.pos, max.cells.per.ident = max.cells.per.ident, 
      random.seed = random.seed, latent.vars = latent.vars, 
      min.cells.feature = min.cells.feature, min.cells.group = min.cells.group, 
      mean.fxn = mean.fxn, fc.name = fc.name, base = base, 
      densify = densify, ...)
    }, error = function(cond) {
      return(cond$message)
    })
    if (is.character(x = genes.de[[i]])) {
      messages[[i]] <- genes.de[[i]]
      genes.de[[i]] <- NULL
    }
  }
  gde.all <- data.frame()
  for (i in 1:length(x = idents.all)) {
    if (is.null(x = unlist(x = genes.de[i]))) {
      next
    }
    gde <- genes.de[[i]]
    if (nrow(x = gde) > 0) {
      if (test.use == "roc") {
        gde <- subset(x = gde, subset = (myAUC > return.thresh | 
                                           myAUC < (1 - return.thresh)))
      }
      else if (is.null(x = node) || test.use %in% c("bimod", 
                                                    "t")) {
        gde <- gde[order(gde$p_val, -abs(gde$pct.1 - 
                                           gde$pct.2)), ]
        gde <- subset(x = gde, subset = p_val < return.thresh)
      }
      if (nrow(x = gde) > 0) {
        gde$cluster <- idents.all[i]
        gde$gene <- rownames(x = gde)
      }
      if (nrow(x = gde) > 0) {
        gde.all <- rbind(gde.all, gde)
      }
    }
  }
  if ((only.pos) && nrow(x = gde.all) > 0) {
    return(subset(x = gde.all, subset = gde.all[, 2] > 0))
  }
  rownames(x = gde.all) <- make.unique(names = as.character(x = gde.all$gene))
  if (nrow(x = gde.all) == 0) {
    warning("No DE genes identified", call. = FALSE, immediate. = TRUE)
  }
  if (length(x = messages) > 0) {
    warning("The following tests were not performed: ", call. = FALSE, 
            immediate. = TRUE)
    for (i in 1:length(x = messages)) {
      if (!is.null(x = messages[[i]])) {
        warning("When testing ", idents.all[i], " versus all:\n\t", 
                messages[[i]], call. = FALSE, immediate. = TRUE)
      }
    }
  }
  if (!is.null(x = node)) {
    gde.all$cluster <- MapVals(vec = gde.all$cluster, from = new.nodes, 
                               to = orig.nodes)
  }
  return(gde.all)
}




#' @method FindMarkers.SuperCell StdAssay
#' @export
#'
FindMarkers.SuperCell.StdAssay <- FindMarkers.SuperCell.Assay

#' @method FoldChange.SuperCell StdAssay
#' @export
#'
FoldChange.SuperCell.StdAssay <- FoldChange.SuperCell.Assay