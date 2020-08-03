#' Preprocess gene expression count matrix
#'
#' This preprocessing function calculates library size, library size 
#' normalization factor, highly variable genes, normalized gene expression 
#' level. These results can be used in training SAVER-CVAE model.
#'
#' This function normalized gene expression count matrix by size factor, then 
#' log scale and standardize it. Size factor is library size diveded by mean 
#' library size. It calculates highly variable genes per batch using 
#' \code{Seurat::FindVariableFeatures} function and takes the union. 
#' 
#' @param x Expression count matrix. The rows correspond to genes and colomns 
#' correspond to cells. 
#'
#' @param B Observed covariates matrix. The rows correspond to cells. 
#' Default is NULL, then scaled log library size will be taken as B.
#'
#' @param n_hvg Number of highly variable genes to select. Default is 3000.
#' 
#' @param hvg_genes_only Whether to subset and normalize expression count matrix
#' using only highly variable genes. TRUE to use only highly variable genes. 
#' FALSE to use genes with non-zero expression in at least perc.exp*100 % of 
#' cells. Default is TRUE.
#' 
#' @param perc.exp Good genes are defined as those with a mean expression of 
#' more than perc.exp. We select highly variable genes from good 
#' genes.
#' 
#' 
#' @return A list with the following components
#' \item{\code{x}}{A subset of expression count matrix. The rows correspond to 
#' cells and colomns correspond to highly variable genes or genes with a mean 
#' expression of more than perc.exp.}
#' \item{\code{x.norm}}{Log normailized and standardized gene expression matrix 
#' (cells by genes).}
#' \item{\code{sf}}{Cell specific size factors used.}
#' \item{\code{B}}{Observed covariates matrix or scaled log library size.}
#' \item{\code{hvg}}{Names of highly variable genes.}
#' 
#' @import tensorflow
#' @import keras
#' @import dplyr
#' @import irlba
#' @import uwot
#' @import cowplot
#' @import ggplot2
#' @import Rtsne
#'
#' @export


preprocess_data <- function(x, B = NULL, n_hvg = 3000, hvg_genes_only = TRUE,
                            perc.exp = 0.01) {
  
  norm_x = FALSE
  hvg_out = TRUE
  
  if (tensorflow::tf$executing_eagerly())
    tensorflow::tf$compat$v1$disable_eager_execution()
  
  libsize <- colSums(x)
  # begin change
  #if (any(libsize == 0)) {
  #  x <- x[, libsize != 0]
  #  libsize <- libsize[libsize != 0]
  #}
  # end change
  n <- ncol(x)
  if (is.null(B)) {
    B <- scale(matrix(log10(libsize)))
  }
  hvg <- NULL
  set.seed(5)
  ind <- sort(sample(1:n, min(n, 20000)))
  good.genes <- which(rowMeans(x[, ind] != 0) >= perc.exp)
  if (norm_x) {
    sf <- libsize/exp(mean(log(libsize)))
  } else {
    sf <- libsize/mean(libsize)
  }
  if (hvg_genes_only & n_hvg < length(good.genes)) {
    cat("Training on", n_hvg, "highly variable genes.\n")
    x.filt <- x[good.genes, ind]
    unique_elements <- apply(B, 2, function(x) length(unique(x)))
    batches <- which(unique_elements == 2)
    if (length(batches) == 0) {
      hvg.ind <- seurat_v3_hvg(x.filt, nfeatures = n_hvg)
    } else if (length(batches) == 1) {
      ind1 <- B[ind, batches] == max(B[ind, batches])
      hvg.ind1 <- seurat_v3_hvg(x.filt[, ind1], nfeatures = n_hvg)
      hvg.ind2 <- seurat_v3_hvg(x.filt[, !ind1], nfeatures = n_hvg)
      hvg.ind <- unique(c(rbind(hvg.ind1, hvg.ind2)))[1:n_hvg]
    } else {
      hvg_list <- vector("list", length(batches))
      for (i in seq_along(batches)) {
        # begin change
        batch_ind <- which(B[ind, batches[i]] == max(B[ind, batches[i]]))
        # end change
        hvg_list[[i]] <- seurat_v3_hvg(x.filt[, batch_ind], nfeatures = n_hvg)
      }
      hvg.ind <- unique(c(do.call("rbind", hvg_list)))[1:n_hvg]
    }
    x.sub <- x[good.genes[hvg.ind], ]
    hvg <- rownames(x.filt)[hvg.ind]
    if (hvg_out) {
      x <- x.sub
    }
  } else {
    cat("Selecting genes with non-zero expression in", perc.exp*100, "% of cells.\n")
    x.sub <- x[good.genes, ]
    if (hvg_out) {
      x <- x[good.genes, ]
    }
  }
  x <- t(x)
  x.sub <- t(x.sub)
  if (class(x) != "matrix") {
    sweep_fun <- sweep_sparse
  } else {
    sweep_fun <- sweep
  }
  x.sub <- sweep_fun(x.sub, 1, sf, "/")
  if (norm_x) {
    x <- sweep_fun(x, 1, sf, "/")
  }
  if (class(x) != "matrix") {
    x <- as.matrix(x)
    gc(verbose = FALSE)
    x.sub@x <- log(x.sub@x+1)
    x.sub <- as.matrix(x.sub)
    gc(verbose = FALSE)
    x.sub <- scale(x.sub)
  } else {
    x.sub <- scale(log(x.sub+1))
  }
  bad.ind <- x.sub > 50
  if (any(bad.ind)) {
    x.sub[bad.ind] <- 50
  }
  gc(verbose = FALSE)
  return(list(x = x, x.norm = x.sub, sf = sf, B = B, hvg = hvg))
}
