#' Calculate SAVER-CVAE corrected counts, SAVER estimates and sampled SAVER
#' estimates.
#'
#' Given trained encoder and decoder model as input, calculate SAVER-CVAE
#' corrected counts, SAVER estimates and sampled SAVER estimates.
#'
#' @param decoder Decoder model.
#'
#' @param x Expression count matrix. The rows correspond to genes and colomns
#' correspond to cells.
#'
#' @param z Mean of posterior normal distribution, which is the low-dimensional
#' representation of the cells.
#'
#' @param B Observed covariates matrix used to train the model.
#' The rows correspond to cells.
#'
#' @param bad.genes Index of genes that could not be predicted accurately
#' by the autoencoder based on cross validation results. Instead, for these
#' genes, denoised SAVER estimate are mean gene expression level across cells.
#'
#' @param subset.genes Index of a subset of genes to be predicted.
#' Default is NULL, then all genes will be predicted.
#'
#' @return A list with the following components
#' \item{\code{estimate}}{Denoised SAVER estimates based on empirical Bayes
#' method. The SAVER estimate and posterior distribution can be used to
#' visualize gene expression and perform trajectory analyses. The SAVER
#' estimate acts as the denoised estimate of the true expression}
#' \item{\code{se}}{Standard error of SAVER posterior distribution.}
#' \item{\code{samp}}{Sampling from the SAVER posterior distribution.
#' Sampling from the SAVER posterior allows investigation on gene expression
#' distribution.}
#' \item{\code{adj.count}}{SAVER-CVAE corrected counts using cdf matching. This
#' preserves the technical noise characteristics of scRNA-seq expression count
#' data after adjusting for the observed covariates. This is useful in
#' downstream scenarios where hypothesis testing on the original count data
#' is desired, such as in differential expression.}
#' \item{\code{time}}{Runtime.}
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



generate_saver <- function(decoder, x, z, B, bad.genes, subset.genes = NULL) {

  if (tensorflow::tf$executing_eagerly())
    tensorflow::tf$compat$v1$disable_eager_execution()

  row.genes = TRUE

  pt <- Sys.time()
  set.seed(5)
  if (row.genes) {
    out_dim <- nrow(x)
    libsize <- colSums(x)
  } else {
    out_dim <- ncol(x)
    libsize <- rowSums(x)
  }
  n <- nrow(z)
  B_blank <- matrix(0, n, ncol(B))
  batches <- split(1:n, ceiling(seq_along(1:n)/1000))
  sf <- libsize/mean(libsize)

  XZ.sum <- rep(0, out_dim)
  XZ.theta.sum <- rep(0, out_dim)
  for (batch in batches) {
    XZ.out <- decoder %>% predict(list(z[batch, , drop = FALSE],
                                       B_blank[batch, , drop = FALSE]))
    XZ.sum <- XZ.sum + colSums(XZ.out[[1]])
    XZ.theta.sum <- XZ.theta.sum + colSums(XZ.out[[2]])
  }

  XZ.mean <- XZ.sum/n
  XZ.theta.mean <- XZ.theta.sum/n

  cat("Calculating SAVER output\n")
  adj.count <- matrix(0, out_dim, n)
  saver.est <- matrix(0, out_dim, n)
  saver.se <- matrix(0, out_dim, n)
  saver.samp <- matrix(0, out_dim, n)

  if (row.genes) {
    mat.names <- dimnames(x)
  } else {
    mat.names <- dimnames(x)[2:1]
  }

  dimnames(adj.count) <- mat.names
  dimnames(saver.est) <- mat.names
  dimnames(saver.se) <- mat.names
  dimnames(saver.samp) <- mat.names
  msg.chunk <- round(quantile(1:length(batches), 0:10/10))
  batch_ind <- 1
  for (batch in batches) {
    if (batch_ind %in% msg.chunk) {
      cat(names(msg.chunk[batch_ind == msg.chunk]), "...")
    }
    X.out <- predict_decoder(decoder, z[batch, , drop = FALSE],
                             B[batch, , drop = FALSE])
    XZ.out <- predict_decoder(decoder, z[batch, , drop = FALSE],
                              B_blank[batch, , drop = FALSE])
    if (row.genes) {
      x.batch <- x[, batch]
    } else {
      x.batch <- t(x[batch, ])
    }
    X.out.sf <- X.out
    X.out.sf[[1]] <- sweep(X.out[[1]], 2, sf[batch], "*")
    xz <- generate_newx(x.batch, X.out.sf, XZ.out)
    XZ.out[[1]][bad.genes, ] <- XZ.mean[bad.genes]
    XZ.out[[2]][bad.genes, ] <- XZ.theta.mean[bad.genes]
    alpha <- xz + XZ.out[[2]]
    beta <- 1 + XZ.out[[2]]/XZ.out[[1]]
    estimate <- alpha/beta
    se <- sqrt(alpha/beta^2)
    samp <- sapply(1:length(batch), function(j) rgamma(out_dim, alpha[, j], beta[, j]))
    ex <- 10^5
    estimate <- floor(estimate*ex)/ex
    se <- ceiling(se*ex)/ex
    samp <- round(samp, 5)
    adj.count[, batch] <- round(xz, 5)
    saver.est[, batch] <- estimate
    saver.se[, batch] <- se
    saver.samp[, batch] <- samp
    batch_ind <- batch_ind + 1
  }
  cat("\n")
  gc(verbose = FALSE)
  if (!is.null(subset.genes)) {
    saver.est = saver.est[subset.genes,]
    saver.se = saver.se[subset.genes,]
    saver.samp = saver.samp[subset.genes,]
    adj.count = adj.count[subset.genes,]
  }
  list(estimate = saver.est, se = saver.se, samp = saver.samp,
       adj.count = adj.count, time = Sys.time() - pt)
}
