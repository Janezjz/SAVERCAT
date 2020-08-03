#' Train a decoder which generates denoised expression matrix 
#'
#' This is the second step of training SAVER-CVAE. 
#' The decoder takes variational posterior distribution parameters estimated in 
#' the first step as well as observed covariates or library size as input,
#' and produces the covariate-adjusted parameters of the negative binomial 
#' distribution. 
#'
#' This step starts by sampling from posterior normal distribution of the 
#' latent variable. Then concatenate with observed covariates or scaled library
#' size as input to the decoder network. The outputs of the decoder are 
#' parameters of the negative binomial distribution.
#' Cross-validation is used to select the epoch which gives the lowest 
#' validation error. 
#' 
#' @param x Expression count matrix. The rows correspond to genes and colomns 
#' correspond to cells. 
#' 
#' @param z.mean Mean of posterior normal distribution of the latent variable,
#' which is the low-dimensional representation of the cells.
#' 
#' @param z.logvar Log of variance of posterior normal distribution of the
#' latent variable.
#'
#' @param B Observed covariates matrix. The rows correspond to cells. 
#' Default is NULL, then scaled log library size will be taken as B.
#'
#' @param sf Vector of cell specific size factors. Default is NULL, then scaled
#' library size will be used.
#' 
#' @param epochs Number of epochs to train the model. Default is 50.
#' 
#' @param alpha Weight of reconstruction loss using only observed covariates B.
#' Default is 0.01.
#' 
#' @param batch_size Number of samples per gradient update. Default is 64.
#' 
#' @param verbose Verbosity mode. 0 is silent, 1 is progress bar, 2 is one line
#' per epoch. Default is 2.
#' 
#' @param save.model Folder specified to save the model.
#' 
#' @param return_saver Whether to calculate corrected expression counts and 
#' SAVER estimates. Default is FALSE.
#' 
#' @return A list with the following components
#' \item{\code{z}}{Mean of posterior normal distribution, which is
#' the low-dimensional representation of the cells.}
#' \item{\code{z.var}}{Log of variance of posterior normal distribution.}
#' \item{\code{z.samp}}{Random sampling from posterior normal distribution of 
#' the latent characterization.}
#' \item{\code{B}}{Covariates matrix used to train the model.}
#' \item{\code{sf}}{Cell specific size factors used.}
#' \item{\code{adj.raw}}{SAVER-CVAE corrected counts using cdf matching. This 
#' preserves the technical noise characteristics of scRNA-seq expression count 
#' data after adjusting for the observed covariates. This is useful in 
#' downstream scenarios where hypothesis testing on the original count data 
#' is desired, such as in differential expression.}
#' \item{\code{saver.est}}{Denoised SAVER estimates based on empirical Bayes 
#' method. The SAVER estimate and posterior distribution can be used to 
#' visualize gene expression and perform trajectory analyses. The SAVER 
#' estimate acts as the denoised estimate of the true expression.}
#' \item{\code{saver.se}}{Standard error of SAVER posterior distribution.}
#' \item{\code{saver.samp}}{Sampling from the SAVER posterior distribution. 
#' Sampling from the SAVER posterior allows investigation on gene expression 
#' distribution.}
#' \item{\code{bad.genes}}{Index of genes that could not be predicted accurately 
#' by the autoencoder based on cross validation results. Instead, for these 
#' genes, denoised SAVER estimate are mean gene expression level across cells.}
#' \item{\code{time}}{Time of decoder training.}
#' \item{\code{train2}}{A list that contains training information and model 
#' details.}
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



run_decoder_sc <- function(x, z.mean, z.logvar, B = NULL, sf = NULL,
                           epochs = 50, alpha = 0.01, 
                           return_saver = FALSE,
                           verbose = 2, 
                           batch_size = 64, 
                           save.model = NULL) {
  
  
  clipnorm = 1
  lr = 1e-3
  nfolds = 6
  loss.thres = 0.1
  clear = TRUE
  optimizer = "Adam"
  dec = c(128, 128)
  bias = TRUE
  nonlinear = rep(TRUE, length(dec))
  act = keras::layer_activation_leaky_relu(alpha = 0.01)
  transpose = TRUE
  
  
  if (tensorflow::tf$executing_eagerly())
    tensorflow::tf$compat$v1$disable_eager_execution()
  
  pt <- Sys.time()
  if (transpose) {
    x <- t(x)
  }
  out_dim <- ncol(x)
  n <- nrow(x)
  k <- ncol(z.mean)
  libsize <- rowSums(x)
  if (is.null(sf)) {
    sf <- libsize/mean(libsize)
  }
  sf.mat <- matrix(log(sf), nrow = nrow(x), ncol = 1)
  if (is.null(B)) {
    B <- scale(matrix(log10(libsize), n, 1))
  }
  l <- ncol(B)
  z_input_mean <- layer_input(shape = k, name = "z_input_mean")
  z_input_logvar <- layer_input(shape = k, name = "z_input_logvar")
  z_input <- layer_input(shape = k, name = "z_input")
  b_input <- layer_input(shape = l, name = "b_input")
  sf_input <- layer_input(shape = 1, name = "sf_input")
  
  n_init <- initializer_orthogonal()
  
  act <- create_activations(nonlinear, act)
  
  sampling <- function(z_args) {
    z_mean <- z_args[[1]]
    z_log_var <- z_args[[2]]
    epsilon <- k_random_normal(k_shape(z_mean))
    return(z_mean + k_exp(z_log_var/2)*epsilon)
  }
  
  z_out <- layer_lambda(f = sampling)(c(z_input_mean, z_input_logvar))
  
  
  decoder_layers <- make_decoder(dec, act, bias,
                                 n_init, out_dim, name = "d", include_last = FALSE)
  
  decoder11 <- layer_concatenate(list(z_out, b_input))
  decoder12 <- layer_concatenate(list(z_input, b_input))
  
  for (i in seq_along(decoder_layers)) {
    decoder11 <- decoder11 %>% decoder_layers[[i]]()
    decoder12 <- decoder12 %>% decoder_layers[[i]]()
  }
  
  last_layer_mu <- layer_dense(units = ncol(x), activation = NULL, 
                               kernel_initializer = n_init,
                               use_bias = bias, name = paste0("mu_out"))
  
  last_layer_theta <- layer_dense(units = ncol(x), activation = NULL, 
                                  kernel_initializer = n_init,
                                  use_bias = bias, name = paste0("theta_out"))
  
  mu_hat <- decoder11 %>% last_layer_mu
  mu_hat_mean <- decoder12 %>% last_layer_mu
  
  theta_min <- log(1e-6)
  theta_max <- log(1e2)
  
  theta_hat <- decoder11 %>% last_layer_theta %>% 
    k_clip(min_value = theta_min, max_value = theta_max)
  theta_hat_mean <- decoder12 %>% last_layer_theta %>% 
    k_clip(min_value = theta_min, max_value = theta_max)
  
  mu_hat_exp <- activation_exponential(mu_hat)
  mu_hat_mean_exp <- activation_exponential(mu_hat_mean)
  
  theta_hat_exp <- activation_exponential(theta_hat)
  theta_hat_mean_exp <- activation_exponential(theta_hat_mean)
  
  mu_hat_sf <- tf$math$add(mu_hat, sf_input)
  mu_hat_mean_sf <- tf$math$add(mu_hat_mean, sf_input)
  
  mu_hat_exp_sf <- activation_exponential(mu_hat_sf)
  mu_hat_mean_exp_sf <- activation_exponential(mu_hat_mean_sf)
  
  out_hat <- layer_concatenate(list(mu_hat_sf, theta_hat), name = "out")
  out_hat_mean <- layer_concatenate(list(mu_hat_mean_sf, theta_hat_mean), name = "out_mean")
  
  decoder_mu <- keras_model(inputs = c(z_input, b_input),
                            outputs = mu_hat_mean_exp)
  
  decoder_mean <- keras_model(inputs = c(z_input, b_input), 
                              outputs = list(mu_hat_mean_exp, theta_hat_mean_exp))
  
  decoder_samp <- keras_model(inputs = c(z_input_mean, z_input_logvar,
                                         b_input),
                              outputs = list(mu_hat_exp, theta_hat_exp))
  
  model <- keras_model(inputs = c(z_input, z_input_mean, z_input_logvar, b_input, sf_input),
                       outputs = list(out_hat, out_hat_mean))
  
  
  model_opt <- opt(optimizer, clipnorm, lr)
  
  nb_loss <- function(y_true, y_pred) {
    log_mu <- mu_hat_sf
    log_theta <- theta_hat
    mu <- mu_hat_exp_sf
    theta <- theta_hat_exp
    f0 <- -tf$math$lgamma(y_true+1)
    f1 <- -tf$math$lgamma(theta)
    f2 <- tf$math$lgamma(y_true + theta)
    f3 <- -(y_true + theta) * tf$math$log(theta + mu)
    f4 <- theta * log_theta
    f5 <- y_true * log_mu
    -k_sum(f0+f1+f2+f3+f4+f5, axis = -1L)
  }
  
  nb_loss0 <- function(y_true, y_pred) {
    log_mu <- mu_hat_mean_sf
    log_theta <- theta_hat_mean
    mu <- mu_hat_mean_exp_sf
    theta <- theta_hat_mean_exp
    f0 <- -tf$math$lgamma(y_true+1)
    f1 <- -tf$math$lgamma(theta)
    f2 <- tf$math$lgamma(y_true + theta)
    f3 <- -(y_true + theta) * tf$math$log(theta + mu)
    f4 <- theta * log_theta
    f5 <- y_true * log_mu
    -k_sum(f0+f1+f2+f3+f4+f5, axis = -1L)
  }
  
  metric_recons_loss <- custom_metric("recons_loss", function(y_true, y_pred) {
    nb_loss(y_true, y_pred)
  })
  
  metric_recons0_loss <- custom_metric("recons0_loss", function(y_true, y_pred) {
    nb_loss0(y_true, y_pred)
  })
  
  Z_blank <- matrix(0, n, k)
  B_blank <- matrix(0, n, l)
  total.loss <- rep(0, epochs)
  lr2 <- rep(0, epochs)
  
  set.seed(5)
  rand_ind <- sample(1:n, n)
  ind_chunk <- split(rand_ind, sort(rand_ind %% nfolds))
  
  w0 <- get_weights(model)
  
  loss <- vector("list", 3)
  test.loss <- matrix(0, out_dim, 3)
  const.loss <- matrix(0, out_dim, 3)
  
  gene.sums <- sapply(1:out_dim, function(i) sum(x[, i]/sf))
  
  for (i in 1:3) {
    test_ind <- ind_chunk[[i]]
    train_ind <- setdiff(rand_ind, test_ind)
    set_weights(model, w0)
    
    test.sum <- sapply(1:out_dim, function(i) sum(x[test_ind, i]/sf[test_ind]))
    gene.means <- (gene.sums - test.sum)/length(train_ind)
    gene.means[gene.means == 0] <- 1e-7
    mu.init <- unname(log(gene.means))
    w1 <- get_weights(get_layer(model, "mu_out"))
    w1[[-1]] <- array(mu.init)
    set_weights(get_layer(model, "mu_out"), w1)
    
    model %>% compile(
      optimizer = model_opt,
      loss = c(nb_loss, nb_loss0),
      loss_weights = c(1-alpha, alpha),
      metrics = c(metric_recons_loss, metric_recons0_loss)
    )
    
    loss[[i]] <- model %>% fit(
      x = list(Z_blank[train_ind, ], z.mean[train_ind, ], z.logvar[train_ind, ],
               B[train_ind, ], sf.mat[train_ind, ]),
      y = list(x[train_ind, ], x[train_ind, ]),
      epochs = epochs,
      batch_size = batch_size,
      verbose = verbose,
      view_metrics = FALSE,
      validation_data = list(list(Z_blank[test_ind, ], z.mean[test_ind, ], 
                                  z.logvar[test_ind, ], B[test_ind, ],
                                  sf.mat[test_ind, ]),
                             list(x[test_ind, ], x[test_ind, ])),
      callbacks = list(callback_early_stopping(patience = 8, verbose = 1,
                                               restore_best_weights = TRUE))
    )
    
    cat("Calculating validation loss...\n")
    
    train_batches <- split(train_ind, ceiling(seq_along(train_ind)/1000))
    test_batches <- split(test_ind, ceiling(seq_along(test_ind)/200))
    
    test.loss.sum <- rep(0, out_dim)
    const.loss.sum <- rep(0, out_dim)
    for (batch in test_batches) {
      X.test <- decoder_mean %>% predict(list(z.mean[batch, , drop = FALSE],
                                              B[batch, , drop = FALSE]))
      XZ.test <- decoder_mean %>% predict(list(z.mean[batch, , drop = FALSE],
                                               B_blank[batch, , drop = FALSE]))
      xz.mu <- colMeans(XZ.test[[1]])
      xz.theta <- colMeans(XZ.test[[2]])
      
      X.test[[1]] <- sweep(X.test[[1]], 1, sf[batch], "*")
      
      xz <- generate_newx(x[batch, ], X.test, XZ.test)
      test.loss.sum <- test.loss.sum + 
        sapply(1:out_dim, function(i) 
          sum(calc_nb_loss(xz[, i], XZ.test[[1]][, i], XZ.test[[2]][, i])))
      const.loss.sum <- const.loss.sum + 
        sapply(1:out_dim, function(i)
          sum(calc_nb_loss(xz[, i], xz.mu[i], xz.theta[i])))
    }
    
    test.loss[, i] <- test.loss.sum/length(train_ind)
    const.loss[, i] <- const.loss.sum/length(train_ind)
  }
  
  val.loss <- lapply(loss, function(x) x$metrics$val_loss)
  min.epochs <- min(sapply(val.loss, function(x) length(x)))
  val.loss <- sapply(val.loss, function(x) x[1:min.epochs])
  val.loss <- as.matrix(val.loss)
  val.epochs <- which.min(rowMeans(val.loss))
  
  const.loss.se <- apply(const.loss, 1, sd)/sqrt(3)
  test.loss.se <- apply(test.loss, 1, sd)/sqrt(3)
  pooled.se <- sqrt((const.loss.se^2 + test.loss.se^2)/2)
  loss.diff <- (rowMeans(const.loss) - rowMeans(test.loss))/pooled.se
  bad.genes <- which(loss.diff < loss.thres)
  
  set_weights(model, w0)
  
  gene.means <- gene.sums/n
  gene.means[gene.means == 0] <- 1e-7
  mu.init <- unname(log(gene.means))
  w1 <- get_weights(get_layer(model, "mu_out"))
  w1[[-1]] <- array(mu.init)
  set_weights(get_layer(model, "mu_out"), w1)
  
  model %>% compile(
    optimizer = model_opt,
    loss = c(nb_loss, nb_loss0),
    loss_weights = c(1-alpha, alpha),
    metrics = c(metric_recons_loss, metric_recons0_loss)
  )
  
  all.loss <- model %>% fit(
    x = list(Z_blank, z.mean, z.logvar, B, sf.mat),
    y = list(x, x),
    epochs = val.epochs,
    batch_size = batch_size,
    verbose = verbose,
    view_metrics = FALSE
  )
  
  cat("Calculating validation loss...")
  
  z.samp <- k_get_value(sampling(list(k_constant(z.mean), k_constant(z.logvar))))
  
  saver.out <- NULL
  
  if (return_saver) {
    saver.out <- generate_saver(decoder_mean, x, z.mean, B, bad.genes)
  }
  
  if (!is.null(save.model)) {
    system(paste("mkdir", save.model))
    save_model_tf(decoder_mean, paste0(save.model, "/decoder"), 
                  include_optimizer = FALSE)
  }
  
  if (clear) {
    k_clear_session()
  }
  
  train.params <- list(loss = all.loss$metrics,
                       test.loss = test.loss,
                       const.loss = const.loss,
                       loss.diff = loss.diff,
                       cv.loss = val.loss,
                       loss.thres = loss.thres,
                       k = k, dec = dec, alpha = alpha,
                       bias = bias,
                       nonlinear = nonlinear,
                       batch_size = batch_size,
                       optimizer = optimizer,
                       time = Sys.time() - pt)
  
  out <- list(z = z.mean, z.var = z.logvar, z.samp = z.samp, B = B,
              sf = sf,
              adj.raw = saver.out$adj.raw,
              saver.est = saver.out$estimate,
              saver.se = saver.out$se,
              saver.samp = saver.out$samp,
              bad.genes = bad.genes,
              time = Sys.time() - pt,
              train2 = train.params)
  out
  
}
