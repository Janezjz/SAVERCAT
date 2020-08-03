#' SAVER-CVAE characterizes latent low-dimensional representation using negative 
#' binomial loss
#'
#' This is the first step training of SAVER-CVAE using negative binomial loss. 
#' This step aims to find a low-dimensional latent representation of the cells 
#' while adjusting for observed covariates
#'
#' The SVAER-CVAE model is trained using the standardized log scale gene 
#' expression level and observed covariates (if NULL, just using scaled 
#' library size) as input to the encoder. The prior distribution of latent 
#' variable is standard normal. The learned multivariate normal mean is a 
#' low-dimensional representation of each cell and can be used in visualization 
#' and clustering analyses.
#' 
#' @param x Expression count matrix. The rows correspond to genes and colomns 
#' correspond to cells. 
#'
#' @param B Observed covariates matrix. The rows correspond to cells. 
#' Default is NULL, then scaled log library size will be taken as B.
#'
#' @param x.norm A log normailized and standardized gene expression matrix 
#' (cells by genes). If \code{preprocess = TRUE}, \code{x.norm} can be 
#' calculated automatically. Default is NULL. 
#'
#' @param sf Vector of cell specific size factors. If \code{preprocess = TRUE}, 
#' \code{sf} can be calculated automatically. Default is NULL. 
#'
#' @param hvg Names of highly variable genes used for training the model. 
#' If \code{preprocess = TRUE}, \code{hvg} can be calculated automatically.
#' Default is NULL.
#'
#' @param preprocess Whether to do preprocess before training the model. 
#' TRUE to apply the default \code{preprocess_data} function to data and then 
#' use the output to train the model. Default is TRUE.
#' 
#' @param hvg_genes_only Whether to train model using only highly variable 
#' genes. TRUE to use only highly variable genes. FALSE to use genes with 
#' non-zero expression in at least perc.exp*100 % of cells. Default is TRUE.
#' 
#' @param n_hvg Number of highly variable genes to select and then to train the 
#' model. Default is 3000.
#' 
#' @param perc.exp Good genes are defined as those with a mean expression of 
#' more than perc.exp. We select highly variable genes from good 
#' genes. If \code{hvg_genes_only = FALSE}, we train model only using good 
#' genes.
#' 
#' @param cell_disp Whether to use cell-specific dispersion parameter. TRUE is 
#' to use cell-specific dispersion parameter. FALSE is to use batch-specific 
#' dispersion parameter. Default is FALSE.
#' 
#' @param k Dimension of the latent representation space produced by SAVER-CVAE. 
#' Default is 30.
#' 
#' @param alpha Weight of reconstruction loss using only observed covariates B.
#' Default is 0.01.
#' 
#' @param epochs Number of epochs to train the model in the first step which 
#' aims to characterize the latent variables.
#' If unspecified, epochs will be 500.
#' 
#' @param batch_size Number of samples per gradient update. Default is 64.
#' 
#' @param verbose Verbosity mode. 0 is silent, 1 is progress bar, 2 is one line
#' per epoch. Default is 2.
#' 
#' @param save.model Folder specified to save the model.
#' 
#' @param validation Whether to run the second step of SAVER-CVAE of 
#' training another denoising decoder. Default is FALSE.
#' 
#' @param val.epochs Number of epochs to train another decoder in the second 
#' step of denoising. Default is 50.
#' 
#' @param return_saver Whether to calculate corrected expression counts and 
#' SAVER estimates. Default is FALSE.
#' 
#' 
#' @return A list with the following components
#' \item{\code{z}}{Mean of posterior normal distribution, which is
#' the low-dimensional representation of the cells.}
#' \item{\code{z.var}}{Log of variance of posterior normal distribution.}
#' \item{\code{z.samp}}{Random generation from posterior normal distribution of 
#' the latent characterization.}
#' \item{\code{B}}{Covariates matrix used to train the model.}
#' \item{\code{sf}}{Cell specific size factors used.}
#' \item{\code{train1}}{A list that contains training information and model 
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


run_vae_sc_nb <- function(x, B = NULL, x.norm = NULL, sf = NULL, hvg = NULL, 
                          preprocess = TRUE, 
                          hvg_genes_only = TRUE, n_hvg = 3000, perc.exp = 0.01,
                          cell_disp = FALSE, 
                          k = 30, alpha = 0.01,
                          epochs = NULL, batch_size = 64, verbose = 2, 
                          save.model = NULL,
                          validation = FALSE, val.epochs = 50,
                          return_saver = FALSE, ...) {
  
  if (tensorflow::tf$executing_eagerly())
    tensorflow::tf$compat$v1$disable_eager_execution()
  
  kl_max = 1
  enc = c(128, 128) 
  dec = c(128, 128)
  bias = TRUE
  B_in = TRUE
  B_disp = TRUE
  B_enc = FALSE
  B_dec = FALSE
  hvg_out = TRUE
  norm_x = FALSE
  nonlinear = c(rep(TRUE, length(enc)), TRUE, rep(TRUE, length(dec)))
  min.lr = 1e-4
  act = keras::layer_activation_leaky_relu(alpha = 0.01)
  bn_center = c(rep(FALSE, length(enc)), FALSE)
  bn_scale = c(rep(TRUE, length(enc)), TRUE)
  anneal_iter = NULL
  optimizer = "Adam"
  clipnorm = 1 
  clipvalue = 5
  lr = 1e-3
  clear = TRUE
  
  
  pt <- Sys.time()
  if (is.null(anneal_iter)) {
    anneal_iter <- min(10000, ceiling(nrow(x)/batch_size)*50)
  }
  if (preprocess) {
    out <- preprocess_data(x, B = B, n_hvg = n_hvg, hvg_genes_only = hvg_genes_only,
                           perc.exp = perc.exp)
    x <- out$x
    x.norm <- out$x.norm
    sf <- out$sf
    B <- out$B
    hvg <- out$hvg
  } else {
    if (is.null(sf)) {
      libsize <- rowSums(x)
      if (norm_x) {
        sf <- libsize/exp(mean(log(libsize)))
      } else {
        sf <- libsize/mean(libsize)
      }
    }
    if (is.null(B)) {
      B <- scale(matrix(log10(libsize)))
    }
  }
  
  input_dim <- ncol(x.norm)
  out_dim <- ncol(x)
  n <- nrow(x)
  
  # begin change
  if (is.null(epochs)) {
    epochs <- max(min(500, ceiling(5e7/n)), 20)
  }
  # end change
  
  x_input <- layer_input(shape = c(input_dim), name = "x_input")
  sf_input <- layer_input(shape = 1, name = "sf_input")
  
  b_input <- layer_input(shape = c(ncol(B)), name = "B")
  
  z_input <- layer_input(shape = c(k), name = "z_input")
  z_input0 <- layer_input(shape = c(k), name = "z_input_samp")
  
  if (B_disp) {
    disp_input <- layer_input(shape = ncol(B), name = "nb_input")
  } else {
    disp_input <- layer_input(shape = 1, name = "nb_input")
  }
  
  n_init <- initializer_orthogonal()
  
  act <- create_activations(nonlinear, act)
  
  if (B_in) {
    z <- layer_concatenate(list(x_input, b_input))
  } else {
    z <- x_input
  }
  
  for (i in seq_along(enc)) {
    z <- z %>% layer_dense(units = enc[i], activation = act[[i]],
                           use_bias = bias, kernel_initializer = n_init,
                           name = paste0("e", i))
    z <- z %>% layer_batch_normalization(center = bn_center[i], 
                                         scale = bn_scale[i], 
                                         name = paste0("be", i))
    
    # begin change
    if(i==1){
      z <- z %>% layer_dropout(rate = 0.1)
    }
    # end change
    
    if (B_enc) {
      z <- layer_concatenate(list(z, b_input))
    }
  }
  
  z_mean <- z %>% layer_dense(units = k, activation = act[[length(enc)+1]],
                              use_bias = bias, kernel_initializer = n_init,
                              name = "z_mean")
  
  z_log_var <- z %>% layer_dense(units = k, activation = NULL,
                                 use_bias = bias, 
                                 kernel_initializer = initializer_orthogonal(gain = 0.01),
                                 name = "z_var")
  
  z_mean <- z_mean %>% layer_batch_normalization(center = bn_center[length(enc) + 1], 
                                                 scale = bn_scale[length(enc) + 1],
                                                 name = "bz")
  
  sampling <- function(z_args) {
    z_mean <- z_args[[1]]
    z_log_var <- z_args[[2]]
    epsilon <- k_random_normal(k_shape(z_mean))
    return(z_mean + k_exp(z_log_var/2)*epsilon)
  }
  
  
  z_out <- layer_lambda(f = sampling)(c(z_mean, z_log_var))
  
  z_out0 <- layer_lambda(f = sampling)(c(z_input0, z_log_var))
  
  if (B_in) {
    mean_out <- keras_model(inputs = list(x_input, b_input), outputs = z_mean)
    var_out <- keras_model(inputs = list(x_input, b_input), outputs = z_log_var)
    samp_out <- keras_model(inputs = list(x_input, b_input), outputs = z_out)
  } else {
    mean_out <- keras_model(inputs = x_input, outputs = z_mean)
    var_out <- keras_model(inputs = x_input, outputs = z_log_var)
    samp_out <- keras_model(inputs = x_input, outputs = z_out)
  }
  
  
  decoder_layers <- make_decoder(dec = dec, act = act[(length(enc) + 2):length(act)], 
                                 bias = bias, n_init = n_init, out_dim = out_dim,
                                 include_last = FALSE)
  
  decoder11 <- layer_concatenate(list(z_out, b_input))
  decoder12 <- layer_concatenate(list(z_out0, b_input))
  decoder12_mean <- layer_concatenate(list(z_input, b_input))
  
  for (i in seq_along(decoder_layers)) {
    decoder11 <- decoder11 %>% decoder_layers[[i]]() 
    decoder12 <- decoder12 %>% decoder_layers[[i]]()
    decoder12_mean <- decoder12_mean %>% decoder_layers[[i]]()
    if (B_dec) {
      decoder11 <- layer_concatenate(list(decoder11, b_input)) 
      decoder12 <- layer_concatenate(list(decoder12, b_input)) 
      decoder12_mean <- layer_concatenate(list(decoder12_mean, b_input)) 
    }
  }
  
  last_layer_mu <- layer_dense(units = ncol(x), activation = NULL, 
                               kernel_initializer = n_init,
                               use_bias = bias, name = paste0("mu_out"))
  
  last_layer_theta <- layer_dense(units = ncol(x), activation = NULL, 
                                  kernel_initializer = n_init,
                                  use_bias = bias, name = paste0("theta_out"))
  
  mu_hat <- decoder11 %>% last_layer_mu
  mu0_hat <- decoder12 %>% last_layer_mu
  mu_hat_mean <- decoder12_mean %>% last_layer_mu
  
  theta_min <- log(1e-6)
  theta_max <- log(max(x)*10)
  
  if (cell_disp) {
    theta_hat <- decoder11 %>% last_layer_theta %>% 
      k_clip(min_value = theta_min, max_value = theta_max)
    theta0_hat <- decoder12 %>% last_layer_theta %>% 
      k_clip(min_value = theta_min, max_value = theta_max)
    theta_hat_mean <- decoder12_mean %>% last_layer_theta %>% 
      k_clip(min_value = theta_min, max_value = theta_max)
  } else {
    theta_hat <- disp_input %>% last_layer_theta %>% 
      k_clip(min_value = theta_min, max_value = theta_max)
    theta0_hat <- disp_input %>% last_layer_theta %>% 
      k_clip(min_value = theta_min, max_value = theta_max)
    theta_hat_mean <- disp_input %>% last_layer_theta %>% 
      k_clip(min_value = theta_min, max_value = theta_max)
  }
  
  if (norm_x) {
    mu_hat_sf <- mu_hat
    mu_hat_mean_sf <- mu_hat_mean
  } else {
    mu_hat_sf <- tf$math$add(mu_hat, sf_input)
    mu_hat_mean_sf <- tf$math$add(mu_hat_mean, sf_input)
  }
  
  mu_hat_exp_sf <- activation_exponential(mu_hat_sf)
  mu_hat_mean_exp_sf <- activation_exponential(mu_hat_mean_sf)
  
  theta_hat_exp <- activation_exponential(theta_hat)
  theta_hat_mean_exp <- activation_exponential(theta_hat_mean)
  
  out_hat <- layer_concatenate(list(mu_hat_sf, theta_hat), name = "out")
  out_hat_mean <- layer_concatenate(list(mu_hat_mean_sf, theta_hat_mean), 
                                    name = "out_mean")
  
  decoder_mean <- keras_model(inputs = c(z_input, b_input, disp_input), 
                              outputs = list(mu_hat_mean, theta_hat_mean))
  
  model <- keras_model(inputs = c(z_input, x_input, b_input, sf_input, disp_input),
                       outputs = c(out_hat, out_hat_mean))
  
  if (norm_x) {
    gene.means <- colMeans(x)
  } else {
    gene.means <- sapply(1:out_dim, function(i) sum(x[, i]/sf))/n
  }
  gene.means[gene.means == 0] <- 1e-7
  mu.init <- unname(log(gene.means))
  
  w1 <- get_weights(get_layer(model, "mu_out"))
  w1[[-1]] <- array(mu.init)
  set_weights(get_layer(model, "mu_out"), w1)
  
  
  model_opt <- opt(optimizer, clipnorm, lr, clipvalue)
  
  kl_weight <- k_variable(0)
  kl_loss <- function() {
    kl_loss = -0.5 * k_sum(1 + z_log_var - k_square(z_mean) - k_exp(z_log_var), 
                           axis = -1L)
    kl_loss
  }
  
  
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
  
  
  vae_loss <- function(y_true, y_pred, kl_weight) {
    recons_loss = nb_loss(y_true, y_pred)
    kl_loss = -0.5 * k_sum(1 + z_log_var - k_square(z_mean) - k_exp(z_log_var), 
                           axis = -1L)
    return(recons_loss + kl_weight*kl_loss)
  }
  
  vae_loss0 <- function(y_true, y_pred, kl_weight) {
    recons_loss = nb_loss0(y_true, y_pred)
    kl_loss = -0.5 * k_sum(1 + z_log_var - k_square(z_mean) - k_exp(z_log_var), 
                           axis = -1L)
    return(recons_loss + kl_weight*kl_loss)
  }
  
  
  metric_kl_loss <- custom_metric("kl_loss", function(y_true, y_pred) {
    kl_loss()
  })
  
  metric_nb_loss <- custom_metric("nb_loss", function(y_true, y_pred) {
    nb_loss(y_true, y_pred)
  })
  
  metric_nb_loss0 <- custom_metric("nb_loss0", function(y_true, y_pred) {
    nb_loss0(y_true, y_pred)
  })
  
  
  metric_kl_weight <- custom_metric("kl_weight", function(y_true, y_pred) {
    kl_weight
  })
  
  model %>% compile(
    optimizer = model_opt,
    loss = c(function(y_true, y_pred) vae_loss(y_true, y_pred, kl_weight),
             function(y_true, y_pred) vae_loss0(y_true, y_pred, kl_weight)),
    loss_weights = c(1-alpha, alpha),
    metrics = c(metric_nb_loss, metric_nb_loss0, metric_kl_loss, metric_kl_weight)
  )
  
  Z_blank <- matrix(0, nrow(x), k)
  B_blank <- matrix(0, nrow(x), ncol(B))
  
  # warm-up KL
  
  anneal_epochs <- ceiling(anneal_iter/(nrow(x.norm)/batch_size))
  if (anneal_epochs > epochs/2) {
    anneal_iter <- floor(anneal_iter*((epochs/2)/anneal_epochs))
    anneal_epochs <- ceiling(epochs/2)
  }
  
  kl_anneal <- KL_Annealing$new(kl_weight, anneal_iter, kl_max)
  lr_earlystop <- LR_EarlyStop$new(min_lr = min.lr, patience = 5)
  reduce_lr <- callback_reduce_lr_on_plateau(monitor = "loss", factor = 0.4,
                                             cooldown = 5,
                                             verbose = verbose, patience = 5,
                                             min_delta = 1e-6)
  callback1 <- list(kl_anneal)
  callback2 <- list(kl_anneal, reduce_lr, lr_earlystop)
  
  sf.mat <- matrix(log(sf))
  if (B_disp) {
    disp.mat <- B
  } else {
    disp.mat <- matrix(1, nrow = n)
  }
  
  cat("KL annealing for", anneal_epochs, "epochs.\n")
  
  loss <- model %>% fit(
    x = list(Z_blank, x.norm, B, sf.mat, disp.mat),
    y = list(x, x),
    epochs = anneal_epochs,
    batch_size = batch_size,
    verbose = verbose,
    view_metrics = FALSE,
    callbacks = callback1
  )
  loss$metrics$lr <- rep(lr, anneal_epochs)
  
  k_set_value(kl_anneal$kl_weight, kl_max)
  
  loss2 <- model %>% fit(
    x = list(Z_blank, x.norm, B, sf.mat, disp.mat),
    y = list(x, x),
    epochs = epochs,
    initial_epoch = anneal_epochs,
    batch_size = batch_size,
    verbose = verbose,
    view_metrics = FALSE,
    callbacks = callback2
  )
  
  if (!is.null(save.model)) {
    system(paste("mkdir", save.model))
    save_model_tf(mean_out, paste0(save.model, "/encoder_mean"), 
                  include_optimizer = FALSE)
    save_model_tf(var_out, paste0(save.model, "/encoder_var"), 
                  include_optimizer = FALSE)
    save_model_tf(decoder_mean, paste0(save.model, "/decoder_train"), 
                  include_optimizer = FALSE)
  }
  
  z.orig <- mean_out %>% predict(list(x.norm, B))
  z.var.orig <- var_out %>% predict(list(x.norm, B))
  z.samp.orig <- samp_out %>% predict(list(x.norm, B))
  
  z.keep <- which(colMeans(z.var.orig) < -0.01)
  z.mean <- z.orig[, z.keep]
  z.var <- z.var.orig[, z.keep]
  z.samp <- z.samp.orig[, z.keep]
  
  train.params <- list(hvg = hvg,
                       hvg_out = hvg_out,
                       norm_x = norm_x,
                       ngenes = out_dim,
                       loss = Map(c, loss$metrics, loss2$metrics),
                       z = z.orig, z.var = z.var.orig,
                       z.samp = z.samp.orig,
                       k = k, enc = enc, dec = dec, alpha = alpha,
                       bias = bias,
                       nonlinear = nonlinear,
                       bn_center = bn_center,
                       bn_scale = bn_scale,
                       batch_size = batch_size,
                       optimizer = optimizer,
                       anneal_iter = anneal_iter,
                       time = Sys.time() - pt)
  
  if (validation) {
    
    out <- run_decoder_sc(x, z.mean, z.var, B, sf = sf, 
                          epochs = val.epochs, alpha = alpha,
                          return_saver = return_saver,
                          verbose = verbose,
                          batch_size = batch_size,
                          save.model = save.model)
    out$train1 <- train.params
  } else {
    out <- list(z = z.mean, z.var = z.var, z.samp = z.samp, B = B,
                sf = sf,
                train1 = train.params)
  }
  
  if (clear) {
    k_clear_session()
  }
  
  out
}

