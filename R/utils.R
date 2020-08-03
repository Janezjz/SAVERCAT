geo_mean <- function(x) {
  exp(mean(log(x)))
}

calc.cor <- function(x1, x2) {
  vec <- rep(0, ncol(x1))
  for (i in 1:ncol(x1)) {
    vec[i] <- cor(x1[, i], x2[, i])
  }
  return(vec)
}

MSE <- function(x1, x2, scale = FALSE) {
  if (scale) {
    mean((scale(x1)-scale(x2))^2)
  } else {
    mean((x1-x2)^2)
  }
}

calc_PCs <- function(x, n = 16) {
  svd <- irlba(x, n)
  return(x%*%svd$v)
}

calc_tsne <- function(x, low_dim = FALSE, n = 16) {
  if (low_dim) {
    pcs <- x
  } else {
    pcs <- calc_PCs(x)
  }
  return(Rtsne(pcs, pca = FALSE, check_duplicates = FALSE)$Y)
}

calc_umap <- function(z) {
  uwot::umap(z, n_neighbors = 30, n_components = 2, metric = "cosine", 
             n_epochs = NULL, learning_rate = 1, min_dist = 0.3,
             spread = 1, set_op_mix_ratio = 1, local_connectivity = 1,
             repulsion_strength = 1, negative_sample_rate = 5,
             a = NULL, b = NULL, fast_sgd = FALSE, verbose = TRUE)
}

seurat_v3_hvg <- function(x, nfeatures = 5000) {
  hvg.df <- Seurat::FindVariableFeatures(x, selection.method = "vst",
                                         verbose = FALSE)
  return(order(hvg.df$vst.variance.standardized, decreasing = TRUE)[1:nfeatures])
} 

sample_z <- function(mu, log_var) {
  sapply(1:ncol(mu), function(i) rnorm(nrow(mu),
                                       mu[, i],
                                       exp(log_var[, i]/2)))
}


opt <- function(optimizer, clipnorm = 1, lr = 0.001, clipvalue = 5) {
  if (optimizer == "Adam" | optimizer == "hybrid2") {
    optimizer_adam(lr = lr, clipnorm = clipnorm, clipvalue = clipvalue)
  } else if (optimizer == "Nadam" | optimizer == "hybrid") {
    optimizer_nadam()
  } else if (optimizer == "SGD_momentum") {
    optimizer_sgd(lr = 0.01, momentum = 0.9, clipnorm = 1)
  } else if (optimizer == "SGD_Nesterov") {
    optimizer_sgd(lr = 0.01, momentum = 0.9, nesterov = TRUE, clipnorm = 1)
  }
}

create_activations <- function(nonlinear, act) {
  if (length(act) == 1) {
    act <- lapply(1:length(nonlinear), function(x) act)
  }
  act.list <- vector("list", length(nonlinear))
  for (i in seq_along(nonlinear)) {
    if (nonlinear[i]) {
      act.list[[i]] <- act[[i]]
    } else {
      act.list[[i]] <- layer_activation(activation = NULL)
    }
  }
  act.list
}

make_decoder <- function(dec, act, 
                         bias, n_init, out_dim, name = "d", last_act = NULL,
                         include_last = TRUE) {
  decoder_layers <- list()
  for (i in seq_along(dec)) {
    decoder_layers[[i]] <- layer_dense(units = dec[i], activation = act[[i]],
                                       use_bias = bias, kernel_initializer = n_init,
                                       name = paste0(name, i))
  }
  if (include_last) {
    decoder_layers[[length(dec)+1]] <- layer_dense(units = out_dim, activation = last_act, 
                                                   kernel_initializer = n_init,
                                                   use_bias = bias, name = paste0(name, "_out"))
  }
  return(decoder_layers)
}

make_bn <- function(dec, batch_norm2 = FALSE) {
  bn_layers <- list()
  for (i in seq_along(dec)) {
    bn_layers[[i]] <- layer_batch_normalization(center = batch_norm2, scale = batch_norm2,
                                                name = paste0("bd", i))
  }
  return(bn_layers)
}

sweep_sparse <- function(x, margin, stats, fun = "*") {
  f <- match.fun(fun)
  if (margin == 1) {
    idx <- x@i + 1
  } else {
    idx <- x@j + 1
  }
  x@x <- f(x@x, stats[idx])
  return(x)
}

MatSD <- function(x, dim = 1, means = NULL, ...) {
  if(dim == 1){
    if (is.null(means)) {
      means <- rowMeans(x, ...)
    }
    sqrt(rowSums((x - means)^2, ...)/(dim(x)[2] - 1))
  } else if (dim == 2) {
    if (is.null(means)) {
      means <- colMeans(x, ...)
    }
    sqrt(rowSums((t(x) - means)^2, ...)/(dim(x)[1] - 1))
  } else stop("Please enter valid dimension")
}

LR_EarlyStop <- R6::R6Class(
  "LR_EarlyStop",
  inherit = KerasCallback, 
  public = list(
    min_lr = NA,
    patience = NA,
    factor = NA,
    
    initialize = function(min_lr = 1e-5, patience = 2) {
      self$min_lr <- min_lr
      self$patience <- patience
    },
    on_epoch_end = function(epoch, logs = list()) {
      current = logs[["lr"]]
      if (current <= self$min_lr) {
        if (self$patience > 0) {
          self$patience <- self$patience - 1
        } else {
          self$model$stop_training <- TRUE
        }
      }
    }
  ))

LR_EarlyStop_val <- R6::R6Class(
  "LR_EarlyStop_val",
  inherit = KerasCallback, 
  public = list(
    min_lr = NA,
    patience = NA,
    factor = NA,
    val_loss = NA,
    
    initialize = function(min_lr = 1e-6, patience = 1) {
      self$min_lr <- min_lr
      self$patience <- patience
      self$val_loss <- c(Inf)
    },
    on_epoch_end = function(epoch, logs = list()) {
      cat(epoch, ":", self$val_loss, "\n")
      self$val_loss <- c(self$val_loss, logs[["val_loss"]])
      cat(self$val_loss, "\n")
      current_loss <- logs[["val_loss"]]
      if (current_loss > self$val_loss[epoch+1]) {
        if (self$patience > 0) {
          self$patience <- self$patience - 1
        } else {
          self$model$stop_training <- TRUE
        }
      }
    }
  ))

Extract_Z <- R6::R6Class(
  "Extract_Z",
  inherit = KerasCallback, 
  public = list(
    z_mean = NA,
    z_logvar = NA,
    mean_out = NA,
    var_out = NA,
    decoder = NA,
    x = NA,
    B = NA,
    
    initialize = function(mean_out = NA, var_out = NA, decoder = NA,
                          x = NA, B = NA) {
      self$z_mean <- NA
      self$z_logvar <- NA
      self$mean_out <- mean_out
      self$var_out <- var_out
      self$decoder <- decoder
      self$x <- x
      self$B <- B
    },
    on_epoch_end = function(epoch, logs = list()) {
      curr_loss = logs[["loss"]]
      if (is.na(curr_loss)) {
        self$model$stop_training <- TRUE
      } else {
        self$z_mean <- self$mean_out %>% predict(list(self$x, self$B))
        self$z_logvar <- self$var_out %>% predict(list(self$x, self$B))
      }
    }
  ))

KL_Annealing <- R6::R6Class(
  "KL_Annealing",
  inherit = KerasCallback,
  public = list(
    kl_weight = 0,
    anneal_iter = 1000,
    kl_max = 1,
    
    initialize = function(kl_weight, anneal_iter, kl_max= 1) {
      self$kl_weight <- kl_weight
      self$anneal_iter <- anneal_iter
      self$kl_max <- kl_max
    },
    on_batch_end = function(batch, logs = list()) {
      step <- 1/self$anneal_iter*self$kl_max
      k_set_value(self$kl_weight, min(k_get_value(self$kl_weight) + step, self$kl_max))
    }
  )
)

kl <- function(z.mean, z.var) {
  (-0.5) * mean(1 + z.var - z.mean^2 - exp(z.var))
}

calc_nbinom_weighted <- Vectorize(function(q0, q1, new_mu, new_theta,
                                           p0, p0_new, p1, p1_new) {
  qvec <- q0:q1
  dvec <- c(0, dnbinom((q0+1):(q1-1), mu = new_mu, size = new_theta), 0)
  dvec[1] <- p0_new - p0
  dvec[length(dvec)] <- p1_new - p1
  dvec <- abs(dvec)
  sum(qvec*dvec)/sum(dvec)
})

get_nbinomx <- function(x, mu, theta, new_mu, new_theta) {
  out <- rep(0, length(x))
  p0 <- pnbinom(x-1, mu = mu, size = theta, log.p = TRUE)
  p1 <- pnbinom(x, mu = mu, size = theta, log.p = TRUE)
  q0 <- qnbinom(p0, mu = new_mu, size = new_theta, log.p = TRUE)
  q1 <- qnbinom(p1, mu = new_mu, size = new_theta, log.p = TRUE)
  
  ind1 <- which(q0 == q1 & q1 != Inf)
  ind2 <- which(q0 != q1 & q1 != Inf)
  ind3 <- which(q1 == Inf)
  
  if (length(ind1) > 0) {
    out[ind1] <- q1[ind1]
  }
  
  if (length(ind2) > 0) {
    p0_new <- pnbinom(q0[ind2], mu = new_mu[ind2], size = new_theta[ind2])
    
    q1away <- abs(q0[ind2] - q1[ind2]) == 1
    
    if (sum(q1away) > 0) {
      d1 <- abs(p0_new[q1away] - exp(p0[ind2][q1away]))
      d2 <- abs(exp(p1[ind2][q1away]) - p0_new[q1away])
      out[ind2[q1away]] <- (q0[ind2[q1away]]*d1 + q1[ind2[q1away]]*d2)/(d1+d2)
    }
    
    if (sum(!q1away) > 0) {
      p1_new <- pnbinom(q1[ind2[!q1away]]-1, mu = new_mu[ind2[!q1away]], 
                        size = new_theta[ind2[!q1away]])
      out[ind2[!q1away]] <- calc_nbinom_weighted(
        q0[ind2[!q1away]], q1[ind2[!q1away]], new_mu[ind2[!q1away]], 
        new_theta[ind2[!q1away]], exp(p0[ind2[!q1away]]), p0_new[!q1away], 
        exp(p1[ind2[!q1away]]), p1_new)
      
    }
  }
  
  if (length(ind3) > 0) {
    x1 <- x[ind3]
    mu1 <- mu[ind3]
    theta1 <- theta[ind3]
    new_mu1 <- new_mu[ind3]
    new_theta1 <- new_theta[ind3]
    sdx <- sqrt(mu1 + mu1^2/theta1)
    sdnewx <- sqrt(new_mu1 + new_mu1^2/new_theta1)
    diff <- (x1 - mu1)/sdx*sdnewx
    lower <- (new_mu1 + diff)*1/2
    upper <- (new_mu1 + diff)*2
    for (i in seq_along(ind3)) {
      k <- ind3[i]
      cand <- floor(seq(lower[i], upper[i], length = 10))
      d1 <- dnbinom(x[k], mu = mu[k], size = theta[k], log = TRUE)
      d2 <- dnbinom(cand, mu = new_mu[k], size = new_theta[k], log = TRUE)
      out[k] <- cand[which.min(abs(d1 - d2))]
    }
  }
  round(out, 5)
}

sample_count <- function(xz) {
  xz.floor <- floor(xz)
  xz.diff <- xz - xz.floor
  out_dim <- nrow(xz)
  xz.samp <- sapply(1:ncol(xz), function(i)
    rbinom(out_dim, size = 1, prob = xz.diff[, i]))
  xz.floor + xz.samp
}


generate_newx <- function(x, X, XZ) {
  newx <- get_nbinomx(c(x), c(X[[1]]), c(X[[2]]),
                      c(XZ[[1]]), c(XZ[[2]]))
  if (is.null(dim(x))) {
    newx
  } else {
    matrix(newx, nrow = nrow(x), ncol = ncol(x))
  }
}


predict_decoder <- function(decoder, z, B, x = NULL, k = NULL, scale = FALSE,
                            clip = FALSE, round = FALSE, digit = 5,
                            transpose = TRUE) {
  out <- decoder %>% predict(list(z, B))
  if (scale | !is.null(k)) {
    if (is.null(k)) {
      k <- rowMeans(x)/colMeans(out[[1]])
    }
    mu <- sweep(out[[1]], 2, k, "*")
  } else {
    mu <- out[[1]]
  }
  theta <- out[[2]]
  if (transpose) {
    mu <- t(mu)
    theta <- t(theta)
  }
  dimnames(mu) <- dimnames(x)
  dimnames(theta) <- dimnames(x)
  if (round) {
    ex <- 10^digit
    mu <- floor(mu*ex)/ex
    theta <- ceiling(theta*ex)/ex
    theta[theta == 0] <- 1/ex
  }
  if (clip) {
    theta[theta > 1e2] <- 1e2
  }
  return(list(mu = mu, theta = theta))
}

calc_disp <- function(mu, y.var) {
  theta <- unname(mu^2/(y.var - mu + 1e-7)) 
  theta[theta <= 0] <- 1e2
  theta[theta < 1e-4] <- 1e-4
  theta[theta > 1e2] <- 1e2
  theta
}

calc_scale <- function(decoder, z, B, x.mean) {
  out <- decoder %>% predict(list(z, B))
  x.mean/colMeans(out[[1]])
}

calc_nb_loss <- function(y, mu, theta) {
  mu[mu == 0] <- 1e-7
  theta[theta == 0] <- 1e-7
  log_mu <- log(mu)
  log_theta <- log(theta)
  f0 <- -lgamma(y+1)
  f1 <- -lgamma(theta)
  f2 <- lgamma(y + theta)
  f3 <- -(y + theta) * log(theta + mu)
  f4 <- theta * log_theta
  f5 <- y * log_mu
  -(f0+f1+f2+f3+f4+f5)
}

calc_pois_loss <- function(y, mu) {
  mu[mu == 0] <- 1e-7
  mu - y * log(mu) + lgamma(y+1)
}



