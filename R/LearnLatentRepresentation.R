#' Characterize latent low-dimensional representation
#'
#' This function generates a low-dimensional latent representation of the cells
#' while adjusting for observed covariates
#'
#' This function takes a Seurat object. By default, raw counts are
#' storeed in "counts" \code{slot}, "RNA" \code{assay}. This function generates
#' a low-dimensional representation using SAVERCAT method and stores it
#' in \code{SeuratObject@reductions$saverc}.
#' Trained model is stored in the specified folder.
#' Before running this function, \code{Preprocess} should be run.
#'
#' @param SeuratObject A Seurat object.
#'
#' @param assay.name Name of assay where count data is stored. Default is 'RNA'.
#'
#' @param slot.name Name of slot where raw count data is stored. Default is 'counts'.
#'
#' @param save.model Folder specified to save the model.
#'
#' @param cell_disp Whether to use cell-specific dispersion parameter. TRUE is
#' to use cell-specific dispersion parameter. FALSE is to use batch-specific
#' dispersion parameter. Default is FALSE.
#'
#' @param k Dimension of the latent representation space produced by SAVERCAT.
#' Default is 30.
#'
#' @param alpha Weight of reconstruction loss using only observed covariates B.
#' Default is 0.01.
#'
#' @param epochs Number of epochs to train the model.
#' If unspecified (\code{Default=NULL}), epochs will be 500.
#'
#' @param batch_size Number of samples per gradient update. Default is 64.
#'
#' @param verbose Verbosity mode. 0 is silent, 1 is progress bar, 2 is one line
#' per epoch. Default is 2.
#'
#'
#' @return A Seurat object with dimensional reduction object stored in
#' \code{SeuratObject@reductions$saverc} and with training details stored in
#' \code{SeuratObject@tools$train}. A list of elements in
#' \code{SeuratObject@tools$train} are:
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
#' @import Seurat
#'
#' @export



LearnLatentRepresentation = function(SeuratObject, assay.name = "RNA",
                                     slot.name = "counts",
                                     save.model = NULL, cell_disp = FALSE,
                                     k = 30, alpha = 0.01,
                                     epochs = NULL, batch_size = 64, verbose = 2) {

  validation = FALSE
  val.epochs = 50
  return_saver = FALSE
  preprocess = FALSE
  hvg_genes_only = TRUE
  n_hvg = 3000
  perc.exp = 0.01

  get_countmatrix = GetAssayData(object = SeuratObject[[assay.name]], slot = slot.name)
  out = SeuratObject@tools$preprocess
  x = out$x
  x.norm = out$x.norm
  sf = out$sf
  B = out$B
  hvg = out$hvg

  suppressWarnings({
    saver.cvae.out = run_vae_sc_nb(x, B = B, x.norm = x.norm, sf = sf, hvg = hvg,
                                   preprocess = preprocess,
                                   hvg_genes_only = hvg_genes_only, n_hvg = n_hvg,
                                   perc.exp = perc.exp,
                                   cell_disp = cell_disp,
                                   k = k, alpha = alpha,
                                   epochs = epochs, batch_size = batch_size, verbose = verbose,
                                   save.model = save.model,
                                   validation = validation, val.epochs = val.epochs,
                                   return_saver = return_saver)
    SeuratObject@tools$train = saver.cvae.out
    embeddings = saver.cvae.out$z
    rownames(embeddings) = colnames(get_countmatrix)
    colnames(embeddings) = paste0("saverc_", 1:ncol(embeddings))
    stdev = as.numeric(apply(saver.cvae.out$z, 2, stats::sd))
    SeuratObject[["saverc"]] = CreateDimReducObject(embeddings = embeddings,
                                                    key = "saverc_",
                                                    stdev = stdev,
                                                    assay = assay.name)
  })
  return(SeuratObject)

}
