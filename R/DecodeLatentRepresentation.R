#' SAVERCAT decodes latent low-dimensional representation
#'
#' This function trains a decoder which takes the low-dimensional
#' representation and produces covariate-adjusted and denoised expression data.
#'
#' This function takes a Seurat object with dimensional reduction object in
#' \code{SeuratObject@reductions$saverc}.
#' This function aims to train a decoder model and store it in the specified
#' folder.
#' Before running this function, \code{LearnLatentRepresentation} should be run.
#'
#' @param SeuratObject A Seurat object.
#'
#' @param assay.name Name of assay where count data is stored. Default is 'RNA'.
#'
#' @param slot.name Name of slot where raw count data is stored. Default is 'counts'.
#'
#' @param save.model Folder specified to save the model.
#'
#' @param alpha Weight of reconstruction loss using only observed covariates B.
#' Default is 0.01.
#'
#' @param epochs Number of epochs to train the model.
#' If unspecified (\code{Default=NULL}), epochs will be 50.
#'
#' @param batch_size Number of samples per gradient update. Default is 64.
#'
#' @param verbose Verbosity mode. 0 is silent, 1 is progress bar, 2 is one line
#' per epoch. Default is 2.
#'
#'
#' @return A Seurat object with training details stored in
#' \code{SeuratObject@tools$train}. Trained model is stored in the specified folder.
#' A list of elements in \code{SeuratObject@tools$train} are:
#' \item{\code{z}}{Mean of posterior normal distribution, which is
#' the low-dimensional representation of the cells.}
#' \item{\code{z.var}}{Log of variance of posterior normal distribution.}
#' \item{\code{z.samp}}{Random sampling from posterior normal distribution of
#' the latent characterization.}
#' \item{\code{B}}{Covariates matrix used to train the model.}
#' \item{\code{sf}}{Cell specific size factors used.}
#' \item{\code{bad.genes}}{Index of genes that could not be predicted accurately
#' by the autoencoder based on cross validation results. Instead, for these
#' genes, denoised SAVER estimate are mean gene expression level across cells.}
#' \item{\code{time}}{Time of decoder training.}
#' \item{\code{train1}}{A list that contains training information and model
#' details when characterizing low-dimensional representation.}
#' \item{\code{train2}}{A list that contains training information and model
#' details when decoding the low-dimensional representation .}
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


DecodeLatentRepresentation = function(SeuratObject, assay.name = "RNA", slot.name = "counts",
                                      save.model = NULL,
                                      alpha = 0.01, epochs = 50, batch_size = 64,
                                      verbose = 2) {


  return_saver = FALSE
  sf = NULL

  get_countmatrix = GetAssayData(object = SeuratObject[[assay.name]], slot = slot.name)
  expr_mat_gc = as.matrix(get_countmatrix) # dim(expr_mat_gc) = (num_gene,num_cell)
  saver.cvae.out = SeuratObject@tools$train
  B = SeuratObject@tools$preprocess$B

  suppressWarnings({
    saver.decoder.out = run_decoder_sc(x = expr_mat_gc, z.mean = saver.cvae.out$z,
                                       z.logvar = saver.cvae.out$z.var,
                                       B = B, sf = sf, epochs = epochs, alpha = alpha,
                                       return_saver = return_saver, verbose = verbose,
                                       batch_size = batch_size, save.model = save.model)

  })

  SeuratObject@tools$train = saver.decoder.out
  SeuratObject@tools$train$train1 = saver.cvae.out$train1
  return(SeuratObject)

}
