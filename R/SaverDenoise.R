#' Calculate SAVER-CVAE corrected counts, SAVER estimates and sampled SAVER
#' estimates.
#'
#' This function is used for calculating SAVERCAT corrected counts,
#' SAVER estimates and sampled SAVER estimates.
#'
#' This function takes a Seurat object and trained decoder model as input, and
#' returns a Seurat object with three new assays. SAVERCAT corrected counts
#' are stored in \code{SeuratObject[["savercovadj"]]}. Denoised SAVER estimates
#' are stored in \code{SeuratObject[["saverdenoised"]]} A sampling from the
#' SAVER posterior distribution is stored in \code{SeuratObject[["saversamp"]]}.
#'
#'
#' @param SeuratObject A Seurat object.
#'
#' @param assay.name Name of assay where count data is stored. Default is 'RNA'.
#'
#' @param slot.name Name of slot where raw count data is stored. Default is 'counts'.
#'
#' @param save.model Folder where decoder model is stored.
#'
#' @param subset.genes Names of a subset of genes to be predicted.
#' Default is NULL, then all genes will be predicted.
#'
#' @return A Seurat object with three new assays:
#' \item{\code{saverdenoised}}{Denoised SAVER estimates
#' based on empirical Bayes
#' method. The SAVER estimate and posterior distribution can be used to
#' visualize gene expression and perform trajectory analyses. The SAVER
#' estimate acts as the denoised estimate of the true expression}
#' \item{\code{saversamp}}{Sampling from the SAVER posterior
#' distribution which allows investigation on gene expression
#' distribution.}
#' \item{\code{savercovadj}}{SAVER-CVAE corrected counts using
#' cdf matching. This object
#' preserves the technical noise characteristics of scRNA-seq expression count
#' data after adjusting for the observed covariates. It is useful in
#' downstream scenarios where hypothesis testing on the original count data
#' is desired, such as in differential expression.}
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


SaverDenoise = function(SeuratObject, assay.name = "RNA", slot.name = "counts",
                        save.model, subset.genes = NULL) {

  decoder = keras::load_model_tf(save.model)

  get_countmatrix = GetAssayData(object = SeuratObject[[assay.name]], slot = slot.name)
  expr_mat_gc = as.matrix(get_countmatrix) # dim(expr_mat_gc) = (num_gene,num_cell)
  all.genes = rownames(x = SeuratObject)
  
  if (is.null(subset.genes)) {
    gene.index = subset.genes
  } else {
    gene.index = match(subset.genes, all.genes)
  }

  vae.out = SeuratObject@tools$train
  B = SeuratObject@tools$preprocess$B
  saver.out = generate_saver(decoder = decoder, x = expr_mat_gc,
                             z = vae.out$z, B = B,
                             bad.genes = vae.out$bad.genes,
                             subset.genes = gene.index)

  SeuratObject@tools$train$denoise.se = saver.out$se
  SeuratObject[["savercovadj"]] = CreateAssayObject(counts = saver.out$adj.count)
  SeuratObject[["saverdenoised"]] = CreateAssayObject(counts = saver.out$estimate)
  SeuratObject[["saversamp"]] = CreateAssayObject(counts = saver.out$samp)

  return(SeuratObject)

}
