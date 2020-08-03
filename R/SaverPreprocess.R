#' SAVERCAT setup
#'
#' This is the preprocessing function used before training SAVERCAT model.
#' It generates covariates matrix, library size normalization factor,
#' highly variable genes and other key information required by SAVERCAT.
#'
#' This function takes a Seurat Object as input. By default, raw counts are
#' storeed in "counts" \code{slot}, "RNA" \code{assay}. Covariates to be
#' adjusted are supposed to be a subset of variable in \code{meta.data}.
#' Names of categorical variables to be adjusted by SAVERCAT are specified by
#' \code{categorical.var} as a vector.
#' Names of continuous variable to be adjusted by SAVERCAT are specified by
#' \code{continous.var} as a vector.
#' Then covariates matrix is constructed for these variables together with log
#' library size and sacled.
#' By default, \code{continous.var=NULL}, \code{categorical.var=NULL}, then only
#' scaled log library size will be taken as covariate matrix.
#' This function calculates highly variable genes per batch using
#' \code{Seurat::FindVariableFeatures} function and takes the union.
#' Size factor is library size diveded by mean library size. Output of this
#' function is stored in "tools" \code{slot}.
#'
#'
#' @param SeuratObject A Seurat object.
#'
#' @param assay.name Name of assay where count data is stored. Default is 'RNA'.
#'
#' @param slot.name Name of slot where raw count data is stored. Default is 'counts'.
#'
#' @param categorical.var Names of categorical variables to be adjusted by
#' SAVERCAT. Default is NULL.
#'
#' @param continous.var Names of continuous variables to be adjusted by
#' SAVERCAT. Default is NULL.
#'
#' @param n_hvg Number of highly variable genes to select. Default is 3000.
#'
#' @param hvg_genes_only Whether to learn latent representation using only
#' highly variable genes. TRUE is to use only highly variable genes. FALSE is
#' to use genes with non-zero expression in at least perc.exp*100 % of cells.
#' Default is TRUE.
#'
#' @param perc.exp Good genes are defined as those with a mean expression of
#' more than perc.exp. We select highly variable genes from good genes.
#'
#' @return A Seurat object with information used to train SAVERCAT stored in
#' \code{@tools$preprocess} as a list with the following components:
#' \item{\code{x}}{A subset of expression count matrix. The rows correspond to
#' cells and colomns correspond to highly variable genes or genes with a mean
#' expression of more than perc.exp.}
#' \item{\code{x.norm}}{Log normailized and standardized gene expression matrix
#' (cells by genes).}
#' \item{\code{sf}}{Cell specific size factors.}
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
#' @import Seurat
#'
#' @export


SaverPreprocess = function(SeuratObject, assay.name = "RNA", slot.name = "counts",
                      categorical.var = NULL, continous.var = NULL,
                      n_hvg = 3000, hvg_genes_only = TRUE, perc.exp = 0.01) {

  get_countmatrix = GetAssayData(object = SeuratObject[[assay.name]], slot = slot.name)
  get_metadata =  SeuratObject[[]]
  expr_mat_gc = as.matrix(get_countmatrix)

  libsize = log10(colSums(expr_mat_gc))
  B.raw = matrix(libsize, ncol = 1)
  colnames(B.raw) = "libsize"
  rownames(B.raw) = colnames(expr_mat_gc)

  if(!is.null(categorical.var)) {
    if(!prod(categorical.var %in% colnames(get_metadata))){
      stop("categorical.var must be a subset of variable names in meta.data")
    }
    for(i in 1:length(categorical.var)) {
      var.ind = factor(get_metadata[, categorical.var[i]])
      var.wide = model.matrix(~ var.ind - 1)
      colnames(var.wide) = levels(var.ind)
      B.raw = cbind(B.raw, var.wide)
    }
  }

  if(!is.null(continous.var)) {
    if(!prod(continous.var %in% colnames(get_metadata))){
      stop("continous.var must be a subset of variable names in meta.data")
    }
    ctsvar.matrix = data.matrix(get_metadata[, continous.var])
    colnames(ctsvar.matrix) = continous.var
    B.raw = cbind(B.raw, ctsvar.matrix)
  }

  B = scale(B.raw)
  out = preprocess_data(expr_mat_gc, B = B, n_hvg = n_hvg,
                        hvg_genes_only = hvg_genes_only, perc.exp = perc.exp)
  out$B.raw = B.raw
  SeuratObject@tools$preprocess = out
  return(SeuratObject)

}
