#' A preprocessed Peripheral Blood Mononuclear Cells dataset.
#'
#' The dataset consists of three 10x PBMC datasets from three healthy donors
#' sequenced using three different 10x technologies - VDJ, V2, and V3. In each
#' dataset, we identified three main cell types - monocytes, T/NK cells, and
#' B cells. Unidentified cells were removed after filtering out genes with
#' non-zero expression in less than 2 cells. We then sub-sampled 10,000 genes
#' and 2,000 cells.
#'
#' @format A list of two vectors
#' \describe{
#'   \item{counts}{Raw gene count matrix}
#'   \item{metadata}{Cell level matadata, including technology and celltype.}
#' }
#' @source \url{https://support.10xgenomics.com/single-cell-vdj/datasets/3.0.0/vdj_v1_hs_pbmc2_5gex_protein}
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_v3}
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc8k}
"pbmc_prep"
