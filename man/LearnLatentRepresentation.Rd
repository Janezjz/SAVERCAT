% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LearnLatentRepresentation.R
\name{LearnLatentRepresentation}
\alias{LearnLatentRepresentation}
\title{Characterize latent low-dimensional representation}
\usage{
LearnLatentRepresentation(
  SeuratObject,
  assay.name = "RNA",
  slot.name = "counts",
  save.model = NULL,
  cell_disp = FALSE,
  k = 30,
  alpha = 0.01,
  epochs = NULL,
  batch_size = 64,
  verbose = 2
)
}
\arguments{
\item{SeuratObject}{A Seurat object.}

\item{assay.name}{Name of assay where count data is stored. Default is 'RNA'.}

\item{slot.name}{Name of slot where raw count data is stored. Default is 'counts'.}

\item{save.model}{Folder specified to save the model.}

\item{cell_disp}{Whether to use cell-specific dispersion parameter. TRUE is
to use cell-specific dispersion parameter. FALSE is to use batch-specific
dispersion parameter. Default is FALSE.}

\item{k}{Dimension of the latent representation space produced by SAVERCAT.
Default is 30.}

\item{alpha}{Weight of reconstruction loss using only observed covariates B.
Default is 0.01.}

\item{epochs}{Number of epochs to train the model.
If unspecified (\code{Default=NULL}), epochs will be 500.}

\item{batch_size}{Number of samples per gradient update. Default is 64.}

\item{verbose}{Verbosity mode. 0 is silent, 1 is progress bar, 2 is one line
per epoch. Default is 2.}
}
\value{
A Seurat object with dimensional reduction object stored in
\code{SeuratObject@reductions$saverc} and with training details stored in
\code{SeuratObject@tools$train}. A list of elements in
\code{SeuratObject@tools$train} are:
\item{\code{z}}{Mean of posterior normal distribution, which is
the low-dimensional representation of the cells.}
\item{\code{z.var}}{Log of variance of posterior normal distribution.}
\item{\code{z.samp}}{Random generation from posterior normal distribution of
the latent characterization.}
\item{\code{B}}{Covariates matrix used to train the model.}
\item{\code{sf}}{Cell specific size factors used.}
\item{\code{train1}}{A list that contains training information and model
details.}
}
\description{
This function generates a low-dimensional latent representation of the cells
while adjusting for observed covariates
}
\details{
This function takes a Seurat object. By default, raw counts are
storeed in "counts" \code{slot}, "RNA" \code{assay}. This function generates
a low-dimensional representation using SAVERCAT method and stores it
in \code{SeuratObject@reductions$saverc}.
Trained model is stored in the specified folder.
Before running this function, \code{Preprocess} should be run.
}
