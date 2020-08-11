# SAVERCAT

SAVERCAT is a method for dimension reduction and denoising of single-cell gene expression data that can flexibly adjust for arbitrary observed sample- and cell-level covariates. With SAVERCAT, you can:

* Obtain a low-dimensional representation of the cells that is corrected for the effects of batch and other confounding covariates.

* Remove the effects of batch and other confounding covariates in the original high-dimensional gene expression matrix.

* Further denoise the data, that is, remove technical variation that is due to the inherent random sampling introduced during the library preparation and sequencing steps of the experiment.

## Installation

You can install SAVERCAT from github.

```R
install.packages("devtools")
devtools::install_github("Janezjz/SAVERCAT")
library(SAVERCAT)
```

SAVERCAT depends on Python versions of Tensorflow 2.0 and Keras for neural network model training. To install these, you can either call functions from R as follows:

```R
tensorflow::install_tensorflow()
keras::install_keras()
```

Or you can install separately via Python.

## Tutorial
Please refer to this [vignette](http://htmlpreview.github.io/https://github.com/Janezjz/SAVERCAT/master/docs/savercat_tutorial.html) for a quick start tutorial.
