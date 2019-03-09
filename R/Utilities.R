#' Generate a tuning series of \code{mtry} for \code{randomForest()}
#'
#' \code{\link{randomForest}()} implements a default value for mtry determined
#' by the number of predictors used. To examine the influence of mtry on
#' prediction accuracy, \code{make_mtry_series()} takes predictor and response
#' matrices and generate a series of mtry for use in a tune grid.
#' Let default value be N, \code{make_mtry_series()} makes a numeric vector
#' consisting 10, N/2, N, 2N, and the number of predictors provided.
#' the number of predictors.
#' @param x Predictor matrix
#' @param y Response matrix
#'
#' @return A numeric factor
#' @export
#'
#' @examples
#' make_mtry_series(x = iris[ , 1:4], y = iris[ , 5])
make_mtry_series <- function(x, y) {
  # default is the default value of mtry of randomForest::randomForest()
  default <- if (!is.null(y) && !is.factor(y)) {max(floor(ncol(x)/3), 1)} else {floor(sqrt(ncol(x)))}

  # Make sure 2N does not exceed the total number of predictors
  if (ncol(x) > 2 * default) {
    mtry_series <- c(floor(1/2 * default), default, 2 * default, ncol(x))
  } else {mtry_series <- c(floor(1/2 * default), default, ncol(x))}

  # If 1/2 * N is larger than 10, add mtry = 10 as a minimum.
  if (1/2 * default > 10) {mtry_series <- c(10, mtry_series)}

  # If 1/2 * N < 1, there would be a 0 in mtry to be removed.
  mtry_series <- setdiff(mtry_series, 0)

  # If the number of predictors are too small, after rounding some duplicates
  # might appear and must be removed.
  mtry_series <- unique(mtry_series)
  return(mtry_series)
}

gene_cor <- function(data, gene.1, gene.2) {
  # Calculate correlation of genes from an expression matrix
  # The genes are selected based on column names
  gene1exp <- as.numeric(data[ , gene.1])
  gene2exp <- as.numeric(data[ , gene.2])
  rcoef <- stats::cor(gene1exp, gene2exp)
  return(rcoef)
}

get_10Xgenenames <- function(object, assay = "RNA") {
  # Extract the names of all detected genes of an assay in a Seuratv3 object
  multiassayobj <- methods::slot(object, name = "assays")
  singleassayobj <- multiassayobj[[assay]]
  assayobj <- GetAssayData(object = singleassayobj, assay.type = assay, slot = "data")
  genenames_10X <- row.names(assayobj)
  return(genenames_10X)
}