#' Generate a Tuning Series of \code{mtry} for \code{randomForest()}
#'
#' \code{\link{randomForest}()} implements a default value for mtry determined
#' by the number of predictors used. To examine the influence of mtry on
#' prediction accuracy, \code{make_mtry_series()} takes predictor and response
#' matrices and generate a series of mtry for use in a tune grid.
#' Let default value be N, \code{make_mtry_series()} makes a numeric vector
#' consisting 10, N/2, N, 2N, and the number of predictors provided.
#' the number of predictors.
#' @param x a predictor matrix
#' @param y a response matrix
#'
#' @return a numeric factor
#' @export
#'
#' @examples
#' make_mtry_series(x = iris[ , 1:4], y = iris[ , 5])
make_mtry_series <- function(x, y) {
  # default is the default value of mtry of randomForest::randomForest()
  default <- if (!is.null(y) && !is.factor(y)) {
    max(floor(ncol(x)/3), 1)
    } else {floor(sqrt(ncol(x)))}

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

#' Calculate Correlation Coefficient from an Expression Matrix
#'
#' \code{calculate_gene_cor()} takes an expression matrix. Two genes will be
#' selected by column/row names and Pearson's correlation will be calculated
#' by \code{\link{cor}()}.
#'
#' @param data a matrix or data.frame
#' @param gene1,gene2 a character string gene name that's used as a column
#' name in the expression matrix
#' @param use.row by default FALSE, if TRUE, calculate_gene_cor will take
#' row names as gene names; if FALSE, column names will be seens as gene names
#' @param method a character string indicating which correlation coefficient to
#' compute (Options: "\code{pearson}" (default), "\code{kendall}",
#' or "\code{spearman}")
#'
#' @return a number
#' @export
#'
#' @examples
#' genedf <- data.frame("a" = c(1:10),
#'                      "b" = c(6:14))
#' calculate_gene_cor(genedf, "a", "b")
calculate_gene_cor <- function(data, gene1, gene2,
                               use.row = FALSE, method = "pearson") {
  # Calculate correlation of genes from an expression matrix
  # The genes are selected based on column names
  if (!is.matrix(data) & !is.data.frame(data)) {
    stop(strwrap("calculate_gene_core only takes a matrix or
                 a data frame as input data."))
  }
  if (use.row) {
    gene_names <- row.names(data)
  } else {
    gene_names <- colnames(data)
  }

  # Make sure the gene names are in the expression matrices
  if (!(gene1 %in% gene_names)) {
    stop(gene1, " is not found in the expression matrix.")
  }
  if (!(gene2 %in% gene_names)) {
    stop(gene2, " is not found in the expression matrix.")
  }

  if (use.row) {
    gene1exp <- as.numeric(data[gene1, ])
    gene2exp <- as.numeric(data[gene2, ])
  } else {
    gene1exp <- as.numeric(data[ , gene1])
    gene2exp <- as.numeric(data[ , gene2])
  }
  rcoef <- stats::cor(gene1exp, gene2exp, method = method)
  return(rcoef)
}

#' Get Gene Names from an Expression Matrix or a Seurat object
#'
#' \code{get_gene_names()} extracts gene names from a matrix, data frame, or
#' Seurat object. For a martix or data frame, you need to tell it whether the
#' gene names are in row names or column names. If it's a Seurat object,
#' you'll need to specify the assay you want.
#' @param x an object to extract gene names from
#'
#' @return a character vector
#' @export
#'
#' @examples
#' # Create a dummy expression matrix
#' expression_df <- data.frame(
#' "sample1" = c(1:10),
#' "sample2" = c(1:10),
#' "sample3" = c(1:10)
#' )
#'
#' row.names(expression_df) <- c(paste("gene", c(1:10), sep = ""))
#'
#' # Create Seurat object from the dummy data
#' seuratobj <- Seurat::CreateSeuratObject(counts = expression_df,
#'                                         assay = "RNA")
#' # Get gene names from a data frame
#' get_gene_names(expression_df, use.row = TRUE)
#'
#' # Get gene names from a Seuratv3 object
#' get_gene_names(seuratobj, assay = "RNA")
get_gene_names <- function(x, ...) {
  UseMethod("get_gene_names", x)
}

#' @param x an expression matrix
#' @param use.row a logical constant indicating whether the gene names are in
#' row names or not
#'
#' @rdname get_gene_names
get_gene_names.default <- function(x, use.row = FALSE) {
  # Type check
  if (!(is.data.frame(x) | is.matrix(x))) {
    stop("Sorry, the object is not supported by get_gene_names().")
  }
  if (use.row) {
    return(row.names(x))
  }
  return(colnames(x))
}

#' @param object a Seuratv3 object
#' @param assay a character string indicating an assay in the Seurat object
#' @param slot a character string indicating a data slot in the Seurat object
#' (Options: \code{counts}, \code{data}, or \code{scale.data})
#'
#' @rdname get_gene_names
get_gene_names.Seurat <- function(x, assay = "RNA") {
  # Extract the names of all detected genes of an assay in a Seuratv3 object
  multiassayobj <- methods::slot(x, name = "assays")
  singleassayobj <- multiassayobj[[assay]]
  assayobj <- Seurat::GetAssayData(object = singleassayobj,
                           assay.type = assay, slot = "data")
  genenames_10X <- row.names(assayobj)
  return(genenames_10X)
}