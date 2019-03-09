mtry_series <- function(x, y) {
  # default is the default value of mtry of randomForest::randomForest()
  default <- if (!is.null(y) && !is.factor(y)) {max(floor(ncol(x)/3), 1)} else {floor(sqrt(ncol(x)))}
  mtry_series <- c(1/2 * default, default, 2 * default, ncol(x))
  if (1/2 * default > 10) {mtry_series <- c(10, mtry_series)}
  return(mtry_series)
}

gene_cor <- function(data, gene.1, gene.2) {
  # Calculate correlation of genes from an expression matrix
  # The genes are selected based on column names
  gene1exp <- as.numeric(data[ , gene.1])
  gene2exp <- as.numeric(data[ , gene.2])
  rcoef <- cor(gene1exp, gene2exp)
  return(rcoef)
}

get_10Xgenenames <- function(object, assay = "RNA") {
  # Extract the names of all detected genes of an assay in a Seuratv3 object
  multiassayobj <- slot(object, name = "assays")
  singleassayobj <- multiassayobj[[assay]]
  assayobj <- GetAssayData(object = singleassayobj, assay.type = assay, slot = "data")
  genenames_10X <- row.names(assayobj)
  return(genenames_10X)
}