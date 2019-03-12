# normalize_path <- function(path) {
#   # Because this script creates a working directory
#   # so relative paths in settings will fail unless prepended "../"
#   if (!grepl("^\\/", path)) {path <- paste0("../", path)}
#   return(path)
# }



#' Loading .RData Files Allowing Assignment To An Variable
#'
#' \code{\link{save}()} allows saving multiple R objects from the environment
#' for later use. One caveat of its default behavior is that when the .RData
#' file is loaded, it retains its variable name, which could collide with
#' something else that is already in the environment. \code{load_rdata()} loads
#' .RData to the function environment, makes sure their is only one object
#' loaded, and return it to allow assignment to another variable name upon
#' loading.
#' @param file a character string specifying the path of the desired Rdata file
#'
#' @return an R object
#' @export
#'
#' @examples
#' temprdata <- tempfile(fileext = ".rData")
#' save(iris, file = temprdata)
#' iris_reload <- load_rdata(file = temprdata)
load_rdata <- function(file) {
  # Allow assignment of .rdata files arbitrarily to a variable
  # in contrast to its original behavior to a fixed name when it was saved
  ## Record the environment before loading
  current_envir <- c(ls(), "current_envir")
  load(file)
  ## Assume the new variable in the environment is the loaded object
  obj_name <- setdiff(ls(), current_envir)
  ## Remind the user if more than one variable is in the file.
  if (length(obj_name) > 1) {
    stop(paste("This .rdata file contains more than 1 object and is likely to",
               "be an environment. Please load it with load()."))
  }
  ## Return the variable to allow assignment
  object <- get(obj_name)
  return(object)
}

#' Load RDS, RData, or CSV Files
#'
#' \code{versaread()} is mainly for loading saved lists of gene names.
#' For .RData files containing one object or an .RDS file, \code{versaread()}
#' loads and returns the object for you to assign. For CSV files, it expects
#' something containing 1 column and with a header, and will ignore everything
#' after the first column.
#' @param type a character string indicating the type of file being loaded
#' (Options: \code{rds}, \code{rdata}, or \code{csv})
#' @param file a character string indicating the path of the desired file
#'
#' @return a character vector or an R object (if it's an RData or RDS file)
#' @export
#'
#' @examples
#' temprds <- tempfile(fileext = ".rds")
#' saveRDS(iris, file = temprds)
#' iris_reload <- versaread(temprds, type = "rds")
#' tempcsv<- tempfile(fileext = ".csv")
#' write.csv(iris$Species, file = tempcsv)
#' iris_species <- versaread(tempcsv, type = "csv")
versaread <- function(file, type) {
  # Read .rds, .rdata, or .csv files
  # It is only expecting a marker list for .csv files though
  ## Check argument type
  if (!type %in% c("rds", "rdata", "csv")) {
    stop("type should be either 'csv', 'rds', or 'rdata'.")
  }
  if (type == "rds") {
    object <- readRDS(file)
  } else if (type == "rdata") {
    object <- load_rdata(file = file)
  } else {
    object <- utils::read.csv(file = file, header = TRUE,
                              stringsAsFactors = FALSE)
    object <- object[ , 1]
  }
  return(object)
}