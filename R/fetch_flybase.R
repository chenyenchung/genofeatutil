#' Retrieving Download URLs for FlyBase ID conversion table and Symbol Synonym Table
#'
#' \code{get_fbase_url()} calls RCurl to check FlyBase FTP release repository to
#' check available version. If not otherwise specified, it will generate dowdload
#' URLs for table for FlyBase ID and Synonym.
#'
#' @param version A character specifying the version desired (e.g., "FB2019_01".)
#'
#' @return A named list with 2 items ("fbid" and "syno")
#'
#' @examples
#' # This retrieves the URLs for the up-to-date version
#' get_fbase_url()
#'
#' # This gets the FB2019_01 version
#' get_fbase_url(version = "FB2019_01")
get_fbase_url <- function(version = NULL) {
  # Find current FBid mapping and synonyms from FlyBase release
  filelist <- RCurl::getURL("ftp://ftp.flybase.net/releases/")
  filelist <- data.table::fread(text = filelist, fill = TRUE, data.table = FALSE)

  # Default version is the newest version
  if (is.null(version)) {
    version <- filelist[nrow(filelist), ncol(filelist)]
  }
  message(paste0("Retrieving URL for ", version, "."))

  # Check if the version number is legal
  if (!grepl("^FB[0-9]{4}_[0-9]{2}", version)) {
    stop("The format for version number seems to be wrong. It is FB[year]_[month] (e.g., 'FB2019_01').")
  }

  # Extract available version on the server
  ver_list <- filelist[ , 9]
  ver_list <- ver_list[!ver_list %in% c("README", "current")]

  if (!version %in% ver_list) {
    msg <- paste0("The version is not available on FlyBase. Please try the follows: ",
                  paste(ver_list, collapse = ", "), ".")
    stop("The version is not available on FlyBase.")
  }

  # Make URL for the specified version
  rldate <- gsub("^FB", "", version)
  ftpurl_prefix <- paste0("ftp://ftp.flybase.net/releases/", version, "/precomputed_files/")
  fbid_postfix <- paste0("genes/fbgn_annotation_ID_fb_", rldate,".tsv.gz")
  syno_postfix <- paste0("synonyms/fb_synonym_fb_", rldate, ".tsv.gz")
  fbid_url <- paste0(ftpurl_prefix, fbid_postfix)
  syno_url <- paste0(ftpurl_prefix, syno_postfix)
  urls <- c(fbid_url, syno_url)
  names(urls) <- c("fbid", "syno")
  return(urls)
}


#' Downloading FlyBase ID conversion table and Symbol Synonym Table
#'
#' \code{fetch_flybase()} downloads conversion tables for different version of
#' Flybase IDs and symbol synonyms from Flybase. If not otherwise specified,
#' current version of the tables will be downloaded, and returned as a list with
#' 2 items.
#' @param version A character specifying the version desired (e.g., "FB2019_01".)
#'
#' @return A list of 2 data frames
#' @export
#'
#' @examples
#' # This retrieves the URLs for the up-to-date version
#' fetch_flybase()
#'
#' # This gets the FB2019_01 version
#' fetch_flybase(version = "FB2019_01")
fetch_flybase <- function (version = NULL) {
  # Get urls for FB id and FB synonyms for conversion
  urls <- get_fbase_url()
  fbtable_list <- list()
  for (nameattr in c("fbid", "syno")) {
    temp_dlfile <- tempfile(fileext = ".gz")
    download.file(url = urls[[nameattr]], destfile = temp_dlfile)
    unzippath <- R.utils::gunzip(temp_dlfile)
    fbtable_list[[nameattr]] <- data.table::fread(input = unzippath, quote = "",
                                                  sep = "\t", data.table = FALSE)
  }
  return(fbtable_list)
}