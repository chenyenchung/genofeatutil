#' Translate SCENIC Motif Names to Corresponding Transcription Factor Names
#'
#' @param x a list from \code{\link{get_motif_info}}
#' @param motif_db a character string containing the path the motif list from
#' SCENIC / RcisTarget
#' @param db a list from \code{make_database()}
#'
#' @return a list with the motifs names in each item of the list converted to
#' their corresponding transcription factor names
motif_to_tf <- function(x, motif_db, db) {
  # It reads the motif table from SCENIC
  # and translate motif names to corresponding gene names
  motif_tbl <- data.table::fread(file = motif_db, quote = "",
                                 data.table = FALSE, sep = "\t")

  translated <- lapply(x, function(y) {
    match_row <- motif_tbl[ , 1] %in% y
    if (length(match_row) > 0) {
      tfnames <- motif_tbl$gene_name[match_row]
      tfnames <- convert_to_genename(x =tfnames, db = db)
      return(tfnames)
    }
  })
  translated <- translated[sapply(translated, function(x) length(x) > 0)]
  return(translated)
}


#' Generate a Top Motif List Adjacent to Input Genes Based on SCENIC /
#' RcisTarget Database
#'
#' @param score_path a character string of the path to the scoring database
#' from SCENIC / RcisTarget
#' @param motif_path a character string of the path to the motif database
#' from SCENIC / RcisTarget
#' @param number a number of motifs with the highest score to leave from the
#' motif list
#' @param threshold a number of score threshold to filter the motifs
#' @param genes.use a character vector to leave only the genes of interest (
#' Usually the targets) in the output motif list
#' @param db a list from \code{make_database()}
#'
#' @return a list of genes and related transcription factors based on motif
#' info
#' @export
get_motif_info <- function(score_path, motif_path, number = NULL,
                           threshold = NULL, genes.use = NULL, db) {
  # It reads the motif score data from SCENIC
  # and return a list of top N TFs or TFs passing assigned threshold
  ### The database is made in Python (type: int64)
  ### so type warning will be triggered when loaded in R
  scoreMat <- suppressWarnings(feather::read_feather(score_path))
  motif_name <- scoreMat$features
  scoreMat <- scoreMat[ , colnames(scoreMat) != "features"]

  gene_motif_list <- lapply(as.list(scoreMat), function(x) {
    motifs <- motif_name[x > threshold]
    if (!is.null(threshold)) {
      motifs <- motif_name[x > threshold]
      x <- x[x > threshold]
    }

    # Order the motifs by score for selecting top N
    motifs <- motifs[order(x, decreasing = TRUE)]
    if (!is.null(number)) {
      motifs <- utils::head(motifs, n = number)
    }
    return(motifs)
  })

  # Removing empty entries
  gene_motif_list <- gene_motif_list[sapply(gene_motif_list,
                                            function(x) length(x) > 0)]
  names(gene_motif_list) <- convert_to_genename(names(gene_motif_list), db = db)

  # Slice the motif list and leave only the targets
  gene_motif_list <- gene_motif_list[names(gene_motif_list) %in% genes.use]
  gene_motif_list <- motif_to_tf(x = gene_motif_list, motif_db = motif_path,
                                 db = db)
  return(gene_motif_list)
}