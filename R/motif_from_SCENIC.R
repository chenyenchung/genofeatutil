motif_to_tf <- function(x, motif_db, flybase_sym, dict10X) {
  # It reads the motif table from SCENIC
  # and translate motif names to corresponding gene names
  motif_tbl <- read.table(motif_db, header = TRUE, sep = "\t",
                          comment.char = "", stringsAsFactors = FALSE)
  translated <- lapply(x, function(y) {
    match_row <- motif_tbl[ , 1] %in% y
    if (length(match_row) > 0) {
      tfnames <- motif_tbl$gene_name[match_row]
      tfnames <- convert_10Xgenename(name = tfnames, flybase_sym = flybase_sym,
                                     dict10X = dict10X)
      return(tfnames)
    }
  })
  translated <- translated[sapply(translated, function(x) length(x) > 0)]
}

get_motif_info <- function(score_path, motif_path, number = NULL, threshold = NULL,
                           flybase_sym, dict10X, genes.use) {
  # It reads the motif score data from SCENIC
  # and return a list of top N TFs or TFs passing assigned threshold
  ### The database is made in Python (type: int64)
  ### so type warning will be triggered when loaded in R
  scoreMat <- suppressWarnings(read_feather(score_path))
  motif_name <- scoreMat$features
  scoreMat <- dplyr::select(scoreMat, -features)
  gene_motif_list <- apply(scoreMat, 2, function(x) {
    if (!is.null(threshold)) {
      motifs <- motif_name[x > threshold]
      x <- x[x > threshold]
    }
    motifs <- motifs[order(x, decreasing = TRUE)]
    if (!is.null(number)) {
      motifs <- head(motifs, n = number)
    }
    return(motifs)
  })
  gene_motif_list <- gene_motif_list[sapply(gene_motif_list, function(x) length(x) > 0)]
  names(gene_motif_list) <- convert_10Xgenename(name = names(gene_motif_list), flybase_sym = flybase_sym,
                                                dict10X = dict10X)
  gene_motif_list <- gene_motif_list[names(gene_motif_list) %in% genes.use]
  gene_motif_list <- motif_to_tf(x = gene_motif_list, motif_db = motif_path,
                                 flybase_sym = flybase_sym, dict10X = dict10X)
  return(gene_motif_list)
}