#' Scoring Predictors with a Selected Set of Criteria
#'
#' \code{score_predictors()} scores the predictor-target pair according to
#' the percentage of increment of MSE, the dominance of the predictive power (
#' represented by a gap in MSE among the predictors), correlation of expression
#' levels, and whether position-weighted matrix-inferred motif score of the
#' predictor is high near the target.
#'
#' To score the predictors, this function takes motif database from
#' SCENIC (motif analysis), an expression matrix (correlation analysis),
#' and a gene name conversion database from \code{\link{make_database}()}.
#'
#' @param x a data frame containing the predictors, targets, and
#' percentage increment of MSE (or other measures of predictor power)
#' @param gap.ratio a number indicating the threshold of "having a gap in MSE";
#' any \%IncMSE beyond gap.ratio * the whole range of \%IncMSE will be seen as
#' a gap and get +1 in score
#' @param mse.threshold a number indicating the threshold of individual
#' \%IncMSE; any \%IncMSE beyond this threshold gets +1 in score
#' @param cor.threshold a number indicating the threshold of individual
#' correlation of expression of the predictor-target pair; any correlation
#' beyond this thresold gets +1 in score
#' @param weight a numeric vector indicating the weights to sum up different
#' scores; the length is 7 with a default of c(1, 1, 1, 1, 1, 1, 1, 1). The
#' order of seven scores is "Raw\%IncMSE", "gap", "MSE_threshold",
#' "cor_threshold", "Pearson_cor", "motif_score_sum", "motif_score_presence".
#' @param expMat a data frame containing an expression matrix or a character
#' string containing a path to a csv / rds / rdata file containing an
#' expression matrix
#' @param expMat.type a character string indicating the expMat file type if
#' \code{expMat} is a path
#' @param use.row a logical value indicating whether the gene names are in the
#' row names or column names of the expression matrix loaded
#' @param score.path a character string indicating the path to a SCENIC motif
#' score dataset
#' @param motif.path a character string indicating the path to a SCENIC motif
#' list
#' @param top.motif.number a numeric value indicating the number of top motifs
#' in terms of score to keep; the default is NULL and keeping all
#' @param motif.threshold a numeric value indicting the threshold score for
#' motifs; only motifs with scores higher than this value will be kept; the
#' default value is NULL and keeping all
#' @param db a list generated by \code{\link{make_database}()}
#' @param motif.list a list generated by \code{\link{get_motif_info}()}; if
#' this argument is not NULL, it has higher priority than \code{score.path} and
#' \code{motif.path}
#'
#' @return a data frame containing the scores of each target-predictor pair
#' @export
#'
#' @examples
#' # Loading dummydb
#' dummypath <- system.file("extdata", "dummy.gtf", package = "genofeatutil")
#' testdb <- make_database(species = "test", gtf.path = dummypath)
#'
#' # Generate dummy result for demo
#' dummyresult <- data.frame("predictor" = c("tf_1", "tf_2", "tf_3"),
#'                           "target" = rep("alias1.2", 3),
#'                           "Raw_%IncMSE" = c(0.3, 0.1, 1e-3),
#'                           row.names = NULL, stringsAsFactors = FALSE)
#'
#' dummyexpMat <- data.frame("sample1" = c(3, 1, 1, 3),
#'                           "sample2" = c(1, 1, 3, 2),
#'                           "sample3" = c(1, 3, 1, 1),
#'                           row.names = c("tf_1", "tf_2", "tf_3", "alias1.2"),
#'                           stringsAsFactors = FALSE)
#' # Example score data frame
#' score_result <- score_predictors(x = dummyresult,
#'                                  expMat = dummyexpMat,
#'                                  use.row = TRUE,
#'                                  score.path = system.file(
#'                                    "extdata",
#'                                    "scoremat.feather",
#'                                    package = "genofeatutil"),
#'                                  motif.path = system.file(
#'                                    "extdata",
#'                                    "dummy_motif.tbl",
#'                                    package = "genofeatutil"),
#'                                  db = testdb)
score_predictors <- function(x,
                             gap.ratio = 0.3,
                             mse.threshold = 0.1,
                             cor.threshold = 0.3,
                             weight = c(1, 1, 1, 1, 1, 1, 1),
                             expMat = NULL,
                             expMat.type = "csv",
                             use.row = TRUE,
                             motif.list = NULL,
                             score.path = NULL,
                             motif.path = NULL,
                             top.motif.number = NULL,
                             motif.threshold = NULL,
                             db = NULL) {
  # Type check
  if (class(x) != "data.frame") {
    stop("score_predictors() only takes a data.frame.\n")
  }

  # In case the input is not arranaged
  x <- x[order(x$Raw_.IncMSE, decreasing = TRUE), ]

  # Normalize gene names to be consistent
  x$predictor <- normalize_genename(x$predictor)
  x$target <- normalize_genename(x$target)

  # Gap score calculation
  full_range <- x$Raw_.IncMSE[1] - x$Raw_.IncMSE[nrow(x)]
  diff_mse <- diff(x$Raw_.IncMSE)
  ## Make everything before the gap TRUE
  score_gap <- diff_mse > gap.ratio * full_range
  is_gap <- match(TRUE, score_gap)
  if (!is.na(is_gap)) {
    score_gap[c(1:is_gap)] <- TRUE
  }
  gap_score <- as.numeric(c(score_gap, FALSE))

  # Check threshold of %IncMSE
  mse_score <- as.numeric(x$Raw_.IncMSE > mse.threshold)

  # Check if the predictor co-expresses with it's target
  if (!is.null(expMat)) {
    if (class(expMat) == "character") {
      avgexp <- versaread(expMat, type = expMat.type)
    }
    if (class(expMat) %in% c("data.frame", "matrix")) {
      avgexp <- expMat
    }
    if (!class(avgexp) %in% c("data.frame", "matrix")) {
      stop(paste("Please check if the path or data frame / matrix provided",
                 "is valid.\n"))
    }
    # Normalize gene names for expMat
    if (use.row) {
      row.names(avgexp) <- normalize_genename(row.names(avgexp))
    } else {
      colnames(avgexp) <- normalize_genename(colnames(avgexp))
    }


    cor_abs <- vector(length = nrow(x))
    for (rown in seq(nrow(x))) {
      cor_abs[rown] <- calculate_gene_cor(avgexp, gene1 = x$predictor[rown],
                                          gene2 = x$target[rown],
                                          use.row = use.row)
    }

    cor_score <- as.numeric(abs(cor_abs) > cor.threshold)
  }


  # Check if the target has a predictor's motif around it
  if (!is.null(score.path) & !is.null(motif.path)) {
    if (is.null(db)) {
      stop("Please make a database first with make_database().\n")
    }
    genes.use <- unique(x$target)

    motif_db <- get_motif_info(score_path = score.path,
                               motif_path = motif.path,
                               threshold = motif.threshold,
                               genes.use = genes.use,
                               number = top.motif.number,
                               db = db)

    motif_sum_score <- vector(length = nrow(x))
    motif_presence_score <- vector(length = nrow(x))

    for (rown in seq(nrow(x))) {
      this_motif <- motif_db[[x$target[rown]]]
      motif_presence_score[rown] <- as.numeric(
        x$predictor[rown] %in% this_motif)
      motif_sum_score[rown] <- sum(this_motif == x$predictor[rown])
    }
  }

  if (!is.null(motif.list)) {
    # Type check
    if (!is.list(motif.list)) {
      stop("motif.list should be a list generated by get_motif_info().\n")
    }

    motif_db <- motif.list

    motif_sum_score <- vector(length = nrow(x))
    motif_presence_score <- vector(length = nrow(x))

    for (rown in seq(nrow(x))) {
      this_motif <- motif_db[[x$target[rown]]]
      motif_presence_score[rown] <- as.numeric(
        x$predictor[rown] %in% this_motif)
      motif_sum_score[rown] <- sum(this_motif == x$predictor[rown])
    }
  }

  # Reverse modification for gene names
  x$predictor <- denormalize_genename(x$predictor)
  x$target <- denormalize_genename(x$target)



  # Generating a result table
  result_df <- data.frame("predictor" = x$predictor,
                          "target" = x$target,
                          "Raw_.IncMSE" = x$Raw_.IncMSE,
                          "gap" = gap_score,
                          "MSE_threshold" = mse_score,
                          "cor_threshold" = cor_score,
                          "Pearson_cor" = cor_abs,
                          "motif_score_sum" = motif_sum_score,
                          "motif_score_presence" = motif_presence_score,
                          row.names = NULL, stringsAsFactors = F
  )

  result_df$weighted_sum <- rowSums(sweep(x = result_df[ , 3:9], MARGIN = 2,
                                          FUN = `*`, STATS = weight))
  return(result_df)
}

#' Integrate Score Tables
#'
#' \code{integrate_score()} takes one specified column from multiple data
#' frames  containing the scores for prediction-target pair, and integrated
#' them into a summary table, in which each column represents one table, and
#' each row represents a prediction-target pair.
#'
#' @param ... named arguments of dataframes; the names will be used as column
#' names in the result.
#' @param column.name a character string indicating the column name containing
#' the score to merge
#' @param na.zero a logical value; if TRUE, NAs will be replaced by 0
#'
#' @return a data frame containing integrated scores of interest
#' @export
#'
#' @examples
#' # Generate dummy data
#' t1 <- data.frame("predictor" = c("tf_1", "tf_2", "tf_3"),
#'                  "target" = c("gene_1", "gene_2", "gene_3"),
#'                  "MSE" = c(1, 1, 0))
#' t2 <- data.frame("predictor" = c("tf_1", "tf_2", "tf_3"),
#'                  "target" = c("gene_1", "gene_2", "gene_3"),
#'                  "MSE" = c(1, 0, 1))
#' t3 <- data.frame("predictor" = c("tf_1", "tf_2", "tf_3", "tf_4"),
#'                  "target" = c("gene_1", "gene_2", "gene_3", "gene_1"),
#'                  "MSE" = c(0, 1, 1, 1))
#'
#' integrated_table <- integrate_score(t1 = t1, t2 = t2, t3 = t3,
#'                                     column.name = "MSE")
integrate_score <- function(..., column.name, na.zero = TRUE) {
  df_list <- list(...)
  # Type check
    # Check the existance of column names
  for (df in df_list) {
    if (class(df) != "data.frame") {
      stop("integrate_score() only takes data frames.\n")
    }
    if (!column.name %in% names(df)) {
      stop("The column.name is not presented in every data frame.\n")
    }
  }

  result_list <- lapply(names(df_list), function(name) {
    df <- df_list[[name]]
    result <- data.frame("ident" = paste(df$predictor, df$target, sep = "@"),
                         stringsAsFactors = FALSE)
    result[[name]] <-  df[[column.name]]
    return(result)
  })
  merged <- Reduce(function(x, y) {merge(x, y, by = "ident", all = TRUE)},
                   result_list)
  if (na.zero) {
    merged[is.na(merged)] <- 0
  }

  # Extract predictors and targets and make them two columns again
  predictor <- sapply(strsplit(merged$ident, "@"), function(x) {x[1]})
  target <- sapply(strsplit(merged$ident, "@"), function(x) {x[2]})

  merged <- merged[ , colnames(merged) != "ident"]
  merged <- data.frame(predictor, target, merged, stringsAsFactors = FALSE)

  return(merged)
}

#' Plot a Summary Line Plot or Heatmap from an Integrated Data Frame
#'
#' \code{plot_score()} takes a summary data frame (presumably generated by
#' \code{\link{integrate_score}()}), and plot a line plot or heatmap to
#' visualize the change of scores in different experiments.
#'
#' @param x a data frame, presumably generated by
#' \code{\link{integrate_score}()}
#' @param facet a character string. Options are "predictor" and "target". If
#' set, the line plot will be separated into multiple panel according to the
#' value of predictor or target
#' @param exp.order a character vector containing the column names provided to
#' \code{plot_score()}, indicating the order of them
#' @param cols.use a character vector containing 3 colors to indicate the color
#' used in a heatmap to represent low, mid, and high value respectively
#' @param predictors.use a character vector indicating the predictors to plot
#' @param targets.use a character vector indicating the predictors to plot
#' @param plot.type a character string indicating whether a lineplot
#' (\code{"line"}), a heatmap (\code{"heatmap"}), or a hierarchy tree
#' (\code{"tree"})should be generated
#' @param title a character string indicating the title of the plot
#' @param k.groups a numeric value indicating the number of clusters to form
#' in hierarchical clustering for the predictor-target pairs
#'
#' @return a \code{ggplot2} object
#' @export
#'
#' @examples
#' # Generating a dummy data frame containing predictor, target, and columns of
#' # scores
#' intedf <- data.frame(
#'   "predictor" = c("tf_1", "tf_2", "tf_3", "tf_4"),
#'   "target" = c("gene_1", "gene_2", "gene_1", "gene_2"),
#'   "t1" = c(0, 5, 1, 5),
#'   "t2" = c(0, 3, 0, 1),
#'   "t3" = c(0, 6, 2, 6),
#'   "t4" = c(2, 7, 8, 9),
#'   "t5" = c(9, 2, 2, 7),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Plot a line plot
#' plot_score(x = intedf, plot.type = "line")
#'
#' # Plot a heatmap
#' plot_score(x = intedf, plot.type = "heatmap")
plot_score <- function(x, facet = NULL, exp.order = NULL, cols.use = NULL,
                       predictors.use = NULL, targets.use = NULL,
                       plot.type = "line", title = NULL, k.groups = NULL) {
  # Check
  ## Type check
  if (class(x) != "data.frame") {
    stop("plot_score() only processes a data frame.\n")
  }

  ## plot.type check
  if (!plot.type %in% c("line", "heatmap", "tree")) {
    stop(paste("plot.type of plot_score() only supports 'line',",
               "'heatmap', or 'tree'.\n"))
  }

  ## k.groups check
  if (!is.null(k.groups) && !is.numeric(k.groups)) {
    stop("k.groups should be a number.\n")
  }

  if (is.null(k.groups) && !is.null(facet) && facet == "k") {
    stop("k.groups needs to be set before faceting by k.\n")
  }

  ## facet check
  if (!is.null(facet) && !facet %in% c("predictor", "target", "k")) {
    stop("facet of plot_score() only supports 'target' or 'predictor'\n")
  }

  # Filter data
  if (!is.null(targets.use)) {
    # If there's wrong names, throw error
    target_exist <- targets.use %in% x$target
    if(!Reduce(`&`, target_exist)) {
      stop(paste("The following targets.use are not found:",
                 paste(targets.use[!target_exist], collapse = ", "),
                 ".\n"))
    }
    x <- x[x$target %in% targets.use, ]
  }

  if (!is.null(predictors.use)) {
    # If there's wrong names, throw error
    predictor_exist <- predictors.use %in% x$predictor
    if(!Reduce(`&`, predictor_exist)) {
      stop(paste("The following predictors.use are not found:",
                 paste(predictors.use[!predictor_exist], collapse = ", "),
                 ".\n"))
    }
    x <- x[x$predictor %in% predictors.use, ]
  }

  # Distance matrix for hierarchical ordering and clustering
  row.names(x) <- paste(x$predictor, x$target)
  dist_pairs <- stats::dist(x[ , !colnames(x) %in% c("predictor", "target")],
                     method = "euclidean")
  hc_pairs <- stats::hclust(dist_pairs)
  pair_order <- hc_pairs$labels[hc_pairs$order]

  if (!is.null(k.groups)) {
    k_ident <- stats::cutree(hc_pairs, k = k.groups)
  }

  if (plot.type == "tree") {
    base_plot <- ggdendro::ggdendrogram(hc_pairs, rotate = TRUE)
    print(base_plot)
    return(base_plot)
  }

  # Prepare a data.frame to plot
  ## R CMD CHECK workaround
  predictor <- ggplot2::sym("predictor")
  target <- ggplot2::sym("target")
  plot_df <- tidyr::gather(x, key = "exp", value = "score",
                           !! -predictor, !! -target)

  ## Order the colnames
  if (is.null(exp.order)) {
    plot_df$exp <- factor(plot_df$exp)
  } else {
    ## When exp.order is not NULL...
    # Type check
    if (!is.character(exp.order)) {
      stop("exp.order should be a character vector.")
    }

    # Check existence of levels
    level_exist <- exp.order %in% unique(plot_df$exp)
    if (!Reduce(`&`, level_exist)) {
      stop(paste("The following exp.order is not found:"),
           paste(exp.order[!level_exist], collapse = ", "),
           ".\n")
    }


    # Check coverage of levels
    level_cover <- unique(plot_df$exp) %in% exp.order
    if (!Reduce(`&`, level_cover)) {
      stop(paste("The following exp is missed from the exp.order:"),
           paste( unique(plot_df$exp)[!level_cover], collapse = ", "),
           ".\n")
    }

    plot_df$exp <- factor(plot_df$exp, levels = exp.order)
  }
  plot_df$ident <- paste(plot_df$predictor, plot_df$target)
  plot_df$ident <- factor(plot_df$ident, levels = pair_order)
  if (!is.null(k.groups)) {
    plot_df$k <- k_ident[plot_df$ident]
  }



  # Plotting
  ## Line plot
  if (plot.type == "line") {
    # Base line plot generation
    ### R CMD CHECK workaround for non-standard evaluation
    exp <- ggplot2::sym("exp")
    score <- ggplot2::sym("score")
    ident <- ggplot2::sym("ident")

    base_plot <- ggplot2::ggplot(plot_df,
                                 ggplot2::aes(x = !! exp, y = !! score,
                                              color = !! ident,
                                              group = !! ident)) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::labs(color = "Predictor-target pair",
                    x = "", y = "Prediction Score")

    if (!is.null(facet)) {
      facet.var <- ggplot2::sym(facet)
      base_plot <- base_plot +
        ggplot2::facet_wrap(ggplot2::vars(!!!facet.var))
    }
  }

  ## Heatmap
  if (plot.type == "heatmap") {
    # Base line plot generation
    base_plot <- ggplot2::ggplot(plot_df,
                                 ggplot2::aes(x = exp, y = ident,
                                              fill = score)) +
      ggplot2::geom_tile() +
      ggplot2::labs(fill = "Prediction Score", x = "", y = "")
  }

  if (!is.null(cols.use)) {
    if(length(cols.use) > 3) {
      warning("cols.use for heatmap only takes 3 colors to form a gradient.")
      base_plot <- base_plot + ggplot2::scale_fill_gradientn(
        colors = cols.use
      )
    }
  } else {
    base_plot <- base_plot + ggplot2::scale_fill_gradientn(
      colors = c("blue", "white", "red")
    )
  }

  if (!is.null(title)) {
    base_plot <- base_plot + ggplot2::ggtitle(title)
  }

  base_plot <- base_plot + ggplot2::theme_classic()

  print(base_plot)
  return(base_plot)
}
