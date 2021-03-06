% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scoring.R
\name{plot_score}
\alias{plot_score}
\title{Plot a Summary Line Plot or Heatmap from an Integrated Data Frame}
\usage{
plot_score(x, facet = NULL, exp.order = NULL, cols.use = NULL,
  predictors.use = NULL, targets.use = NULL, plot.type = "line",
  title = NULL, k.groups = NULL)
}
\arguments{
\item{x}{a data frame, presumably generated by
\code{\link{integrate_score}()}}

\item{facet}{a character string. Options are "predictor" and "target". If
set, the line plot will be separated into multiple panel according to the
value of predictor or target}

\item{exp.order}{a character vector containing the column names provided to
\code{plot_score()}, indicating the order of them}

\item{cols.use}{a character vector containing 3 colors to indicate the color
used in a heatmap to represent low, mid, and high value respectively}

\item{predictors.use}{a character vector indicating the predictors to plot}

\item{targets.use}{a character vector indicating the predictors to plot}

\item{plot.type}{a character string indicating whether a lineplot
(\code{"line"}), a heatmap (\code{"heatmap"}), or a hierarchy tree
(\code{"tree"})should be generated}

\item{title}{a character string indicating the title of the plot}

\item{k.groups}{a numeric value indicating the number of clusters to form
in hierarchical clustering for the predictor-target pairs}
}
\value{
a \code{ggplot2} object
}
\description{
\code{plot_score()} takes a summary data frame (presumably generated by
\code{\link{integrate_score}()}), and plot a line plot or heatmap to
visualize the change of scores in different experiments.
}
\examples{
# Generating a dummy data frame containing predictor, target, and columns of
# scores
intedf <- data.frame(
  "predictor" = c("tf_1", "tf_2", "tf_3", "tf_4"),
  "target" = c("gene_1", "gene_2", "gene_1", "gene_2"),
  "t1" = c(0, 5, 1, 5),
  "t2" = c(0, 3, 0, 1),
  "t3" = c(0, 6, 2, 6),
  "t4" = c(2, 7, 8, 9),
  "t5" = c(9, 2, 2, 7),
  stringsAsFactors = FALSE
)

# Plot a line plot
plot_score(x = intedf, plot.type = "line")

# Plot a heatmap
plot_score(x = intedf, plot.type = "heatmap")
}
