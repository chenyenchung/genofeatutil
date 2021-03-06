% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scoring.R
\name{score_predictors}
\alias{score_predictors}
\title{Scoring Predictors with a Selected Set of Criteria}
\usage{
score_predictors(x, gap.ratio = 0.3, mse.threshold = 0.1,
  cor.threshold = 0.3, weight = c(1, 1, 1, 1, 1, 1, 1),
  expMat = NULL, expMat.type = "csv", use.row = TRUE,
  motif.list = NULL, score.path = NULL, motif.path = NULL,
  top.motif.number = NULL, motif.threshold = NULL, db = NULL)
}
\arguments{
\item{x}{a data frame containing the predictors, targets, and
percentage increment of MSE (or other measures of predictor power)}

\item{gap.ratio}{a number indicating the threshold of "having a gap in MSE";
any \%IncMSE beyond gap.ratio * the whole range of \%IncMSE will be seen as
a gap and get +1 in score}

\item{mse.threshold}{a number indicating the threshold of individual
\%IncMSE; any \%IncMSE beyond this threshold gets +1 in score}

\item{cor.threshold}{a number indicating the threshold of individual
correlation of expression of the predictor-target pair; any correlation
beyond this thresold gets +1 in score}

\item{weight}{a numeric vector indicating the weights to sum up different
scores; the length is 7 with a default of c(1, 1, 1, 1, 1, 1, 1, 1). The
order of seven scores is "Raw\%IncMSE", "gap", "MSE_threshold",
"cor_threshold", "Pearson_cor", "motif_score_sum", "motif_score_presence".}

\item{expMat}{a data frame containing an expression matrix or a character
string containing a path to a csv / rds / rdata file containing an
expression matrix}

\item{expMat.type}{a character string indicating the expMat file type if
\code{expMat} is a path}

\item{use.row}{a logical value indicating whether the gene names are in the
row names or column names of the expression matrix loaded}

\item{motif.list}{a list generated by \code{\link{get_motif_info}()}; if
this argument is not NULL, it has higher priority than \code{score.path} and
\code{motif.path}}

\item{score.path}{a character string indicating the path to a SCENIC motif
score dataset}

\item{motif.path}{a character string indicating the path to a SCENIC motif
list}

\item{top.motif.number}{a numeric value indicating the number of top motifs
in terms of score to keep; the default is NULL and keeping all}

\item{motif.threshold}{a numeric value indicting the threshold score for
motifs; only motifs with scores higher than this value will be kept; the
default value is NULL and keeping all}

\item{db}{a list generated by \code{\link{make_database}()}}
}
\value{
a data frame containing the scores of each target-predictor pair
}
\description{
\code{score_predictors()} scores the predictor-target pair according to
the percentage of increment of MSE, the dominance of the predictive power (
represented by a gap in MSE among the predictors), correlation of expression
levels, and whether position-weighted matrix-inferred motif score of the
predictor is high near the target.
}
\details{
To score the predictors, this function takes motif database from
SCENIC (motif analysis), an expression matrix (correlation analysis),
and a gene name conversion database from \code{\link{make_database}()}.
}
\examples{
# Loading dummydb
dummypath <- system.file("extdata", "dummy.gtf", package = "genofeatutil")
testdb <- make_database(species = "test", gtf.path = dummypath)

# Generate dummy result for demo
dummyresult <- data.frame("predictor" = c("tf_1", "tf_2", "tf_3"),
                          "target" = rep("alias1.2", 3),
                          "Raw_\%IncMSE" = c(0.3, 0.1, 1e-3),
                          row.names = NULL, stringsAsFactors = FALSE)

dummyexpMat <- data.frame("sample1" = c(3, 1, 1, 3),
                          "sample2" = c(1, 1, 3, 2),
                          "sample3" = c(1, 3, 1, 1),
                          row.names = c("tf_1", "tf_2", "tf_3", "alias1.2"),
                          stringsAsFactors = FALSE)
# Example score data frame
score_result <- score_predictors(x = dummyresult,
                                 expMat = dummyexpMat,
                                 use.row = TRUE,
                                 score.path = system.file(
                                   "extdata",
                                   "scoremat.feather",
                                   package = "genofeatutil"),
                                 motif.path = system.file(
                                   "extdata",
                                   "dummy_motif.tbl",
                                   package = "genofeatutil"),
                                 db = testdb)
}
