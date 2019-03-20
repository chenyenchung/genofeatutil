context("test-scoring")

# score_predictors() ------------------------------------------------------
test_that("score_predictors()", {
  # Loading dummydb
  dummypath <- system.file("extdata", "dummy.gtf",
                           package = "genofeatutil")
  testdb <- make_database(species = "test",
                          gtf.path = dummypath)

  # Dummy result
  dummyresult <- data.frame("predictor" = c("tf_1", "tf_2", "tf_3"),
                            "target" = rep("alias1.2", 3),
                            "Raw_%IncMSE" = c(0.3, 0.1, 1e-3),
                            row.names = NULL, stringsAsFactors = FALSE)

  dummyexpMat <- data.frame("sample1" = c(3, 1, 1, 3),
                            "sample2" = c(1, 1, 3, 2),
                            "sample3" = c(1, 3, 1, 1),
                            row.names = c("tf_1", "tf_2", "tf_3", "alias1.2"),
                            stringsAsFactors = FALSE)
  # saveRDS(dummyexpMat, "inst/extdata/dummyexpmat.rds")

  # Test type error
  expect_error(
    score_predictors(x = TRUE, expMat = dummyexpMat, use.row = TRUE,
                     score.path = system.file("extdata",  "scoremat.feather",
                                              package = "genofeatutil"),
                     motif.path = system.file("extdata", "dummy_motif.tbl",
                                              package = "genofeatutil"),
                     db = testdb)
  )

  # Generate dummy result
  score_result <- score_predictors(x = dummyresult,
                                   expMat = dummyexpMat, use.row = TRUE,
                                   score.path = system.file(
                                     "extdata", "scoremat.feather",
                                     package = "genofeatutil"),
                                   motif.path = system.file(
                                     "extdata", "dummy_motif.tbl",
                                     package = "genofeatutil"),
                                   db = testdb)

  # Column numbers and output type
  expect_equal(dim(score_result)[2], 10)
  expect_output(str(score_result), "data.frame")

  # expMat type check
  expect_error(
    score_predictors(x = dummyresult,
                     expMat = 10, use.row = TRUE,
                     score.path = system.file(
                       "extdata", "scoremat.feather",
                       package = "genofeatutil"),
                     motif.path = system.file(
                       "extdata", "dummy_motif.tbl",
                       package = "genofeatutil"),
                     db = testdb)
  )

  # expMat read rds file
  expect_equal(
    score_predictors(x = dummyresult,
                     expMat = system.file("extdata", "dummyexpmat.rds",
                                          package = "genofeatutil"),
                     expMat.type = "rds",
                     use.row = TRUE,
                     score.path = system.file(
                       "extdata", "scoremat.feather",
                       package = "genofeatutil"),
                     motif.path = system.file(
                       "extdata", "dummy_motif.tbl",
                       package = "genofeatutil"),
                     db = testdb),
    score_result
  )

  # Transposed expMat
  expect_equal(
    score_predictors(x = dummyresult,
                     expMat = t(dummyexpMat),
                     use.row = FALSE,
                     score.path = system.file(
                       "extdata", "scoremat.feather",
                       package = "genofeatutil"),
                     motif.path = system.file(
                       "extdata", "dummy_motif.tbl",
                       package = "genofeatutil"),
                     db = testdb),
    score_result
  )

  # Premade motif list
  motif_list <-get_motif_info(
    score_path = system.file("extdata",  "scoremat.feather",
                             package = "genofeatutil"),
    motif_path = system.file("extdata", "dummy_motif.tbl",
                             package = "genofeatutil"),
    db = testdb
  )
  expect_equal(
    score_predictors(x = dummyresult,
                     expMat = dummyexpMat,
                     use.row = TRUE,
                     score.path = system.file(
                       "extdata", "scoremat.feather",
                       package = "genofeatutil"),
                     motif.list = motif_list),
    score_result
  )
})



# integrate_score() -------------------------------------------------------
test_that("integrate_score()", {
  # Generate dummy data
  t1 <- data.frame("predictor" = c("tf_1", "tf_2", "tf_3"),
                   "target" = c("gene_1", "gene_2", "gene_3"),
                   "MSE" = c(1, 1, 0))
  t2 <- data.frame("predictor" = c("tf_1", "tf_2", "tf_3"),
                   "target" = c("gene_1", "gene_2", "gene_3"),
                   "MSE" = c(1, 0, 1))
  t3 <- data.frame("predictor" = c("tf_1", "tf_2", "tf_3", "tf_4"),
                   "target" = c("gene_1", "gene_2", "gene_3", "gene_1"),
                   "MSE" = c(0, 1, 1, 1))

  integrated_table <- integrate_score(t1 = t1, t2 = t2, t3 = t3,
                                      column.name = "MSE")
  # Type check
  expect_error(integrate_score(t1 = as.matrix(t1), t2 = t2, t3 = t3,
                               column.name = "MSE"))
  # Check fake column name
  expect_error(integrate_score(t1 = t1, t2 = t2, t3 = t3,
                               column.name = "fake"))
  # Check dims and class
  expect_equal(class(integrated_table), "data.frame")
  expect_equal(ncol(integrated_table), 5)

  # Check na.processing
  expect_equivalent(
    is.na(integrate_score(t1 = t1, t2 = t2, t3 = t3,
                          column.name = "MSE", na.zero = FALSE)),
    matrix(data = c(rep(FALSE, 17), TRUE, TRUE, FALSE), ncol = 5, byrow = TRUE)
  )
})


# plot_score() ------------------------------------------------------------

test_that("plot_score()", {
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

  # Type check
  expect_error(plot_score(TRUE))

  # Argument sanity check
  expect_error(plot_score(intedf, plot.type = "boxplot"))
  expect_error(plot_score(intedf, facet = "sky"))
  expect_error(plot_score(intedf, exp.order = c(1, 2, 3)))

  # Filter sanity check
  expect_error(plot_score(intedf, predictors.use = "nothere"))
  expect_error(plot_score(intedf, targets.use = "nothere"))

  # General output check
  expect_equal_to_reference(plot_score(x = intedf, plot.type = "line"),
                            system.file("extdata", "lineplot.rds",
                                        package = "genofeatutil"))
  expect_equal_to_reference(plot_score(x = intedf, plot.type = "heatmap"),
                            system.file("extdata", "heatmap.rds",
                                        package = "genofeatutil"))

})
