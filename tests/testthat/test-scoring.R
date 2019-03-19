context("test-scoring")

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
