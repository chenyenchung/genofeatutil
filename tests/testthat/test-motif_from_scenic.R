context("test-motif_from_scenic")


# motif_to_tf -------------------------------------------------------------

test_that("motif_to_tf", {
  dummypath <- system.file("extdata", "dummy.gtf",
                           package = "genofeatutil")
  testdb <- make_database(species = "test",
                          gtf.path = dummypath)

  candidate <- list("gene_1" = c("motif_1", "motif_2", "motif_18"),
                    "gene_2" = c("motif_2", "motif_3", "motif_6"))
  expect_equal(
    motif_to_tf(x = candidate,
                motif_db = system.file("extdata", "dummy_motif.tbl",
                                       package = "genofeatutil"),
                db = testdb),
    list("gene_1" = c("gn_tf_1", "gn_tf_1", "gn_tf_5"),
         "gene_2" = c("gn_tf_1", "gn_tf_1", "gn_tf_2"))
  )

  expect_equal_to_reference(
    motif_to_tf(x = list("gene1" = paste("motif", c(1:20), sep = "_")),
                motif_db = system.file("extdata", "dummy_motif.tbl",
                                       package = "genofeatutil"),
                db = testdb),
    file = system.file("extdata", "motif_to_tf.rds",
                       package = "genofeatutil")
  )

})

# get_motif_info() --------------------------------------------------------

test_that("get_motif_info()", {
  dummypath <- system.file("extdata", "dummy.gtf",
                           package = "genofeatutil")
  testdb <- make_database(species = "test",
                          gtf.path = dummypath)
  motif_list <-get_motif_info(
    score_path = system.file("extdata",  "scoremat.feather",
                             package = "genofeatutil"),
    motif_path = system.file("extdata", "dummy_motif.tbl",
                             package = "genofeatutil"),
    db = testdb
  )

  # Baseline type check
  expect_true(
    is.list(motif_list)
  )

  # Consistency check without filtering
  expect_equal_to_reference(motif_list,
                            system.file("extdata", "all_motif_info.rds",
                                        package = "genofeatutil"))

  # Check thresholded of 5
  motif5_list <-get_motif_info(
    score_path = system.file("extdata",  "scoremat.feather",
                             package = "genofeatutil"),
    motif_path = system.file("extdata", "dummy_motif.tbl",
                             package = "genofeatutil"),
    db = testdb, threshold = 5
  )

  expect_equal_to_reference(motif5_list ,
                            system.file("extdata", "t5_motif_info.rds",
                                        package = "genofeatutil"))

  # Check top 5
  motift5_list <-get_motif_info(
    score_path = system.file("extdata",  "scoremat.feather",
                             package = "genofeatutil"),
    motif_path = system.file("extdata", "dummy_motif.tbl",
                             package = "genofeatutil"),
    db = testdb, number = 5
  )

  expect_equal_to_reference(motift5_list ,
                            system.file("extdata", "top5_motif_info.rds",
                                        package = "genofeatutil"))
  })
