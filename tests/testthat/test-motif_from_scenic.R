context("test-motif_from_scenic")

test_that("motif_to_tf", {
  dummypath <- system.file("extdata", "dummy.gtf",
                           package = "genofeatutil")
  testdb <- make_database(species = "test",
                          gtf.path = dummypath)

  candidate <- list("gene_1" = c("motif_1", "motif_2"),
                    "gene_2" = c("motif_2", "motif_3"))
  motif_to_tf(x = candidate, motif_db = system.file("extdata", "dummy_motif"),
              db = testdb)
})
