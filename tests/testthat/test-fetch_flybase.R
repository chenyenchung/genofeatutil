context("test-fetch_flybase")

test_that("URL retrieval OK", {
  # To make sure the URLs retrieved are correct
  urls <- get_fbase_url()
  fbidsep <- strsplit(urls[["fbid"]], "\\/")[[1]]
  synosep <- strsplit(urls[["syno"]], "\\/")[[1]]

  # Check if the former part of the URLS are the same
  expect_equivalent(fbidsep[1:6], synosep[1:6])

  # Check if the URLs are for .gz files
  expect_true(grepl("\\.gz$", fbidsep[8]))
  expect_true(grepl("\\.gz$", synosep[8]))

  # Check if the URLs are for FTP
  expect_equivalent(fbidsep[1], "ftp:")
  expect_equivalent(synosep[1], "ftp:")
})

test_that("Successfully captured wrong version number", {
  expect_error(get_fbase_url("FB1000_AA"))
  expect_error(get_fbase_url("FB1800_13"))
})

test_that("The table structure is correct", {
  test_table <- fetch_flybase()
  fbid_header <- c("##gene_symbol", "organism_abbreviation", "primary_FBgn#",
                   "secondary_FBgn#(s)", "annotation_ID", "secondary_annotation_ID(s)")
  syno_header <- c("##primary_FBid", "organism_abbreviation", "current_symbol",
                   "current_fullname", "fullname_synonym(s)", "symbol_synonym(s)")
  expect_equal(colnames(test_table[["fbid"]]), fbid_header)
  expect_equal(colnames(test_table[["syno"]]), syno_header)
})