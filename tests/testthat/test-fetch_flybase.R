context("test-fetch_flybase")

# get_fbase_url() ---------------------------------------------------------
test_that("URL retrieval", {
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

  # Check if wrong version numbers will be captured
  expect_error(get_fbase_url("FB1000_AA", dirstr = urls))
  expect_error(get_fbase_url("FB1800_13", dirstr = urls))
})

# fetch_flybase() ---------------------------------------------------------
test_that("Custom URL and the retrieved table structure", {
  # To prevent multiple attempts to connect FlyBae FTP, the test is fetching
  # local dummy files
  test_table <- fetch_flybase(
    paths = c("extdata/dummy_fbgn.tsv",
              "extdata/dummy_syno.tsv")
  )
  fbid_header <- c("##gene_symbol", "organism_abbreviation",
                   "primary_FBgn#", "secondary_FBgn#(s)",
                   "annotation_ID", "secondary_annotation_ID(s)")
  syno_header <- c("##primary_FBid", "organism_abbreviation",
                   "current_symbol", "current_fullname",
                   "fullname_synonym(s)", "symbol_synonym(s)")

  expect_equal(colnames(test_table[["fbid"]]), fbid_header)
  expect_equal(colnames(test_table[["syno"]]), syno_header)
})
