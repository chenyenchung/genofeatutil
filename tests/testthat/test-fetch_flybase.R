context("test-fetch_flybase")

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

test_that("Custom URL and the retrieved table structure", {
  # Provide URLs directly to prevent multiple access to FlyBase FTP site
  # This test script will call get_fbase_url() for 4 times if we don't
  # specify URL for fetch_flybase() here.
  # As a result, FlyBase FTP declines the connection and our test will fail.
  fburl <- "ftp://ftp.flybase.net/releases/FB2019_01/precomputed_files/"
  test_table <- fetch_flybase(
    urls = c(paste0(fburl, "genes/fbgn_annotation_ID_fb_2019_01.tsv.gz"),
             paste0(fburl, "synonyms/fb_synonym_fb_2019_01.tsv.gz"))
  )
  fbid_header <- c("##gene_symbol", "organism_abbreviation", "primary_FBgn#",
                   "secondary_FBgn#(s)", "annotation_ID", "secondary_annotation_ID(s)")
  syno_header <- c("##primary_FBid", "organism_abbreviation", "current_symbol",
                   "current_fullname", "fullname_synonym(s)", "symbol_synonym(s)")
  expect_equal(colnames(test_table[["fbid"]]), fbid_header)
  expect_equal(colnames(test_table[["syno"]]), syno_header)
})
