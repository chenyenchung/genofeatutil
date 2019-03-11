context("test-name_conversion")

# * Dummy GTF is in extdata/dummy.gtf

test_that("normalize and denormalize gene names", {

  # normalize_genename() --------------------------------------------------
  # denormalize_genename() ------------------------------------------------
  genename <- "128-gene(F)/"
  normalized <- normalize_genename(genename)
  denormalized <- denormalize_genename(normalized)
  expect_equal(denormalized, substr(start = 2,
                                    stop = nchar(make.names(genename)),
                                    make.names(genename)))
  expect_match(normalized, "^gn_")
  expect_match(denormalized, "^128")
})

test_that("Things that require the database prepared", {
  # prepare_database() --------------------------------------------------

  # Unsupported species error
  expect_error(prepare_database(species = "cele",
                                gtf.path = "extdata/dummy.gtf"),
               regexp = "prepare_database does not support")
  testdb <- prepare_database(species = "test",
                             gtf.path = "extdata/dummy.gtf")
  # Output structure check
  expect_equal(sort(names(testdb)), c("fbid", "gtf", "metadata", "syno"))


  # generate_flybase_sym() ------------------------------------------------

  testsym <- generate_flybase_sym(testdb)
  # Check output structure
  expect_output(str(testsym), "List of 5")
  expect_true(Reduce(`&`, c("symbol_dict", "alias_dict") %in% names(testsym)))
  expect_match(testsym[["symbol_dict"]][1], "FBgn")
  expect_match(testsym[["alias_dict"]][1], "FBgn")


  # generate_gene_mapping() ----------------------------------------------

  # Prepare dummy data
  genes <- c("CR41571", "CG45784")
  ids <- c("FBgn0085804", "FBgn0267431")
  notgtf <- testdb[names(testdb) != "gtf"]

  # Unmapped gene warning
  expect_error(generate_gene_mapping(db = notgtf),
                 regexp = "The database list seems to be wrong.")
  test_mapping <- generate_gene_mapping(db = testdb)
  expect_true(!"gtf" %in% names(test_mapping))
  expect_true("to_name_dict" %in% names(test_mapping))
})

