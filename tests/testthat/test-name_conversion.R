context("test-name_conversion")

# * Dummy GTF is in extdata/dummy.gtf


test_that("Things that require the database prepared", {
  # prepare_database() --------------------------------------------------

  # Unsupported species error
  expect_error(prepare_database(species = "cele",
                                gtf.path = "extdata/dummy.gtf"),
               regexp = "prepare_database does not support")
  testdb <- prepare_database(species = "test",
                             gtf.path = "extdata/dummy.gtf")
  # Output structure check
  expect_equal(sort(names(testdb)), c("fbid", "gtf", "species", "syno"))


  # generate_flybase_sym() ------------------------------------------------

  testsym <- generate_flybase_sym(testdb)
  # Check output structure
  expect_output(str(testsym), "List of 2")
  expect_equal(names(testsym), c("symbol", "alias"))
  expect_match(testsym[["symbol"]][1], "FBgn")
  expect_match(testsym[["alias"]][1], "FBgn")


  # generate_gene_mapping() ----------------------------------------------

  # Prepare dummy data
  genes <- c("CR41571", "CG45784")
  ids <- c("FBgn0085804", "FBgn0267431")

  # Unmapped gene warning
  expect_warning(generate_gene_mapping(c(genes, "not_here"), db = testdb),
                 regexp = "Some gene names are not found in the GTF file")
  test_mapping <- generate_gene_mapping(genes, db = testdb)
  expect_equivalent(test_mapping, genes)
  expect_equal(names(test_mapping), ids)

})

