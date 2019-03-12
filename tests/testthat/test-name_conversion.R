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
  dummypath <- system.file("extdata", "dummy.gtf",
                           package = "genofeatutil")
  # Unsupported species error
  expect_error(prepare_database(species = "cele",
                                gtf.path = dummypath),
               regexp = "prepare_database does not support")
  testdb <- prepare_database(species = "test",
                             gtf.path = dummypath)
  # Output structure check
  expect_equal(sort(names(testdb)), c("fbid", "gtf", "metadata", "syno"))


  # generate_flybase_sym() ------------------------------------------------

  testsym <- generate_flybase_sym(testdb)
  # Check output structure
  expect_output(str(testsym), "List of 5")
  expect_true(Reduce(`&`, c("symbol_dict", "alias_dict") %in% names(testsym)))
  expect_match(testsym[["symbol_dict"]][1], "FBgn")
  expect_match(testsym[["alias_dict"]][1], "FBgn")
})

test_that("generate_gene_mapping", {
  # prepare_database() --------------------------------------------------
  dummypath <- system.file("extdata", "dummy.gtf",
                           package = "genofeatutil")
  testdb <- prepare_database(species = "test",
                             gtf.path = dummypath)
  testdb <- generate_fbid_version_table(testdb)

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

test_that("update_fbgn", {
  # prepare_database() --------------------------------------------------
  dummypath <- system.file("extdata", "dummy.gtf",
                           package = "genofeatutil")
  testdb <- make_database(species = "test",
                          gtf.path = dummypath)
  rawdb <- prepare_database(species = "test",
                            gtf.path = dummypath)


  # update_fbgn -----------------------------------------------------------
  expect_error(update_fbgn("FBgn0032045", db = rawdb),
               regexp = "The database list seems to be wrong.")
  single_convert <- update_fbgn("FBgn0032045", db = testdb)
  expect_equal(single_convert, "FBgn0262029")
  vector_convert <- update_fbgn(
    c("FBgn0032045", "FBgn0086896", "FBgn0000410", "FBgn0025975",
      "FBgn0032046", "FBgn0051610", "FBgn0069196"),
    db = testdb
  )
  expect_equal(vector_convert, rep("FBgn0262029", 7))
})


test_that("convert_gene_to_fbgn", {
  # prepare_database() --------------------------------------------------
  dummypath <- system.file("extdata", "dummy.gtf",
                           package = "genofeatutil")
  testdb <- make_database(species = "test",
                          gtf.path = dummypath)
  rawdb <- prepare_database(species = "test",
                            gtf.path = dummypath)

  # Reject incorrect db
  expect_error(convert_gene_to_fbgn("CG31610", db = rawdb))

  # Dummy query
  query <- c("d", "CG31610", "mt:dummy", "robo")
  exp_id <- c("FBgn0262029", "FBgn0262029")

  expect_equal(convert_gene_to_fbgn(query[1:3], db = testdb),
               exp_id)

  # Test detection of multiple alias mapping
  dummymulti <- c("FBgn0000001", "FBgn0000002")
  names(dummymulti) <- c("a", "a")
  testdb[["alias_dict"]] <- c(testdb[["alias_dict"]], dummymulti)
  expect_warning(convert_gene_to_fbgn("a", db = testdb))
})

