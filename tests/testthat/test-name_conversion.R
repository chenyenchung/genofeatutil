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


# prepare_database() ---------------------
test_that("prepare_database()", {

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

  # GTF sanity check
  expect_true(Reduce(`&`, testdb[["gtf"]]$type == "gene"))

  # FBid sanity check
  colnames_fbid <- c("primary_FBgn#", "secondary_FBgn#(s)")
  expect_true(Reduce(`&`, colnames_fbid %in% colnames(testdb[["fbid"]])))

  # Synonym sanity check
  expect_true(Reduce(`&`, testdb[["syno"]]$organism_abbreviation == "Dmel"))
  colnames_syno <- c("##primary_FBid", "current_symbol", "symbol_synonym(s)")
  expect_true(Reduce(`&`, colnames_syno %in% colnames(testdb[["syno"]])))
})


# generate_flybase_sym() ------------------------------------------------
test_that("generate_flybase_sym()", {
  dummypath <- system.file("extdata", "dummy.gtf",
                           package = "genofeatutil")
  testdb <- prepare_database(species = "test",
                             gtf.path = dummypath)

  # Check proper input
  expect_error(generate_flybase_sym(testdb[names(testdb) != "syno"]))

  testsym <- generate_flybase_sym(testdb)

  # Check output structure
  expect_true(Reduce(`&`, c("symbol_dict", "alias_dict") %in% names(testsym)))
  expect_match(testsym[["symbol_dict"]][1], "FBgn")
  expect_match(testsym[["alias_dict"]][1], "FBgn")
})

# generate_fbid_version_table() -----------------------------------------
test_that("generate_fbid_version_table()", {
  dummypath <- system.file("extdata", "dummy.gtf",
                           package = "genofeatutil")
  testdb <- prepare_database(species = "test",
                             gtf.path = dummypath)

  # Check proper input
  expect_error(generate_fbid_version_table(testdb[names(testdb) != "fbid"]))

  fbtable_plus <- generate_fbid_version_table(testdb)
  # Check output structure
  expect_true("id_dict" %in% names(fbtable_plus))
  expect_false("fbid" %in% names(fbtable_plus))

  # Output sanity check
  value_fbgn <- grepl("^FBgn", fbtable_plus[["id_dict"]])
  names_fbgn <- grepl("^FBgn", names(fbtable_plus[["id_dict"]]))
  expect_true(Reduce(`&`, value_fbgn))
  expect_true(Reduce(`&`, names_fbgn))
})


# generate_gene_mapping() -------------------------------------------------

test_that("generate_gene_mapping", {
  # prepare_database()
  dummypath <- system.file("extdata", "dummy.gtf",
                           package = "genofeatutil")
  testdb <- prepare_database(species = "test",
                             gtf.path = dummypath)
  testdb <- generate_fbid_version_table(testdb)

  # Prepare dummy data
  notgtf <- testdb[names(testdb) != "gtf"]
  notiddict <- testdb[names(testdb) != "id_dict"]

  # Unmapped gene warning
  expect_error(generate_gene_mapping(db = notgtf),
               regexp = "The database list seems to be wrong.")
  expect_error(generate_gene_mapping(db = notiddict),
               regexp = "An id conversion table is required.")
  test_mapping <- generate_gene_mapping(db = testdb)
  expect_true(!"gtf" %in% names(test_mapping))
  expect_true("to_name_dict" %in% names(test_mapping))

  # Result sanity check
  genemap_names <- names(test_mapping[["to_name_dict"]])

  ## The names are expected to be FlyBase gene numbers
  expect_true(Reduce(`&`, grepl("^FBgn", genemap_names)))
})

# update_fbgn() -----------------------------------------------------------

test_that("update_fbgn", {
  # prepare_database()
  dummypath <- system.file("extdata", "dummy.gtf",
                           package = "genofeatutil")
  testdb <- make_database(species = "test",
                          gtf.path = dummypath)
  rawdb <- prepare_database(species = "test",
                            gtf.path = dummypath)


  # update_fbgn
  expect_error(update_fbgn("FBgn0000001", db = rawdb),
               regexp = "The database list seems to be wrong.")
  single_convert <- update_fbgn("FBgn0000012", db = testdb)
  expect_equal(single_convert, "FBgn0000001")
  vector_convert <- update_fbgn(
    c("FBgn0000012", "FBgn0000013", "FBgn0000014", "FBgn0000015",
      "FBgn0000016", "FBgn0000017", "FBgn0000018"),
    db = testdb
  )
  expect_equal(vector_convert, paste0("FBgn000000", seq(7)))
})

# convert_gene_to_fbgn() --------------------------------------------------

test_that("convert_gene_to_fbgn", {
  # prepare_database()
  dummypath <- system.file("extdata", "dummy.gtf",
                           package = "genofeatutil")
  testdb <- make_database(species = "test",
                          gtf.path = dummypath)
  rawdb <- prepare_database(species = "test",
                            gtf.path = dummypath)

  # Reject incorrect db
  expect_error(convert_gene_to_fbgn("gene_1", db = rawdb))

  # Dummy query
  query <- c("gene_1", "alias_2.1", "mt:dummy", "alias_3.3", "gene_21")
  exp_id <- c("FBgn0000001", "FBgn0000002", "FBgn0000003")

  expect_equal(convert_gene_to_fbgn(query[1:4], db = testdb),
               exp_id)

  # Test detection of multiple alias mapping
  expect_warning(convert_gene_to_fbgn("alias_confuse", db = testdb))
})

# distinguish_fbgn() ------------------------------------------------------


test_that("distinguish_fbgn", {
  test <- c("FBgn1",
            "FBgn2",
            "Cha")
  expect_equal(distinguish_fbgn(test), c(TRUE, TRUE, FALSE))
  expect_error(distinguish_fbgn(c(1:10)))
})


# convert_to_gene_name() --------------------------------------------------


test_that("convert_to_genename", {
  # prepare_database()
  dummypath <- system.file("extdata", "dummy.gtf",
                           package = "genofeatutil")
  testdb <- make_database(species = "test",
                          gtf.path = dummypath)
  rawdb <- prepare_database(species = "test",
                            gtf.path = dummypath)

  # Check corrupted database
  expect_error(convert_to_genename("gene_1", db = rawdb))

  # Check null input
  expect_error(convert_to_genename(x = NULL , db = testdb))

  # Check gene name translation
  expect_equivalent(convert_to_genename("FBgn0000020", db = testdb),
               normalize_genename("tf_3"))
  expect_equivalent(convert_to_genename("alias_6.3", testdb,
                                        normalize = FALSE),
               "gene_6")

})