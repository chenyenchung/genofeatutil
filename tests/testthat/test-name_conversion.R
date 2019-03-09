context("test-name_conversion")
library(dplyr)

# Raw data that ought to be used in multiple tests
fbpath <- "data/fb_synonym_fb_2019_01.tsv"
fbgn_raw <- data.table::fread(file = fbpath, header = TRUE, sep = "\t",
                              skip = 4, data.table = F, stringsAsFactors = F, quote = "")
fbgn_raw <- dplyr::filter(fbgn_raw, organism_abbreviation == "Dmel")
fbgn_raw <- fbgn_raw[grepl("^FBgn", fbgn_raw$`##primary_FBid`), ]


test_that("generate_flybase_sym", {
  fbgn <- generate_flybase_sym(fbpath)

  # The generated object should be a list of two items ("symbol" and "alias")
  expect_output(str(fbgn), "List of 2")

  # Sample 5 symbols from the database
  # and query the raw table to see if result is consistent
  symbols <- sample(names(fbgn[["symbol"]]), size = 5)

  # Normalize the gene names from the raw data
  fbgn_raw$current_symbol <- normalize_genename(fbgn_raw$current_symbol)

  # Query the generated object
  qobj_symbol <- unname(sort(fbgn[["symbol"]][symbols]))

  # Query the raw table
  qraw_symbol <- sort(fbgn_raw$`##primary_FBid`[fbgn_raw$current_symbol %in% symbols])

  expect_equal(qobj_symbol, qraw_symbol)

  # Sample 5 fbids from the raw table
  # and query the raw table to seeif result is consistent
  aliases_id <- sample(fbgn_raw$`##primary_FBid`, size = 5)

  # Prepare raw aliases from the raw table
  raw_aliases <- fbgn_raw$`symbol_synonym(s)`[fbgn_raw$`##primary_FBid` %in% aliases_id]
  raw_aliases <- lapply(raw_aliases, function (x) {
    strsplit(x, ",")[[1]]
  })
  raw_aliases <- unlist(raw_aliases, use.names = FALSE)
  qraw_alias <- normalize_genename(raw_aliases)

  # Query the generated object
  qobj_alias <- names(fbgn[["alias"]])[fbgn[["alias"]] %in% aliases_id]
  expect_equal(qraw_alias, qobj_alias)
})

# The object for query
fbgn <- generate_flybase_sym(fbpath)

test_that(desc = "categorize_gene", {
  # Sample symbols and aliases from raw table
  symbol <- sample(fbgn_raw$current_symbol, 5)
  aliases <- sample(fbgn_raw$`symbol_synonym(s)`, 5)
  aliases <- lapply(aliases, function(x) {
    strsplit(x, ",")[[1]]
  })
  aliases <- unlist(aliases, use.names = F)

  # Some genes that are not there that should trigger a warning.
  not_there <- c("nyu")

  test_case <- c(symbol, aliases, not_there)

  # Expecting a warning because "nyu" is not a gene name
  expect_warning(test_query <- categorize_gene(test_case, flybase_sym = fbgn))
  test_query <- categorize_gene(test_case, flybase_sym = fbgn)

  # Expecting genes to be separated into those match current symbols and those matches aliases
  # Note that some aliases are also current gene names, so the length of categorization results might
  # be different than the test input (some aliases, if also a symbol, will be counted as a symbol)
  # The bottom line is that all those categorized as symbol should be found in symbols, and vice versa for aliases
  expect_output(str(test_query), "List of 2")
  expect_true(Reduce(`&`, test_query[["symbol"]] %in% names(fbgn$symbol)))
  expect_true(Reduce(`&`, test_query[["alias"]] %in% names(fbgn$alias)))
})
