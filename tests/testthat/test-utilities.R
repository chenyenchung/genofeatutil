context("test-utilities")

# Test make_mtry_series() --------------
test_that("mtry series", {
  df <- data.frame("x1" = c(1:10),
                   "x2" = c(1:10),
                   "x3" = c(1:10),
                   "x4" = c(1:10),
                   "x5" = c(1:10),
                   "x6" = c(1:10),
                   "x7" = c(1:10),
                   "x8" = c(1:10),
                   "x9" = c(1:10),
                   "yf" = factor(c(1:10)),
                   "yn" = c(1:10)
                   )
  expect_equal(make_mtry_series(x = df[ , 1:9], y = df[ , 10]), c(1, 3, 6, 9))
  expect_equal(make_mtry_series(x = df[ , 1:4], y = df[ , 11]), c(1, 2, 4))
})

# Test calculate_gene_cor() ------------------
test_that("calculate_gene_cor", {
  genedf <- data.frame("a" = c(1:10),
                       "b" = c(6:15))
  genedf_r <- t(genedf)

  # Check consistency with cor()
  ## Pearson
  expect_equal(calculate_gene_cor(genedf, "a", "b"),
               cor(genedf$a, genedf$b))
  ## Kendall
  expect_equal(calculate_gene_cor(genedf, "a", "b", method = "kendall"),
               cor(genedf$a, genedf$b, method = "kendall"))
  ## Spearman
  expect_equal(calculate_gene_cor(genedf, "a", "b", method = "spearman"),
               cor(genedf$a, genedf$b, method = "spearman"))

  ## Check gene names as row names
  expect_equal(calculate_gene_cor(genedf_r, "a", "b", use.row = TRUE),
               cor(genedf$a, genedf$b))
  ## Kendall
  expect_equal(calculate_gene_cor(genedf_r, "a", "b",
                                  method = "kendall", use.row = TRUE),
               cor(genedf$a, genedf$b, method = "kendall"))
  ## Spearman
  expect_equal(calculate_gene_cor(genedf_r, "a", "b",
                                  method = "spearman", use.row = TRUE),
               cor(genedf$a, genedf$b, method = "spearman"))

  # Check non-existent gene name handling
  expect_error(calculate_gene_cor(genedf, "a", "c"))

  # Check data type detection
  expect_error(calculate_gene_cor(genedf$a, "a", "b"))

})

# Test get_gene_names() (S3 generic) ------------------
test_that("get_gene_names", {
  # Create a dummy expression matrix
  expression_df <- data.frame(
    "sample1" = c(1:10),
    "sample2" = c(1:10),
    "sample3" = c(1:10)
  )
  row.names(expression_df) <- c(paste("gene", c(1:10), sep = ""))

  # Create Seurat object from the dummy data
  seuratobj <- Seurat::CreateSeuratObject(counts = expression_df,
                                          assay = "RNA")
  # Check type
  expect_error(get_gene_names(list()))
  # Check data.frame
  expect_equal(get_gene_names(expression_df, use.row = T),
               c(paste("gene", c(1:10), sep = ""))
               )
  # Check matrix
  expect_equal(get_gene_names(t(expression_df)),
               c(paste("gene", c(1:10), sep = ""))
               )

  # Check Seurat object
  expect_equal(get_gene_names(seuratobj,
                              assay = "RNA"),
               c(paste("gene", c(1:10), sep = ""))
               )
})
