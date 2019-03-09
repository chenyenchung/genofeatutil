context("test-file_process")

test_that("Loading works", {
  # Create a temporary .RData file to load
  dummy_df <- data.frame("a" = c(1:5),
                         "b" = c("A", "B", "C", "D", "E"))
  df_load <- dummy_df
  tempobj <- tempfile()
  save(df_load, file = tempobj)
  rm(df_load)
  load(tempobj)
  expect_equivalent(load_rdata(file = tempobj), df_load)
})

test_that("Stopping at multiple objects works", {
  # Create a temporary .RData file to load
  tempobj <- tempfile()
  dummy_df <- data.frame("a" = c(1:5),
                         "b" = c("A", "B", "C", "D", "E"))
  dummy_vec <- c(1:100)
  save(dummy_df, dummy_vec, file = tempobj)
  expect_error(load_rdata(file = tempobj))
})

test_that("Versaread calls load_rdata()", {
  # Create a temporary .RData file to load
  tempobj <- tempfile()
  dummy_df <- data.frame("a" = c(1:5),
                         "b" = c("A", "B", "C", "D", "E"))
  save(dummy_df, file = tempobj)
  expect_equivalent(versaread(tempobj, type = "rdata"), load_rdata(tempobj))
})

test_that("Versaread calls readRDS()", {
  # Create a temporary .RData file to load
  tempobj <- tempfile()
  dummy_df <- data.frame("a" = c(1:5),
                         "b" = c("A", "B", "C", "D", "E"))
  saveRDS(dummy_df, file = tempobj)
  expect_equivalent(versaread(tempobj, type = "rds"), readRDS(tempobj))
})

test_that("Versaread calls read.csv()", {
  # Create a temporary .RData file to load
  tempobj <- tempfile()
  dummy_df <- data.frame("a" = c(1:5),
                         "b" = c("A", "B", "C", "D", "E"))
  write.table(dummy_df$b, file = tempobj, sep = ",", row.names = FALSE)
  expect_equivalent(versaread(tempobj, type = "csv"),
                    read.csv(tempobj, header = TRUE, stringsAsFactors = FALSE)[ , 1])
})