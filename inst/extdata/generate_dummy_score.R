dummy_score <- data.frame("features" = paste("motif", seq(20), sep = "_"),
                          "gene_1" = c(rep(0, 19), rep(10, 1)),
                          "gene_2" = c(rep(0, 18), rep(10, 2)),
                          "gene_3" = c(rep(0, 17), rep(10, 3)),
                          "gene_4" = c(rep(0, 16), rep(10, 4)),
                          "gene_5" = c(rep(0, 15), rep(10, 5)),
                          "gene_6" = c(rep(0, 14), rep(10, 6)),
                          stringsAsFactors = FALSE)
feather::write_feather(dummy_score, path = "inst/extdata/scoremat.feather")
