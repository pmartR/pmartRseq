context("test-imports.R")

test_that("import file checks work", {
  expect_error(import_seqData(e_data_filepath = "doesn't matter", f_data_filepath = "doesn't matter"),
               "Unsupported biom format. Must end in .biom, .csv, or .txt")
  expect_error(import_seqData(f_data_filepath = "doesn't matter", e_data_filepath = "something.biom"),
               "Unsupported sample information format. Must end in .fasta, .csv, or .txt")
})

test_that("import works", {
  imported_data <- import_seqData(e_data_filepath = "test_OTU.csv", f_data_filepath = "test_Meta.csv",
                                  e_meta_filepath = "test_Taxa.csv")
  expect_equal(imported_data$guessed_edata_cname, "OTU.ID")
  expect_length(imported_data, 6)
  expect_true(all(c("e_data", "f_data", "e_meta") %in% names(imported_data)))
})
