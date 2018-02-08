context("test-structure.R")

test_that("converting import works", {
  imported_data <- import_seqData(e_data_filepath = "test_OTU.csv", f_data_filepath = "test_Meta.csv",
                                  e_meta_filepath = "test_Taxa.csv")
  structured_data <- as.seqData(e_data = imported_data$e_data,
                               f_data = imported_data$f_data,
                               e_meta = imported_data$e_meta,
                               edata_cname = imported_data$guessed_edata_cname,
                               fdata_cname = imported_data$guessed_fdata_cname,
                               taxa_cname = imported_data$guessed_taxa_cname,
                               data_type = "rRNA")
  expect_output(str(structured_data), "List")
  expect_is(structured_data, "seqData")
})
