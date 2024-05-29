test_that("mergeBarcodeOres works", {

  col_ncbi_tax <- get_ncbi_taxonomy("Colossendeis megalonyx", ask = FALSE)
  col_bold_tax <- get_bold_taxonomy("Colossendeis megalonyx", ask = FALSE)

  rec_ncbi <- download_ncbi(col_ncbi_tax, ask = FALSE)
  rec_bold <- download_bold(col_bold_tax, ask = FALSE)

  all_rec <- mergeBarcodeOres(rec_ncbi, rec_bold, resolve.conflicts = TRUE)

  expect_true(is(all_rec, "data.frame"))

})
