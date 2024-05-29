test_that("buildSequences quick testing of records on the ncbi with multiple copies", {
  tax <- get_ncbi_taxonomy("Bathyteuthis abyssicola")
  refdb <- download_ncbi(ncbi_tax = tax, ask = FALSE)

  expect_true(is(refdb, "data.frame"))

})
