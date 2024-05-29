test_that("download_ncbi outputs a refdb data.frame", {

  tv_tax <- get_ncbi_taxonomy("Thouarella variabilis")

  refdb_ncbi <- download_ncbi(tv_tax, ask = FALSE)

  expect_true(is(refdb_ncbi, "data.frame"))

})
