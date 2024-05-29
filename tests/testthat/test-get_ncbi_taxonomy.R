test_that("get_ncbi_taxonomy returns dataframe", {
  Sys.sleep(1/3)
  tax <- get_ncbi_taxonomy("Thouarella variabilis")
  Sys.sleep(1/3)
  tax2 <- get_ncbi_taxonomy("1593383")

  expect_equal(class(tax), "data.frame")
  expect_equal(class(tax2), "data.frame")

})
