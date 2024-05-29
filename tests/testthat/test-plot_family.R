test_that("plot functions work", {

  col_ncbi_tax <- get_ncbi_taxonomy("Mycale acerata", ask = FALSE)
  rec_ncbi <- download_ncbi(col_ncbi_tax, ask = FALSE)

  x1 <- plot_length(rec_ncbi)
  x2 <- plot_primers(rec_ncbi)

  expect_true(is(x1, "ggplot"))
  expect_true(is(x2, "ggplot"))

})
