test_that("test bold main functions", {

  prim_tax <- get_bold_taxonomy("Fannyella", descend = TRUE, ask = FALSE)

  expect_true(is(prim_tax, "data.frame"))
  expect_equal(dim(prim_tax)[1], 6)
  expect_equal(dim(prim_tax)[2], 4)

  bold_tab <- download_bold(prim_tax, ask = FALSE)

  expect_true(is(bold_tab, "data.frame"))

})
