test_that("connection_handler works and gives errors", {

  t0 <- Sys.time()

  expect_equal(connection_handler(function() {"a" * 1}, attempts = 5),
               "ERROR: non-numeric argument to binary operator")

  diff <- difftime(Sys.time(), t0, units = "secs")

  expect_true(5 < diff && diff < 6)
})

test_that("This test should work in both multisession and sequential future plans", {

  # this lauches the first iterations, allowing to set up the workers and taking
  # the following overhead time required to opne the first background processes
  ncbi_limit_handler(seq(1, 8, 1), api_rate = 10, function(x) {
    Sys.sleep(0.1)
    1
  })

  t0s <- Sys.time()

  ncbi_limit_handler(seq(1, 7, 1), api_rate = 10, function(x) {
    Sys.sleep(1)
    1
  })

  diffs <- round(difftime(Sys.time(), t0s, units = "secs"))

  if (is(future::plan(), "sequential")) {

    expect_true(diffs >= 8)

  } else if (is(future::plan(), "multisession")) {

    expect_true(diffs <= 5)

  }
})

test_that("check functioning of warnings catcher in connection_handler", {

  f <- function() {connection_handler(fun = function() {
    warning("warning")}, attempts = 1)}

  expect_warning(f())

})
