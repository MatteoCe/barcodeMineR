test_that("web_history_parameter export 2 elements list and gives correct values", {
  l <- web_history_parameter(0, 20, 20)

  expect_true(is.list(l))

  l1 <- web_history_parameter(0, 1, 20)

  expect_true(is.null(l1[[1]]))

  l2 <- web_history_parameter(80, 86, 20)

  expect_true(l2[[2]] == 6)

})

test_that("web_history_splitter returns vector with correct values and corrects special cases", {

  splits <- web_history_splitter(110, 20)

  expect_true(all((length(splits) == 6) & (c(0, 20, 40, 60, 80, 100) %in% splits)))

  splits <- web_history_splitter(100, 20)

  expect_true(all((length(splits) == 5) & (c(0, 20, 40, 60, 80) %in% splits)))

})

test_that("ask_user returns all values if parameter ask is disabled", {

  df <- data.frame(num = c(1,2,3,4,5,6),
                   value = c("Carota", "Carota", "Pomodoro", "Patata", "Patata", "Patata"),
                   rank = c("Radice", "Radice", "Frutto", "Tubero", "Tubero", "Tubero"))

  df_res <- ask_user(df, "value", ask = FALSE)

  expect_equal(df, df_res)

})

test_that("select_accession returns all values if parameter ask is disabled", {

  df <- data.frame(GBInterval_accession = c("AAA001.1", "AAA002.1", "AAA003.1", "AAA004.1", "AAA005.1", "AAA006.1"),
                   Gene_value = c("COI", "CO1", NA, NA, "18s", NA),
                   Product_value = c(NA, NA, "Cytochrome oxydase c 1", "16s", NA, "18s product"))

  df_res <- select_accessions(df, ask = FALSE)

  expect_equal(dim(df)[1], length(df_res))

})

test_that("get_lower_tax_rank gives a character vector and return NULL with 'species'", {

  expect_null(get_lower_tax_rank("subspecies"))

  expect_true(is(get_lower_tax_rank("genus"), "character"))

  expect_equal(get_lower_tax_rank("species"), "subspecies")
  expect_equal(get_lower_tax_rank("genus"), "species")
  expect_equal(get_lower_tax_rank("subfamily"), "genus")
  expect_equal(get_lower_tax_rank("family"), "subfamily")
  expect_equal(get_lower_tax_rank("order"), "family")
  expect_equal(get_lower_tax_rank("class"), "order")

  expect_equal(get_lower_tax_rank("subspecies", upper = TRUE), "species")
  expect_equal(get_lower_tax_rank("species", upper = TRUE), "genus")
  expect_equal(get_lower_tax_rank("genus", upper = TRUE), "subfamily")
  expect_equal(get_lower_tax_rank("subfamily", upper = TRUE), "family")
  expect_equal(get_lower_tax_rank("family", upper = TRUE), "order")
  expect_equal(get_lower_tax_rank("order", upper = TRUE), "class")

})

test_that("textToDNAStringSet returns a DNAStringSet from both data.frame and character", {

  expect_true(is(textToDNAStringSet(as.character(">CODICE\nACGTACGTACGTACGT")), "DNAStringSet"))

  sequences <- data.frame(sourceID = "CODICE", markerCode = "MARKER", DNA_seq = "AGCTAGCTAGCTAGCT")

  expect_true(is(textToDNAStringSet(sequences), "DNAStringSet"))

})
