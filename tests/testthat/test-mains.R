# CONTEXT

tax_ncbi <- get_ncbi_taxonomy("Clio pyramidata", ask = FALSE)
tax_bold <- get_bold_taxonomy("Clio pyramidata", descend = TRUE, ask = FALSE)

rec_ncbi <- download_ncbi(tax_ncbi, ask = FALSE)
rec_bold <- download_bold(tax_bold, ask = FALSE)
rec_custom1 <- loadBarcodeOre(example_record, example_sequence)

dim1 <- nrow(rec_ncbi)
dim2 <- nrow(rec_bold)
dim3 <- nrow(rec_custom1)

merged1 <- mergeBarcodeOres(rec_ncbi, rec_bold, resolve.conflicts = FALSE)
merged2 <- mergeBarcodeOres(rec_ncbi, rec_bold, resolve.conflicts = TRUE)
merged3 <- mergeBarcodeOres(rec_ncbi, rec_bold, rec_custom1, resolve.conflicts = FALSE)

# TESTS

test_that("loadBarcodeOre works on loaded objects", {

  expect_true(is(rec_custom1, "data.frame"))
  expect_equal(1, nrow(rec_custom1))

})

test_that("loadBarcodeOre works on tsv and fasta files", {

  write.table(example_record, file = paste(tempdir(), "temp_example.tsv", sep = "/"), sep = "\t")
  Biostrings::writeXStringSet(example_sequence, paste(tempdir(), "temp_example.fa", sep = "/"))

  path_rec <- file.path(tempdir(), "temp_example.tsv")
  path_seq <- file.path(tempdir(), "temp_example.fa")

  rec_custom2 <- loadBarcodeOre(path_rec, path_seq)

  expect_true(is(rec_custom2, "data.frame"))

})

test_that("mergeBarcodeOres works", {

  expect_true(is(merged1, "data.frame"))
  expect_true(is(merged2, "data.frame"))
  expect_true(is(merged3, "data.frame"))
  expect_equal(nrow(merged1), sum(dim1, dim2))
  expect_equal(nrow(merged3), sum(dim1, dim2, dim3))

})

test_that("plot functions work", {

  x1 <- plot_length(rec_ncbi)
  x2 <- plot_primers(rec_ncbi)
  x3 <- plot_length(rec_bold)
  x4 <- plot_primers(rec_bold)
  x5 <- plot_primers(rec_ncbi, select = "Clio", size_range = c(1, 2))
  x6 <- plot_length(rec_ncbi, select = "Clio", level = "genus")

  expect_true(is(x1, "ggplot"))
  expect_true(is(x2, "ggplot"))
  expect_true(is(x3, "ggplot"))
  expect_true(is(x4, "ggplot"))
  expect_true(is(x5, "ggplot"))
  expect_true(is(x6, "ggplot"))

})
