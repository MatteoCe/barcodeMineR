test_that("loadBarcodeOre works on loaded objects", {

  bo <- loadBarcodeOre(example_record, example_sequence)
  expect_true(is(bo, "data.frame"))

})

test_that("loadBarcodeOre works on tsv and fasta files", {

  write.table(example_record, file = paste(tempdir(), "temp_example.tsv", sep = "/"), sep = "\t")
  Biostrings::writeXStringSet(example_sequence, paste(tempdir(), "temp_example.fa", sep = "/"))

  path_rec <- file.path(tempdir(), "temp_example.tsv")
  path_seq <- file.path(tempdir(), "temp_example.fa")

  bo <- loadBarcodeOre(path_rec, path_seq)

  expect_true(is(bo, "data.frame"))

})
