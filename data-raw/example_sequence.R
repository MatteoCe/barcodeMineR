## code to prepare `example_sequence` dataset goes here

# fake sequence used as example for the loadBarcodeOre function. To use in
# combination with the 'example_record' object, including the corresponding
# (fake) metadata.
example_sequence <-
  Biostrings::DNAStringSet("AAACTCAAAGAGAAATTTAAAAAACTCAAAGAGAAATTTAAAAAACTCAAA
                           GAGAAATTTAAAAAACTCAAAGAGAAATTTAAAAAACTCAAAGAGAAATTTA
                           AAAAACTCAAAGAGAAATTTAAAAAACTCAAAGAGAAATTTAAAAAACTCAA
                           AGAGAAATTTAAAAAACTCAAAGAGAAATTTAAAAAACTCAAAGAGAAATTT
                           AAAAAACTCAAAGAGAAATTTAAAAAACTCAAAGAGAAATTTAAAAAACTCA
                           AAGAGAAATTTAAAAAACTCAAAGAGAAATTTAAAAAACTCAAAGAGAAATT
                           TAAAAAACTCAAAGAGAAATTTAAAAAACTCAAAGAGAAATTTAAAAAACTC
                           AAAGAGAAATTTAAAAAACTCAAAGAGAAATTTAAAAAACTCAAAGAGAAAT
                           TTAAAAAACTCAAAGAGAAATTTAAAAAACTCAAAGAGAAATTTAAAAAACT
                           CAAAGAGAAATTTAAAAAACTCAAAGAGAAATTTAAAAAACTCAAAGAGAAA
                           TTTAAAAAACTCAAAGAGAAATTTAAAAAACTCAAAGAGAAATTTAAAAAAC
                           TCAAAGAGAAATTTAAAAAACTCAAAGAGAAATTTAAAAAACTCAAAGAGAA
                           ATTTAAAAAACTCAAAGAGAAATTTAAATGATGAC")

example_sequence@ranges@NAMES <- "SEQ_01|COI"

usethis::use_data(example_sequence, overwrite = TRUE)
