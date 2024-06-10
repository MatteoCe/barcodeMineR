## code to prepare `example_record` dataset goes here

# record corresponding to fake "Dissostichus mawsoni" sequence as example for
# loadBarcodeOre function
example_record <- structure(list(sourceID = "SEQ_01", source = "ACRONYM", phylum = "Chordata",
                                 class = "Actinopteri", order = "Perciformes", family = "Nototheniidae",
                                 genus = "Dissostichus", species = "Dissostichus mawsoni",
                                 lengthSource = 658L, sampleID = "SAMPLE01", identified_by = "John Smith",
                                 taxNotes = NA_character_, db_xref = NA_character_, NCBI_ID = NA_character_,
                                 institutionStoring = "INSTITUTION", collected_by = "John Smith",
                                 collection_date = "23-Aug-1992", altitude = NA_character_,
                                 depth = 300, country = "Southern Ocean", lat = -74.6789,
                                 lon = 164.2315, directionPrimers = "F|R", PCR_primers = "newPrimF|newPrimR",
                                 note = NA_character_), row.names = c(NA, -1L), class = c("tbl_df",
                                                                                          "tbl", "data.frame"))

usethis::use_data(example_record, overwrite = TRUE)
