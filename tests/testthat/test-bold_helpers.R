# CONTEXT

tax_table1 <- structure(list(queryName = c("Clio", "Clio pyramidata", "Clio pyramidata antarctica",
"Limacina helicina antarctica", "Maldane sarsi"), taxid = c(193218L,
202265L, 655201L, 193227L, 204051L), taxon = c("Clio", "Clio pyramidata",
"Clio pyramidata antarctica", "Limacina helicina antarctica",
"Maldane sarsi"), rank = c("genus", "species", "subspecies",
"subspecies", "species")), row.names = c(NA, 5L), errors = structure(list(), names = character(0)), params = list(
    fuzzy = FALSE, tax_division = NULL, tax_rank = NULL), class = "data.frame")

bold_count1 <- bold_record_counter(tax_table1, 0.06)
tax_table2 <- dplyr::full_join(tax_table1, bold_count1, by = "taxon")
bold_group1 <- bold_record_grouper(tax_table2, 100)
Sys.sleep(1)
bold_fetch1 <- bold_fetcher(bold_group1[[2]], tax_table2)
# tests cleanMined
tax_table3 <- structure(list(queryName = "Kali indica", taxid = 59484L, taxon = "Kali indica",
    rank = "species", records = 9L), errors = structure(list(), names = character(0)), params = list(
    fuzzy = FALSE, tax_division = NULL, tax_rank = NULL), row.names = 202L, class = "data.frame")
bold_group2 <- bold_record_grouper(tax_table3, 100)
Sys.sleep(1)
bold_fetch2 <- ncbi_limit_handler(bold_group2, api_rate = 1, function(id) {
  bold_fetcher(bold_group2[[id]], tax_table3)
}, message = "Downloading BOLD records", seed = NULL) %>% purrr::compact() %>% do.call(rbind, .)
bold_fetch3 <- extractRecordsTabBOLD(bold_fetch2, tax_table3)
sequences1 <- textToDNAStringSet(bold_fetch3)
#

# TESTS

test_that("bold_record_counter removes children taxonomies already present in taxonomic table", {

  expect_s3_class(bold_count1, "data.frame")
  expect_equal(nrow(bold_count1), 3)

})

test_that("bold_record_grouper groups species based on rate number", {

    expect_true(is(bold_group1, "list"))

    for (rate in c(100, 10)) {
        bold_group2 <- bold_record_grouper(tax_table2, rate)
        expect_equal(length(bold_group2), if (rate == 10) {3} else {2})
    }

})

test_that("bold_fetcher gives data.frame of correct dimensions", {

    expect_true(nrow(bold_fetch1) < sum(tax_table2[tax_table2$taxon %in% bold_group1[[2]], ]$records))
    expect_equal(ncol(bold_fetch1), 80)
    expect_s3_class(bold_fetch1, "data.frame")

})

test_that("extractRecordsTabBOLD gives data.frame of correct dimensions", {

    bold_extract1 <- extractRecordsTabBOLD(bold_fetch1, tax_table2)

    expect_s3_class(bold_extract1, "data.frame")
    expect_true(nrow(bold_extract1) < sum(tax_table2[tax_table2$taxon %in% bold_group1[[2]], ]$records))
    expect_equal(ncol(bold_extract1), 29)

})

test_that("cleanMinedFromGenBankRecords behaves correctly", {

    # lacking examples that refer to records Mined from GenBank more than one
    # time for now it only tests normal records
    accn1 <- paste(bold_fetch3[bold_fetch3$species %in% "Kali indica", ]$sourceID,
                   bold_fetch3[bold_fetch3$species %in% "Kali indica", ]$markerCode,
                   sep = "|")
    sequences1 <- textToDNAStringSet(bold_fetch3)
    accn2 <- cleanMinedFromGenBankRecords(bold_fetch3, sequences1, accn1)

    expect_equal(accn1, accn2)

    # fake record and sequence
    record_tab1 <- structure(list(source = c("BOLD", "BOLD"), sourceID = c("FNZB069-08",
"GBMTG5038-16"), DNA_seq = c("CCTGTATCTGGTATTTGGTGCATGAGCTGGTATAGTAGGTACAGCCTTAAGCCTTCTAATCCGGGCTGAATTAAGCCAACCTGGTGCCCTACTTGGGGATGACCAAATTTACAATGTTATCGTTACAGCCCACGCCTTTGTAATGATTTTCTTTATAGTAATACCAATCATGATTGGAGGTTTCGGTAATTGACTCATCCCATTAATGATTGGGGCCCCCGATATAGCATTCCCTCGAATAAACAACATGAGCTTTTGACTTTTACCTCCTTCCTTCCTGCTTCTTTTAGCATCTTCTGGTGTTGAAGCTGGAGCCGGGACTGGGTGAACAGTTTATCCCCCTTTAGCTGGTAATCTAGCCCACGCCGGAGCATCAGTAGATTTAACTATCTTTTCACTACATCTGGCAGGGGTTTCATCTATCCTTGGGGCTATCAATTTTATCACAACCATTATT",
"CCTGTATCTGGTATTTGGTGCATGAGCTGGTATAGTAGGTACAGCCTTAAGCCTTCTAATCCGGGCTGAATTAAGCCAACCTGGTGCCCTACTTGGGGATGACCAAATTTACAATGTTATCGTTACAGCCCACGCCTTTGTAATGATTTTCTTTATAGTAATACCAATCATGATTGGAGGTTTCGGTAATTGACTCATCC"
), markerCode = c("COI-5P", "COI-5P"), lengthGene = c(457, 200
), sampleID = c("FJ969274", "FJ969274"), QueryName = c("Ammothea carolinensis",
"Ammothea carolinensis"), phylum = c("Arthropoda", "Arthropoda"
), class = c("Pycnogonida", "Pycnogonida"), order = c("Pantopoda",
"Pantopoda"), family = c("Ammotheidae", "Ammotheidae"), genus = c("Ammothea",
"Ammothea"), species = c("Ammothea carolinensis", "Ammothea carolinensis"
), identified_by = c("Andrew Stewart", NA), taxNotes = c(NA,
NA), db_xref = c("taxon:59484", "taxon:59484"), NCBI_ID = c("FJ969274",
NA), institutionStoring = c("Museum of New Zealand Te Papa Tongarewa",
NA), collected_by = c("RV Tangaroa", NA), collection_date = c(NA,
NA), lat = c(-66.832, NA), lon = c(171.262, NA), altitude = c(NA,
NA), depth = c(50, NA), country = c("Southern Ocean", NA), directionPrimers = c("R|F",
NA), lengthSource = c(457, 457), PCR_primers = c("M13R|M13F",
NA), note = c("Mined from GenBank, NCBI", "Mined from GenBank, NCBI"
)), row.names = 1:2, class = "data.frame")

    sequences2 <- textToDNAStringSet(record_tab1)
    sequences2[1] <- Biostrings::subseq(sequences2[1], start = 1, end = 457)
    sequences2[2] <- Biostrings::subseq(sequences2[2], start = 1, end = 200)
    accn3 <- paste(record_tab1$sourceID, record_tab1$markerCode, sep = "|")
    clean_accn <- cleanMinedFromGenBankRecords(record_tab1, sequences2, accn3)

    expect_equal("FNZB069-08|COI-5P", clean_accn)

})
