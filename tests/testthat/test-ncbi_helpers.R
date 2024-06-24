# CONTEXT
Sys.sleep(1/2.8)
search1 <- ncbi_searcher("Thouarella", db = "taxonomy")
search1_list <- list("1" = search1)
rate_xml <- 20
id_list1 <- ncbi_id_extract(search1_list, rate_xml = rate_xml, api_rate = 2.8)
Sys.sleep(1/2.8)
fetch1 <- ncbi_xml_fetcher(search1$web_history, db = "taxonomy", retstart = 0, rate = rate_xml)
# below: offline but used in two tests
taxon_nodes <- XML_Taxon_extract_nodes(fetch1)
#
Sys.sleep(1/2.8)
search2 <- ncbi_searcher("283554", db = "nucleotide")
if (search2$count >= 20) {
    rate <- 20
} else {
    rate <- search2$count
}
Sys.sleep(1/2.8)
fetch2 <- ncbi_xml_fetcher(search2$web_history, db = "nucleotide", retstart = 0, rate = rate)
# below: offline but used in two tests
nodes_source <- XML_extract_nodes("source", fetch2)
nodes_CDS <- XML_extract_nodes("CDS", fetch2)
nodes_rRNA <- XML_extract_nodes("rRNA", fetch2)
nodes_feats <- XML_extract_nodes(c("CDS", "rRNA"), fetch2)
#
Sys.sleep(1/2.8)
fetch3 <- ncbi_fasta_fetcher(search2$web_history, retstart = 0, rate = rate)
# below: all offline, provide taxonomy tap manually
tax_tab <- structure(list(queryName = "Polymastia invaginata", 
                          taxid = "283554", 
                          rank = "species", 
                          scientificName = "Polymastia invaginata", 
                          phylum = "Porifera", 
                          class = "Demospongiae", 
                          order = "Polymastiida",
                          family = "Polymastiidae", 
                          genus = "Polymastia", 
                          species = "Polymastia invaginata"), 
                     row.names = c(NA, -1L), 
                     class = "data.frame")

rec_tab <- extractRecordsTab(nodes_source, tax_tab)
sel_tab <- extractSelectionTab(nodes_feats, rec_tab$sourceID)
names <- select_accessions(sel_tab, ask = FALSE) %>% unique() %>% .[!is.na(.)]
#
# below: offline
sequences <- textToDNAStringSet(fetch3)

accn <- unique(stringr::str_remove_all(names, "\\|.*"))
seq_groups <- split(accn, ceiling(seq_along(accn) / rate))

names_group <- purrr::keep(stringr::str_split(names, "\\|"), function(name) {
    name[1] %in% seq_groups[[1]]
  }) %>%
    purrr::map_chr(., function(name) {
      stringr::str_c(name[1], name[2], sep = "|")
    }) %>%
    unique()

seq_tab_tot <- lapply(names_group, 
                        buildSequences, 
                        DNAString = sequences, 
                        selection_tab = sel_tab) %>% purrr::compact() %>% do.call(c, .)
rec_tab_tot <- lapply(names, function(name) {
    buildRecord(name, records = rec_tab, sequences = seq_tab_tot)
  }) %>% purrr::compact() %>% do.call(rbind, .)
final_tab <- buildBarcodeOre(rec_tab_tot)
final_tab_prefix <- buildBarcodeOre(rec_tab_tot, prefix = "RECORD")
#

# TESTS
test_that("clean_taxonomy gives a data.frame with collapsed duplicated names", {

    # offline function, using made up data:

    clean_taxonomy_test1 <- structure(list(queryName = c("Laevilacunaria antarctica", "Laevilitorina antarctica"
), taxid = c("2841070", "2841070"), rank = c("species", "species"
), scientificName = c("Laevilacunaria antarctica", "Laevilacunaria antarctica"
), phylum = c("Mollusca", "Mollusca"), class = c("Gastropoda", 
"Gastropoda"), order = c("Littorinimorpha", "Littorinimorpha"
), family = c("Littorinidae", "Littorinidae"), genus = c("Laevilacunaria", 
"Laevilacunaria"), species = c("Laevilacunaria antarctica", "Laevilacunaria antarctica"
)), row.names = c(NA, -2L), class = "data.frame")
  
  clean_taxonomy_test1_result <- structure(list(queryName = "Laevilacunaria antarctica|Laevilitorina antarctica", 
                                              taxid = "2841070", rank = "species", scientificName = "Laevilacunaria antarctica", 
                                              phylum = "Mollusca", class = "Gastropoda", order = "Littorinimorpha", 
                                              family = "Littorinidae", genus = "Laevilacunaria", species = "Laevilacunaria antarctica"), row.names = c(NA, 
                                                                                                                                                       -1L), class = "data.frame")

  result <- clean_taxonomy(clean_taxonomy_test1)

  expect_equal(result, clean_taxonomy_test1_result)

})

test_that("ncbi_searcher returns esearch/list object", {

  expect_s3_class(search1, c("esearch", "list"))
  expect_s3_class(search2, c("esearch", "list"))

})

test_that("ncbi_xml_fetcher returns XMLInternalDocument/XMLAbstractDocument object", {

  expect_s3_class(fetch1, c("XMLInternalDocument", "XMLAbstractDocument"))
  expect_s3_class(fetch2, c("XMLInternalDocument", "XMLAbstractDocument"))

})

test_that("ncbi_fasta_fetcher returns character object", {

  expect_true(is(fetch3, "character"))

})

test_that("ncbi_id_extract returns list of character vectors", {

  expect_true(is(id_list1, "list"))
  expect_true(is(id_list1[[1]], "character"))

})

test_that("ncbi_id_extract returns of count/rate_xml number of elements", {

  expect_equal(length(id_list1), search1_list[[1]]$count/rate_xml)

})

test_that("xml_Taxon_extract_nodes extracts a list of rate_xml elements", {

  expect_s3_class(taxon_nodes, "XMLNodeSet")
  expect_equal(length(taxon_nodes), rate_xml)

})

test_that("XML_extract_nodes gives a XMLNodeSet class for different outputs and rate number of sources", {

  expect_s3_class(nodes_source, "XMLNodeSet")
  expect_s3_class(nodes_CDS, "XMLNodeSet")
  expect_s3_class(nodes_rRNA, "XMLNodeSet")
  expect_equal(length(nodes_source), rate)

})

test_that("XML_extract gives different outputs based on argument", {

  # offline
  expect_true(all(sapply(nodes_source, function(x) {
    accn <- XML_extract(x, "accession", quals = NULL)
    stringr::str_detect(accn, "^[[:alnum:]|\\_]*\\.[0-9]$")})))

  expect_true(all(sapply(nodes_source, function(x) {
    location <- XML_extract(x, "location", quals = NULL)
    stringr::str_detect(location, "[\\<0-9]+\\.\\.[\\>0-9]+")})))

  expect_true(all(sapply(nodes_source, function(x) {
    from <- XML_extract(x, "location", quals = "from")
    to <- XML_extract(x, "location", quals = "to")
    stringr::str_detect(from, "^[:alnum:]*$")
    stringr::str_detect(to, "^[:alnum:]*$")})))

  expect_true(all(sapply(nodes_CDS, function(x) {
    accn <- XML_extract(x, "qualifier", quals = "gene")
    stringr::str_detect(accn, "^[:alnum:]*$")})))

})

test_that("XML_root gives correct class XMLInternalElementNode/XMLInternalNode/XMLAbstractNode", {

  xml_root <- XML_root(fetch2)

  expect_s3_class(xml_root, c("XMLInternalElementNode", "XMLInternalNode", "XMLAbstractNode"))

})

test_that("extractTaxonomyTab gives a data.frame of 9 columns and rate_xml rows", {

  # offline
  taxon_nodes_tab <- extractTaxonomyTab(taxon_nodes)

  expect_s3_class(taxon_nodes_tab, "data.frame")
  expect_equal(dim(taxon_nodes_tab)[2], 9)
  expect_equal(dim(taxon_nodes_tab)[1], rate_xml)

})

test_that("extractRecordsTab give data.frame and 'rate' records", {

  expect_equal(class(rec_tab), "data.frame")
  expect_true(all(dim(rec_tab) == c(rate, 29)))

})

test_that("extractSelectionTab give data.frame and correct dimension", {

  expect_equal(class(sel_tab), "data.frame")
  expect_equal(dim(sel_tab)[1], length(nodes_feats))
  expect_equal(dim(sel_tab)[2], 4)

})

test_that("textToDNAStringSet returns a DNAStringSet from both data.frame and character", {

  sequences <- data.frame(sourceID = "CODICE", markerCode = "MARKER", DNA_seq = "AGCTAGCTAGCTAGCT")

  expect_true(is(textToDNAStringSet(fetch3), "DNAStringSet"))
  expect_true(is(textToDNAStringSet(sequences), "DNAStringSet"))

})

test_that("select_accession returns all values if parameter ask is disabled", {

  expect_equal(length(names), dim(rec_tab)[1])

})

test_that("buildSequences gives DNAStringSet with correct names", {

  expect_true(is(seq_tab_tot, "DNAStringSet"))
  expect_equal(seq_tab_tot@ranges@NAMES, names_group)

})

test_that("buildRecords gives data.frame with correct dimensions", {

  expect_true(is(rec_tab_tot, "data.frame"))
  expect_equal(nrow(rec_tab_tot), rate)
  expect_equal(ncol(rec_tab_tot), 29)

})

test_that("buildBarcodeOre gives data.frame with correct dimensions", {

  new_col1 <- colnames(final_tab)[!(colnames(final_tab) %in% colnames(rec_tab_tot))]
  new_col2 <- colnames(final_tab_prefix)[!(colnames(final_tab_prefix) %in% colnames(rec_tab_tot))]

  expect_equal(new_col1, "recordID")
  expect_equal(new_col2, "recordID")

  expect_true(is(final_tab, "data.frame"))
  expect_equal(nrow(final_tab), rate)
  expect_equal(ncol(final_tab), 30)

  expect_true(is(final_tab_prefix, "data.frame"))
  expect_equal(nrow(final_tab_prefix), rate)
  expect_equal(ncol(final_tab_prefix), 30)

})
