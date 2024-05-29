test_that("xml_Taxon_extract_nodes and XML_Taxon_extract return correct class and number of elements", {
  Sys.sleep(1/3)
  search <- ncbi_searcher("Thouarella", db = "taxonomy")
  Sys.sleep(1/3)
  fetch <- ncbi_xml_fetcher(search$web_history, db = "taxonomy", retstart = 0, rate = 20)
  nodes <- XML_Taxon_extract_nodes(fetch)

  expect_s3_class(nodes, "XMLNodeSet")
  expect_equal(length(nodes), 20)

  taxonomies <- extractTaxonomyTab(nodes)

  expect_equal(class(taxonomies), "data.frame")
  expect_equal(dim(taxonomies)[2], 9)
})

test_that("XML_extract_nodes returns XMLNodeSet class and XML_extract returns the desired different ouputs", {
  Sys.sleep(1/3)
  search <- ncbi_searcher("360011", db = "nucleotide")
  Sys.sleep(1/3)
  fetch <- ncbi_xml_fetcher(search$web_history, db = "nucleotide", retstart = 0, rate = 20)
  nodes_source <- XML_extract_nodes("source", fetch)
  nodes_CDS <- XML_extract_nodes("CDS", fetch)
  nodes_rRNA <- XML_extract_nodes("rRNA", fetch)

  expect_s3_class(nodes_source, "XMLNodeSet")
  expect_s3_class(nodes_CDS, "XMLNodeSet")
  expect_s3_class(nodes_rRNA, "XMLNodeSet")
  expect_equal(length(nodes_source), 20)

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

test_that("extractSelectionTab and extractRecordsTab give correct class and dimensions", {
  Sys.sleep(1/3)
  abatus_tax <- get_ncbi_taxonomy("Abatus", ask = FALSE)
  Sys.sleep(1/3)
  search <- ncbi_searcher("271674", db = "nucleotide")
  Sys.sleep(1/3)
  fetch <- ncbi_xml_fetcher(search$web_history, db = "nucleotide", retstart = 0, rate = 10)

  nodes_source <- XML_extract_nodes("source", fetch)
  nodes_feat <- XML_extract_nodes(c("CDS", "rRNA"), fetch)

  rec_tab <- extractRecordsTab(nodes_source, abatus_tax)
  sel_tab <- extractSelectionTab(nodes_feat, rec_tab$sourceID)

  expect_equal(class(rec_tab), "data.frame")
  expect_true(all(dim(rec_tab) == c(10, 29)))

  expect_equal(class(sel_tab), "data.frame")
  expect_true(all(dim(sel_tab) > 3))

})

