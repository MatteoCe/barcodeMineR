test_that("searcher in taxonomy returns esearch/list class", {
  Sys.sleep(1/3)
  expect_s3_class(ncbi_searcher("Thouarella", db = "taxonomy"), c("esearch", "list"))
})

test_that("searcher in nucleotide returns esearch/list class", {
  Sys.sleep(1/3)
  expect_s3_class(ncbi_searcher("283554", db = "nucleotide"), c("esearch", "list"))
})

test_that("fetcher in taxonomy returns XMLInternalDocument/XMLAbstractDocument class", {
  Sys.sleep(1/3)
  search <- ncbi_searcher("Thouarella", db = "taxonomy")
  Sys.sleep(1/3)
  fetch <- ncbi_xml_fetcher(search$web_history, db = "taxonomy", retstart = 0, rate = 20)
  expect_s3_class(fetch, c("XMLInternalDocument", "XMLAbstractDocument"))
})

test_that("fetcher in nucleotide returns XMLInternalDocument/XMLAbstractDocument class", {
  Sys.sleep(1/3)
  search <- ncbi_searcher("283554", db = "nucleotide")
  Sys.sleep(1/3)
  fetch <- ncbi_xml_fetcher(search$web_history, db = "nucleotide", retstart = 0, rate = search$count)
  expect_s3_class(fetch, c("XMLInternalDocument", "XMLAbstractDocument"))
})

test_that("fasta fetcher returns character class", {
  Sys.sleep(1/3)
  search <- ncbi_searcher("360011", db = "nucleotide")
  Sys.sleep(1/3)
  fetch <- ncbi_fasta_fetcher(search$web_history, retstart = 0, rate = 20)
  expect_equal(class(fetch), "character")
})

test_that("ncbi_id_extract returns list of character vectors numbered as total count / rate_xml", {

  search_list <- list()

  search_list[[1]] <- ncbi_searcher("Thouarella", db = "taxonomy")
  id_list <- ncbi_id_extract(search_list, rate_xml = 20, api_rate = 3)

  expect_equal(length(id_list), search_list[[1]]$count/20)
  expect_true(is(id_list, "list"))
  expect_true(is(id_list[[1]], "character"))

})
