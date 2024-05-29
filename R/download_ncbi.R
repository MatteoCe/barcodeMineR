#' Download records from the NCBI
#'
#' @param ncbi_tax a data.frame, as returned from the get_ncbi_taxonomy function.
#' @param ncbi_ids a character vector including ncbi accession numbers.
#' @param rate_xml the number of xml to be downloaded at a time. Defaults to
#' 200, but can be lowered for bad internet connections.
#' @param rate_fasta the number of fasta sequences to be downloaded at a time.
#' Same as rate_xml, but defaults to 100. Many fasta can correspond to
#' mitogenomes and chromosomes, which may lead to an error if downloaded in
#' great numbers. Can be increased if needed.
#' @param default.filter logical, whether to filter the records excluding whole
#' genome shotgun sequences and transcribed shotgun assembly. Defaults to TRUE.
#' @param filter an additional query filter in the form of a/multiple string/s
#' to add to every searched taxid. This will allow any user to specifically
#' filter every search with a custom query. Multiple strings should be provided
#' in the form of a character vector of single query filters (see description
#' for details).
#' @param api_rate the time rate adopted for each request sent to the NCBI.
#' Defaults to NULL, which will let the function check whether the key parameter
#' is also NULL or not. In the former case, the api rate will be set to a
#' request every 1/3 seconds, otherwise every 1/10 seconds. Works only in a
#' multisession future plan. Check the "Description" for more information.
#' @param ask defaults to TRUE. If disabled the selection process will be
#' skipped and all feature keys from each accession number will be retrieved.
#' @param prefix defaults to NULL. Character string that will be used to create
#' numbered custom ids for each record in ascending order.
#'
#' @return a refdb tibble data.frame, including the DNA sequence as a field.
#'
#' @importClassesFrom Biostrings DNAStringSet
#' @importFrom methods as
#'
#' @export
#'
#' @description
#' A short description...
download_ncbi <- function(ncbi_tax = NULL, ncbi_ids = NULL, rate_xml = 200, rate_fasta = 100, default.filter = TRUE, filter = NULL, api_rate = NULL, ask = TRUE, prefix = NULL) {

  ### Search (step) each taxid in the nucleotide database with entrez_search

  # check api limit can be adopted based on presence or not of an ncbi API key
  api_rate <- set_ncbi_rate(api_rate)

  # if no taxonomic table is provided, skip the first steps of the function
  # and use ncbi_ids to search the records
  if (!is.null(ncbi_tax)) {

    # get list of txids
    taxids <- unique(sort(ncbi_tax$taxid)) %>% split(., ceiling(seq_along(.) / 7))

    # search taxa or taxid in the NCBI
    search_list <- ncbi_limit_handler(taxids, api_rate = api_rate, function(id) {

      search <- ncbi_searcher(taxids[[id]], db = "nucleotide", retmax = rate_xml, default.filter = default.filter, filter = filter)

      if (search$count <= rate_xml) {

        search$ids

      } else {

        search

      }

    }, message = "Searching taxa on nucleotide database")

    # get all ids from the previous search
    ids_groups <- ncbi_id_extract(search_list, db = "nucleotide", rate_xml = rate_xml, api_rate = api_rate)

  } else {

    ids_groups <- split(ncbi_ids, ceiling(seq_along(ncbi_ids) / rate_xml))

  }

  # Download xml for each web_history
  xml_tables <- ncbi_limit_handler(ids_groups, api_rate = api_rate, function(id) {

    parameter <- web_history_parameter(0, length(ids_groups[[id]]), rate_xml)

    # download the xml for that specific web_history
    fetch_xml <- ncbi_xml_fetcher(ids = ids_groups[[id]],
                                  db = "nucleotide",
                                  retstart = parameter[[1]],
                                  rate = parameter[[2]])

    # create list to store selection and records tabs
    storing_list <- list()

    # append xml to new custom root
    fetch_xml_doc <- XML_root(fetch_xml)

    # extract all information as a selection tab and records_tab

    # extract sources from the whole XML object
    featureSource_nodes <- XML_extract_nodes("source", fetch_xml_doc)

    # create the records table
    records_tab <- extractRecordsTab(featureSource_nodes, ncbi_tax)

    # extract all CDS and rRNA features from the whole XML object
    feature_nodes <- XML_extract_nodes(c("CDS", "rRNA"), fetch_xml_doc)

    # create selection table and records table
    selection_tab <- extractSelectionTab(feature_nodes, records_tab$sourceID)

    # merge both types of tabs as a single list
    storing_list[[1]] <- selection_tab
    storing_list[[2]] <- records_tab

    storing_list

  }, message = "Processing records information from nucleotide database")

  # extract selection tab from list of tables
  selection_tab <- lapply(xml_tables, function(tab) {

    tab[[1]]

  }) %>% do.call(rbind, .)

  # extract records tab from list of tables and update formats as required in barcodeOre objects
  records_tab <- lapply(xml_tables, function(tab) {

    tab[[2]]

  }) %>% do.call(rbind, .)

  ### Select (step)

  names <- select_accessions(selection_tab, ask = ask) %>% unique() %>% .[!is.na(.)]

  ### Download sequences (step) from the NCBI based on the vector of accession numbers selected earlier

  # split the vector of accession numbers to groups based on rate chosen
  accn <- unique(stringr::str_remove_all(names, "\\|.*"))

  # get all ids from the previous search
  seq_groups <- split(accn, ceiling(seq_along(accn) / rate_fasta))

  # get web_history for each group of "rate" number of accession numbers
  fasta_list <- ncbi_limit_handler(seq_groups, api_rate = api_rate, function(id) {

    parameter <- web_history_parameter(0, length(seq_groups[[id]]), rate_fasta)

    # download fasta for each group of accession numbers
    fastaResults <- ncbi_fasta_fetcher(ids = seq_groups[[id]],
                                       retstart = parameter[[1]],
                                       rate = parameter[[2]])

    # change format to fasta and then to DNAStringSet
    dnaStringFormat <- textToDNAStringSet(fastaResults)

    # xxx
    names_group <- purrr::keep(stringr::str_split(names, "\\|"), function(name) {

      name[1] %in% seq_groups[[id]]

    }) %>%

      purrr::map_chr(., function(name) {

        stringr::str_c(name[1], name[2], sep = "|")

      }) %>%

      unique()

    # apply buildSequences to all chosen accession numbers|gene_name
    sequences_sel <- lapply(names_group, buildSequences, DNAString = dnaStringFormat, selection_tab = selection_tab) %>% purrr::compact()

    # merge all into a single object
    sequences_total <- do.call(c, sequences_sel)

    sequences_total

  }, message = "Processing fasta from nucleotide database")

  ### Build records table and barcodeOre object (step)

  # merge all DNAstringSet into a single object
  sequences_total <- do.call(c, fasta_list)

  # apply buildRecords to all chosen accession numbers|gene_name, using the new
  # total sequences as they contain the correct region of the chosen CDS or rRNA
  p <- progressr::progressor(steps = length(names), on_exit = TRUE)

  records_total <- future.apply::future_lapply(names, function(name) {

    p(message = sprintf("Final processing of records and sequences"))

    buildRecord(name, records = records_tab, sequences = sequences_total)

  }) %>% purrr::compact() %>% do.call(rbind, .)

  # create refdb data.frame
  buildBarcodeOre(records_total, prefix)

}
