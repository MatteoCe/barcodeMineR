#' Collapse query names that corresponds to the same taxonomy on the ncbi.
#'
#' @param taxonomy_df the ncbi taxonomy table as obtained in get_ncbi_taxonomy.
#'
#' @return a taxonomy table with duplicates collapsed.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom rlang .data
#'
#' @description
#' Sometimes it may happen that different taxonomic names correspond to only one name on the ncbi. This will collapse the rows with different query names into a single one with "name1|name2" in the field "QueryName".
#'
clean_taxonomy <- function(taxonomy_df) {

  # set binding visible variables
  queryName <- NULL

  taxids <- unique(taxonomy_df$taxid)

  purrr::map(taxids, function(id) {

    # extract the taxonomic table info corresponding to the species that is processed now
    single_tax <- taxonomy_df %>% dplyr::filter(., .data$taxid == id)

    # if multiple queryName species correspond to the same taxonomy in the NCBI tax database
    # then the input species if modified in order to contain a collapsed string
    # with both names.
    if (nrow(single_tax) > 1) {

      # collapse input species
      coll_input <- paste(single_tax$queryName, collapse = "|")

      # filter duplicate of rest columns
      rem_dup_tab <-  single_tax %>% dplyr::select(., -queryName) %>% unique()

      # bind input collapsed and table with removed duplicates
      single_tax <- single_tax %>% dplyr::mutate(., queryName = coll_input) %>%
        dplyr::select(queryName) %>% dplyr::bind_cols(., rem_dup_tab) %>% unique()

    }

    single_tax

  }) %>% do.call(rbind, .)


}

#' Search a taxon or taxid in the NCBI 'taxonomy' or 'nucleotide' databases
#'
#' @param id A taxon or taxid.
#' @param db The NCBI database that will be searched. Either 'nucleotide' or 'taxonomy'.
#' @param retmax the number of maximum ids to show in the esearch list object.
#' Defaults to 200.
#' @param default.filter logical, whether to filter the records excluding whole
#' genome shotgun sequences and transcribed shotgun assembly. Defaults to TRUE.
#' @param filter Additional text string/s that will be included in the query.
#' Works only with the 'nucleotide' option in the 'db' parameter.
#'
#' @return An esearch/list object, including a web_history object. The latter
#' will be used by the 'fetcher' functions to retrieve xml or fasta objects.
#'
#' @keywords internal
#' @noRd
#'
#' @description
#' This internal function will return an esearch/list object that will be used
#' by the internal 'fetcher' functions to retrieve the wanted data.
#' It can work only on a single 'id' object (a single string character) and for
#' this reason is usually run inside the functions reported in "requests_handlers.R"
#' like connection_handler() and ncbi_limit_handler().
#'
ncbi_searcher <- function(id, db = "taxonomy", retmax = 200, default.filter = TRUE, filter = NULL) {

  if (db == "taxonomy") {

    if (set_id_type(id) == TRUE) {
      term <- paste(id, "[ORGN]", sep = " ")

    } else {
      term <- paste0("txid", id, "[ORGN]")

    }

  } else if (db == "nucleotide") {

    term <- paste0(stringr::str_flatten(paste0("txid", id), collapse = "[ORGN] OR "), "[ORGN]")

    if (default.filter) {

      default.filters <- c('NOT wgs[Keyword]',
                           'NOT tsa[Keyword]',
                           'AND biomol_genomic[PROP]',
                           'AND (cds[Feature key] OR rrna[Feature key] OR "misc rna"[Feature key])')

      for (filt in default.filters) {

        term <- paste0("(", term, ") ", filt)

      }

    }

    if (!is.null(filter)) {

      for (filt in filter) {

        term <- paste0("(", term, ") ", filt)

      }

    }
  }

    search <- rentrez::entrez_search(
      db = db,
      term = term,
      retmax = retmax,
      use_history = TRUE
    )

    search
}

#' Extract ids from an esearch objects list.
#'
#' @param search_list the output of the function ncbi_searcher
#' @param db the database against which the ids were searched.
#' @param rate_xml the number of xml to be downloaded at a time.
#' @param api_rate the time rate adopted for each request sent to the NCBI.
#'
#' @return a list with character vectors of ncbi IDs to be used by the function postAndCheck for creating web_history objects.
#'
#' @keywords internal
#' @noRd
#' @importFrom rlang .data
#'
#' @description
#' Searching with esearch the ncbi database will give web_history objects and ids. However, the number of hits corresponding to that "search" may exceed the retrieved number of ids, based on the chosen retstart and retmax (the paramete rate_xml). This will run entrez_summary in order to retrieve all ids from each esearch object.
#'
ncbi_id_extract <- function(search_list, db = "nucleotide", rate_xml, api_rate = NULL) {

  if (purrr::some(search_list, is.list) &
      any(lapply(purrr::keep(search_list, is.list), function(search) {length(search)}) > 0)) {

    ids <- purrr::map(search_list, function(search) {

      if (is(search, "list")) {
        if (length(search) > 0) {

          splits <- web_history_splitter(search$count, 200)

          ncbi_limit_handler(splits, api_rate = api_rate, function(id) {

            split <- splits[id]

            parameter <- web_history_parameter(split, search$count, 200)

            sum <- rentrez::entrez_summary(db = db,
                                           web_history = search$web_history,
                                           retstart = parameter[[1]],
                                           retmax = parameter[[2]])

            # when the summary result is of class "esummary", it means that it refers
            # to only one gi/uid number, thus names(sum) will not report the gi/uid number but
            # the name of each element of the esummary record
            if (all(class(sum) == c("esummary", "list"))) {

              if (db == "taxonomy") {
                return(sum$uid)
              } else {
                return(sum$gi)
              }

            } else if (all(class(sum) == c("esummary_list", "list"))) {

              return(names(sum))

            } else {

              stop("issue with esummary report in ncbi_limit_handler")

            }

          }, message = "Retrieving additional ids") %>% do.call(c, .)

        } else {
          return(NULL)
        }

      } else if (is(search, "character")) {
        return(search)
      }

    }) %>% purrr::compact() %>% do.call(c, .)

  } else {

    ids <- do.call(c, purrr::keep(search_list, is.character))

  }

  if (db == "nucleotide") {

    ids <- unique(ids)

  }

  split(ids, ceiling(seq_along(ids) / rate_xml))

}

#' Fetch the taxonomic classification or the sequence metadata from the NCBI 'taxonomy' or 'nucleotide' databases
#'
#' @param searchWebHist A web_history object.
#' @param ids a character vector of accession/GI ids.
#' @param db The NCBI database from which to fetch data. Either 'nucleotide' or 'taxonomy'.
#' @param retstart The initial numeric position of the ids included in the web_history object and that will be fetched. For web_history objects with 1 ids only, 'NULL' is preferable to 0.
#' @param rate The final numeric position of the ids included in the web_history object and that will be fetched.
#'
#' @return Both 'taxonomy' and 'nucleotide' db searches will return an xml object of 'XMLInternalDocument/XMLAbstractDocument' class. This object should be processed with the functions at 'xml_helpers.R' to extract the wanted data.
#'
#' @keywords internal
#' @noRd
#'
#' @description
#' A web_history object is used to download xml objects from the NCBI database. The parameters retstart and rate are changed internally to download all ids from the web_history object, usually identifying 20 ids at a time.
#'
#'
#' @examples \dontrun{
#'   search <- ncbi_searcher("Thouarella", db = "taxonomy")
#'   fetch <- ncbi_xml_fetcher(search$web_history, retstart = 0, rate = 20)
#' }
ncbi_xml_fetcher <- function(searchWebHist = NULL, ids = NULL, db = "taxonomy", retstart, rate) {

  fetch <- rentrez::entrez_fetch(
    db = db,
    id = ids,
    web_history = searchWebHist,
    rettype = "xml",
    retmax = rate,
    retstart = retstart,
    parsed = TRUE
  )

  # get number of retrieved xmls
  if (db == "taxonomy") {
    numxml <- length(XML_Taxon_extract_nodes(fetch))

  } else if (db == "nucleotide") {
    numxml <- length(XML_extract_nodes("source", fetch))

  }

  # some records may have multiple sources (EU831279.1), thus now it stops only
  # when numxml is lower than rate
  if (numxml < rate) {
    stop("...")
  }

  fetch
}

#' Fetch fasta sequences from the NCBI 'nucleotide' database
#'
#' @param searchWebHist A web_history object.
#' @param ids a character vector of accession/GI ids.
#' @param retstart The initial numeric position of the ids included in the web_history object and that will be fetched. For web_history objects with 1 ids only, 'NULL' is preferable to 0.
#' @param rate The final numeric position of the ids included in the web_history object and that will be fetched.
#'
#' @keywords internal
#' @noRd
#'
#' @return A character string as fasta sequence/s.
#'
#' @description
#' Very similar to the ncbi_xml_fetcher, except that it retrieves fasta sequences.
#'
#' @examples \dontrun{
#'  search <- ncbi_searcher("360011", db = "nucleotide")
#'  fetch <- ncbi_fasta_fetcher(search$web_history, retstart = 0, rate = 20)
#' }
ncbi_fasta_fetcher <- function(searchWebHist = NULL, ids = NULL, retstart, rate) {

  seqFetch <- rentrez::entrez_fetch(
    db = "nucleotide",
    id = ids,
    web_history = searchWebHist,
    rettype = "fasta",
    retmax = rate,
    retstart = retstart
  )

  # get number of retrieved sequences
  numseq <- stringr::str_count(seqFetch, ">")

  if (numseq != rate) {
    stop("...")
  }

  seqFetch
}

#' Get a web_history object from a vector of accession numbers.
#'
#' @param names a vector of accession numbers.
#'
#' @return a list with a web_history object as first element and the corresponding accessions as second element.
#'
#' @description
#' Old version. Also checks that the web_history object contains the correct number of records.
#'
#'
#' @keywords internal
#' @noRd
#'
postAndCheck_OLD <- function(names) {

  post <- rentrez::entrez_post(db = "nucleotide", id = names)
  postSumm <- rentrez::entrez_summary(db = "nucleotide", web_history = post)

  if (length(postSumm) != length(names)) {# repeat if the summary
    stop("...")                           # does not contain same
  }                                       # numb of acc Num queried

  post

}

# old post

# create a list with web_history objects corresponding to groups of ncbi ids
#webenv_list <- ncbi_limit_handler(seq(1:length(ids_groups)), api_rate = api_rate, function(id) {

  # create web_history objects with the grouped accession numbers
  #rentrez::entrez_post(db = "nucleotide", id = ids_groups[[id]])

#})
