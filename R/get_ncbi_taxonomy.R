#' Get taxonomic classification of the NCBI taxonomy database for a list of taxa
#' or taxid.
#'
#' @param ids Character string of species names or taxid from the NCBI.
#' @param rate the number of maximum ids to show in the esearch list object.
#' @param api_rate The rate with which to perform the function on each element.
#' Must be a number between 3 and 10 which will translate in a rate of
#' '1 / api_rate' seconds.
#' @param ask Should the function ask the user whether to filter the final
#' output for taxonomic ranks. Default to TRUE.
#'
#' @return A data.frame object with the searched taxonomic classification.
#'
#' @export
#' @importFrom rlang .data
#'
#' @description
#' This function is the first step of the pipeline. It allows to download the
#' taxonomic classification of either taxonomic names or taxid from the NCBI.
#' A character vector of either type must be supplied to the ids parameter.
#' The api_rate parameter can be left NULL, as default, allowing the function
#' to automatically detect the best option related to the current future plan
#' and the NCBI http requests limit. In fact, if key is also left NULL, the
#' standard https requests limit will be adopted (3 requests in any 1 second
#' window), whereas if the NCBI API key is provided, up to 10 requests per
#' second will be performed. This speed up the recovery of data from the NCBI,
#' especially for time-consuming requests.
#'
#' The future plan framework allows to send requests every 1 / api_rate seconds,
#' thus ensuring that no more that 'api_rate' requests will be sent. However,
#' due to fluctuations in the internet connection, still more than that number
#' of requests will arrive to the NCBI server, causing errors. The function can
#' handle up to 5 consecutive errors per request, but too many errors might
#' block the whole process. Apart from reducing the number of taxonomic names
#' or taxid searched, the api_rate parameter can be modified in order to slow
#' down the requests sent per second. It overrides the automatic selection of
#' the optimal parameter (either 3 or 10) and accepts one decimal degree number
#' between 1.0 and 10.0, so if the internet connection is particularly bad, it
#' can be set to 2.5/2.0, for example, in order to slow the number of requests
#' per second and reduce the possibility of errors.
#'
#'
#' @examples \dontrun{
#' get_ncbi_taxonomy("Polymastia invaginata")
#'
#' get_ncbi_taxonomy("283554")
#'
#' }
get_ncbi_taxonomy <- function(ids, rate = 200, api_rate = NULL, ask = TRUE) {

  # check api limit can be adopted based on presence or not of an ncbi API key
  api_rate <- set_ncbi_rate(api_rate)

  # ensure that all ids are unique
  ids <- unique(ids)

  # search taxa or taxid in the NCBI
  search_list <- ncbi_limit_handler(ids, api_rate = api_rate, function(id) {

    search <- ncbi_searcher(ids[id], db = "taxonomy", retmax = rate)

    if (search$count <= rate) {

      search$ids

    } else {

      search

    }

  }, message = "Searching species on taxonomy database")

  # if no species name is found, stop now
  if (length(purrr::compact(search_list)) == 0) {

    stop(paste("No taxonomy for", ids, "was found in the NCBI database"))

  }

  # create a tab including the query species name and the taxid
  id_tab <- purrr::imap(search_list, \(x, y) if (length(x) > 0) {

    # define if the ids are taxonomic names or taxid
    if (set_id_type(ids) == TRUE) {
      query <- ids[y]

    } else {
      query <- NA

    }

    if ("esearch" %in% class(x)) {
      number <- x$count

    } else {
      number <- length(x)

    }

    data.frame(queryName = rep(query, number))

  }) %>% purrr::compact() %>% do.call(rbind, .)

  # get all ids from the previous search
  ids_groups <- ncbi_id_extract(search_list, db = "taxonomy", rate_xml = rate, api_rate = api_rate)

  # search the full taxonomy for the taxid
  taxonomies <- ncbi_limit_handler(ids_groups, api_rate = api_rate, function(id) {

    parameter <- web_history_parameter(0, length(ids_groups[[id]]), rate)

    taxonomy_xml <- ncbi_xml_fetcher(ids = ids_groups[[id]],
                                     db = "taxonomy",
                                     retstart = parameter[[1]],
                                     rate = parameter[[2]])

    xml_node <- XML_Taxon_extract_nodes(taxonomy_xml)

    extractTaxonomyTab(xml_node)

  }, message = "Searching taxonomic classification") %>% do.call(rbind, .) %>% dplyr::bind_cols(id_tab, .)

  if ((length(unique(taxonomies$rank)) > 1) && (ask == TRUE)) {
    filt_taxonomies <- ask_user(taxonomies, "rank", ask = ask)

    clean_taxonomies <- clean_taxonomy(filt_taxonomies)

  } else {
    clean_taxonomies <- clean_taxonomy(taxonomies)

  }

  clean_taxonomies

}

#' Collapse query names that corresponds to the same taxonomy on the ncbi.
#'
#' @param taxonomy_df the ncbi taxonomy table as obtained in get_ncbi_taxonomy.
#'
#' @return a taxonomy table with duplicates collapsed.
#'
#' @keywords internal
#' @importFrom rlang .data
#'
#' @description
#' Sometimes it may happen that different taxonomic names correspond to only one name on the ncbi. This will collapse the rows with different query names into a single one with "name1|name2" in the field "QueryName".
#'
clean_taxonomy <- function(taxonomy_df) {

  taxids <- unique(sort(taxonomy_df$taxid))

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
      rem_dup_tab <-  single_tax %>% dplyr::select(., -.data$queryName) %>% unique()

      # bind input collapsed and table with removed duplicates
      single_tax <- single_tax %>% dplyr::mutate(., queryName = coll_input) %>%
        dplyr::select(.data$queryName) %>% dplyr::bind_cols(., rem_dup_tab) %>% unique()

    }

    single_tax

  }) %>% do.call(rbind, .)


}
