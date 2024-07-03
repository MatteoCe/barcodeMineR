#' Get taxonomic classification of the NCBI taxonomy database for a list of taxa
#' or taxids
#'
#' @param ids `character` A character string of species names or taxid from the
#'   NCBI.
#' @param api_rate `integer` The API rate with which to iterate each separate
#'   request. Must be a number between 3 and 10 which will translate in a rate
#'   of `1 / api_rate` seconds.
#' @param ask `logical` Should the function ask the user whether to filter the
#'   final output for taxonomic ranks. Default `TRUE`.
#'
#' @return `data.frame` A data.frame object with the searched taxonomic
#'   classification.
#'
#' @export
#' @importFrom rlang .data
#'
#' @description
#' The taxonomy functions are used to define exactly which records will be
#' retrieved by the download functions. More precisely, only records that
#' include the taxid/scientificName, retrieved using the taxonomy functions,
#' as the lowest taxonomic identification will be downloaded.
#'
#' For more details see the 'Searching taxonomy' vignette:
#' \code{vignette("Searching taxonomy", package = "barcodeMineR")}
#'
#' @seealso [get_bold_taxonomy()]
#'
#' @examples
#' get_ncbi_taxonomy("Polymastia invaginata")
#'
#' get_ncbi_taxonomy("283554")
#'
get_ncbi_taxonomy <- function(ids, api_rate = NULL, ask = TRUE) {

  # set the number of maximum ids to show in the esearch list object
  rate <- 200

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

  }, message = "Searching taxa on taxonomy database")

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
