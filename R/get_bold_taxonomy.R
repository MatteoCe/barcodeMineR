#' Get taxonomic classification of the BOLD taxonomy database for a list of taxa.
#'
#' @param ids character vector of species names
#' @param api_rate the rate with which to perform the function on each element.
#' For the BOLD api it can be higher, thus it will be defined by the number of
#' workers available in the chosen future framework, unless specified otherwise
#' @param descend logical. If set to TRUE (default), use the taxize package to
#' retrieve lower level taxonomies.
#' @param group grouping factor for ids names. Defaults to 10.
#' @param ask should the function ask the user whether to filter the final
#' output for taxonomic ranks. Default to TRUE.
#'
#' @return a data.frame with taxid from the taxonomy database of the BOLD
#' website and the corresponding species name(s) retrieved
#'
#'
#' @export
#' @importFrom methods is
#' @importFrom rlang .data
#'
#' @description
#' Due to the functioning of the bold package, random number generation might
#' set warning messages when this function is performed in multisession plan.
#' For this reason, the seed parameter in the future::future function called
#' internally is set to NULL, in order to suppress these warnings.
#'
#'
#' @examples \dontrun{
#' (x <- c("Achelia assimilis", "Alcyonium antarcticum"))
#' get_bold_taxonomy(x)
#' }
get_bold_taxonomy <- function(ids, group = 1, descend = TRUE, api_rate = 0.06, ask = TRUE) {

  # set visible binding to variable
  taxon <- tax_rank <- input <- NULL

  if (!is(future::plan(), "sequential")) {

    stop("BOLD data retrieval currently do not support parallelization")

  }

  # set the api rate equal to the number of workers available
  if (is.null(api_rate)) {

    api_rate <- future::nbrOfWorkers()

  }

  ids <- unique(sort(ids)) %>% split(., ceiling(seq_along(.) / group))

  # search taxa or taxid in the NCBI
  taxa_list <- ncbi_limit_handler(ids, api_rate = api_rate, function(id) {

    taxid <- bold::bold_tax_name(ids[[id]], fuzzy = FALSE)

    if (all(is.na(taxid$taxid))) {

      return(NULL)

    } else {

      dplyr::select(taxid,
                    input,
                    taxid,
                    taxon,
                    tax_rank) %>%
        dplyr::rename(queryName = "input",
                      rank = "tax_rank")

    }

  }, message = "Searching species on taxonomy database", seed = NULL)

  # remove empty elements (no species on bold taxonomy database) and unite dataframes
  taxonomies <- taxa_list %>%
    purrr::compact() %>%
    do.call(rbind, .) %>%
    dplyr::filter(!is.na(.data$taxid))

  # get downstream taxonomic names to the species level (including intermediate
  # levels)
  if (!all(unique(taxonomies$rank) == "species") & descend) {

    nonspecies <- taxonomies[taxonomies$rank != "species",] %>% dplyr::select(c("taxon", "taxid", "rank"))

    # add a lock for when a higher taxonomic name has no children
    nonspecies$lock <- 0

    while (any(nonspecies$lock == 0)) {

      taxids <- unique(sort(nonspecies[nonspecies$lock == 0, "taxid"])) %>% split(., ceiling(seq_along(.) / 1))

      nonspecies <- ncbi_limit_handler(taxids, api_rate = api_rate, function(id) {

        downto <- nonspecies[nonspecies$taxid == taxids[[id]], "rank"] %>% get_lower_tax_rank()

        downstream <- suppressMessages(taxize::bold_downstream(taxids[[id]], downto = downto))

        while (all(downto != "species" & all(dim(downstream) %in% 0))) {

            downto <- get_lower_tax_rank(downto)

            downstream <- suppressMessages(taxize::bold_downstream(taxids[[id]], downto = downto))

        }

        if (all(dim(downstream) %in% 0)) {

          downstream <- nonspecies[nonspecies$taxid == taxids[[id]], ]
          colnames(downstream) <- c("taxon", "taxid", "rank", "lock")
          downstream$lock <- 1

        } else {

          colnames(downstream) <- c("taxon", "taxid", "rank")
          downstream$lock <- 0

          parent <- nonspecies[nonspecies$taxid == taxids[[id]], ]
          colnames(parent) <- c("taxon", "taxid", "rank", "lock")
          parent$lock <- 1

          downstream <- rbind(parent, downstream)

          if ("species" %in% unique(downstream$rank)) {

            downstream[downstream$rank == "species", ]$lock <- 1

          }

        }

        downstream

      }, message = "Searching downstream taxonomies", seed=NULL) %>% purrr::compact() %>% do.call(rbind, .)

    }

    # reorder new table with the same order as the original
    nonspecies$queryName <- nonspecies$taxon
    taxonomies_add <- nonspecies[, c("queryName", "taxid", "taxon", "rank")]

    # merge dataframes from species and non species
    taxonomies <- rbind(taxonomies, taxonomies_add) %>% unique()

  }

  # ask user to filter if there are multiple ranks
  if ((length(unique(taxonomies$rank)) > 1) && (ask == TRUE)) {

    filt_taxonomies <- ask_user(taxonomies, "rank", ask = ask)

    return(filt_taxonomies)

  } else {

    return(taxonomies)

  }

}

