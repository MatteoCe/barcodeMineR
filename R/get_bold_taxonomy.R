#' Get taxonomic classification of the BOLD taxonomy database for a list of taxa
#'
#' @param ids `character` A character string of species names.
#' @param api_rate `integer` The API rate with which to iterate each separate
#'   request. If left to default `NULL`, it will be set to one request every 16
#'   seconds approximately. Use caution when overriding the default.
#' @param descend `logical` Use the taxize package to retrieve lower level
#'   taxonomies. Default `TRUE`.
#' @inheritParams get_ncbi_taxonomy
#'
#' @return `data.frame` A data.frame object with the searched taxa, the taxid
#'   and the rank.
#'
#' @export
#' @importFrom methods is
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
#' @seealso [get_ncbi_taxonomy()]
#'
#' @examples
#' get_bold_taxonomy("Achelia assimilis", descend = FALSE)
#'
get_bold_taxonomy <- function(ids, api_rate = NULL, ask = TRUE, descend = TRUE) {

  # set grouping factor for ids names
  rate <- 1

  # set visible binding to variable
  taxon <- tax_rank <- input <- NULL

  if (!is(future::plan(), "sequential")) {

    stop("BOLD data retrieval currently do not support parallelization")

  }

  # set the api rate equal to the number of workers available
  if (is.null(api_rate)) {

    api_rate <- 0.06

  }

  ids <- unique(sort(ids)) %>% split(., ceiling(seq_along(.) / rate))

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

  }, message = "Searching taxa on taxonomy database", seed = NULL)

  if (purrr::every(taxa_list, is.null)) {
    stop(paste("No taxonomy for", ids, "was found in the BOLD database"))
  }

  # remove empty elements (no species on bold taxonomy database) and unite dataframes
  taxonomies <- taxa_list %>%
    purrr::compact() %>%
    do.call(rbind, .) %>%
    dplyr::filter(!is.na(.data$taxid))

  # get downstream taxonomic names to the subspecies level (including intermediate
  # levels)
  if (!all(unique(taxonomies$rank) == "subspecies") & descend) {

    nonspecies <- taxonomies[taxonomies$rank != "subspecies",] %>% dplyr::select(c("taxon", "taxid", "rank"))

    # add a lock for when a higher taxonomic name has no children
    nonspecies$lock <- 0

    while (any(nonspecies$lock == 0)) {

      nonspecies_done <- nonspecies[nonspecies$lock == 1, ]

      taxids <- unique(sort(nonspecies[nonspecies$lock == 0, "taxid"])) %>% split(., ceiling(seq_along(.) / 1))

      nonspecies <- ncbi_limit_handler(taxids, api_rate = api_rate, function(id) {

        downto <- nonspecies[nonspecies$taxid == taxids[[id]], "rank"] %>% get_lower_tax_rank()

        if (downto == "subspecies") {

          check_subspecies <- bold::bold_tax_id2(taxids[[id]], dataTypes = "stats") %>% t()

          if (length(check_subspecies[rownames(check_subspecies) %in% "publicsubspecies", ]) == 0) {

            downstream <- nonspecies[nonspecies$taxid == taxids[[id]], ]
            colnames(downstream) <- c("taxon", "taxid", "rank", "lock")
            downstream$lock <- 1

            return(downstream)

          }

        }

        downstream <- suppressMessages(taxize::bold_downstream(taxids[[id]], downto = downto))

        while (all(downto != "subspecies" & all(dim(downstream) %in% 0))) {

            downto <- get_lower_tax_rank(downto)

            downstream <- suppressMessages(taxize::bold_downstream(taxids[[id]], downto = downto))

        }

        if (all(downto == "subspecies" & all(dim(downstream) %in% 0))) {

          downstream <- suppressMessages(taxize::bold_downstream(taxids[[id]], downto = downto, intermediate = TRUE))

          if (is(downstream, "data.frame")) {

          } else {

            downstream <- downstream$intermediate[[1]]
            downstream$rank <- "subspecies"

            check_subspecies <- bold::bold_tax_id2(downstream$id, dataTypes = "stats") %>% dplyr::filter(., .data$publicrecords > 0)
            downstream <- dplyr::filter(downstream, id %in% check_subspecies$input)

          }

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

          if ("subspecies" %in% unique(downstream$rank)) {

            downstream[downstream$rank == "subspecies", ]$lock <- 1

          }

        }

        downstream %>% dplyr::filter(., !(.data$lock == 0 & .data$taxid %in% taxids))

      }, message = "Searching downstream taxonomies", seed=NULL) %>% purrr::compact() %>% do.call(rbind, .) %>% rbind(., nonspecies_done)

    }

    # reorder new table with the same order as the original
    nonspecies$queryName <- nonspecies$taxon
    taxonomies_add <- nonspecies[, c("queryName", "taxid", "taxon", "rank")]

    # merge dataframes from species and non species
    taxonomies <- rbind(taxonomies, taxonomies_add) %>% unique()

  }

  # reset row numbers
  rownames(taxonomies) <- seq(1, nrow(taxonomies))

  # ask user to filter if there are multiple ranks
  if ((length(unique(taxonomies$rank)) > 1) && (ask == TRUE)) {

    taxonomies <- ask_user(taxonomies, "rank", ask = ask)

  }

  # use bold_stats to count the number of records corresponding to each taxon
  bold_count <- bold_record_counter(taxonomies, api_rate)

  # full join the two tables
  dplyr::full_join(taxonomies, bold_count, by = dplyr::join_by("taxon", "taxid"))

}

