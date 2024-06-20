#' Merge multiple refdb objects to create a single one
#'
#' @param ... `data.frame` Any number of refdb-formatted data frames, as those
#'   obtained with `download_ncbi`, `download_bold` and `loadBarcodeOre`.
#' @param resolve.conflicts `logical` If set to FALSE, the script will merge
#'   the refdb objects from different sources (BOLD, NCBI and custom) without
#'   any quality control steps. Otherwise, it searches for "mining duplicates",
#'   records that appear in both sources due to internal mining, and returns
#'   only the original record. Defaults to `TRUE`.
#' @return `data.frame` A refdb object, including the records and sequences from
#'   all refdb objects provided as arguments.
#'
#' @export
#' @importFrom rlang .data
#'
#' @description
#' This function merges the output of multiple data frames obtained using the
#' `download_bold` and `download_ncbi` functions or from private data obtained
#' with the function `loadBarcodeOre`. It resolves conflicts originated from
#' internal mining performed by the online database, which would result in
#' duplicates in the final object. This operation can be avoided using the
#' argument `resolve.conflicts`.
#'
#' @examples
#' # search and download Maldane sarsi records:
#' tax_ncbi <- get_ncbi_taxonomy("Maldane sarsi", ask = FALSE)
#' tax_bold <- get_bold_taxonomy("Maldane sarsi", ask = FALSE)
#' rec_ncbi <- download_ncbi(tax_ncbi, ask = FALSE)
#' rec_bold <- download_bold(tax_bold, ask = FALSE)
#'
#' # merge all results into one
#' mergeBarcodeOres(rec_ncbi, rec_bold)
#'
mergeBarcodeOres <- function(..., resolve.conflicts = TRUE) {

  args <- list(...)

  refdbs <- do.call(rbind, args)

  # get vector of unique identifiers in merged record tabs
  refdb_ids <- paste(refdbs$recordID, refdbs$markerCode, sep = "|")

  # check if there are duplicate unique identifiers
  if (nrow(refdbs) != length(unique(refdb_ids))) {

    message("Duplicate records were found:\n")
    duplicates <- refdbs %>% tidyr::unite("ids", .data$recordID, .data$markerCode, sep = "|") %>%
      dplyr::count(.data$ids) %>%
      dplyr::filter(.data$n > 1)
    print(duplicates)

  }

  # conditional crossroad, is the user intends to remove "duplicate" records
  # usually obtained from a concurrent search of species on both NCBI and BOLD
  # and wants to merge the products. Basically, all records from BOLD that were
  # mined from the NCBI by the database mining process will be removed if also
  # here represented by NCBI records. Similarly, all those NCBI records that
  # originated from BOLD will be searched in the BOLD records and eventually
  # removed from the former if found. This allows to have only one representative
  # (from either NCBI or BOLD) of any occurrence.

  if (resolve.conflicts) {

    if (!(all(c("NCBI", "BOLD") %in% unique(refdbs$source)))) {

      message("The resolve.conflicts parameters should be set only for refdb objects obtained following the suggested pipeline/workflow.")

    } else {

      procs_BOLD_remove <- procs_NCBI_remove <- NULL

      if (nrow(refdbs %>% dplyr::filter(.data$note == "Mined from GenBank, NCBI")) > 0) {

        procs_BOLD_remove <- boldDuplicateCheck(refdbs)

      }

      if (length(grep("BOLD:", refdbs$db_xref)) > 0) {

        procs_NCBI_remove <- ncbiDuplicateCheck(refdbs)

      }

      procs_remove <- c(procs_BOLD_remove, procs_NCBI_remove)

      if (is.null(procs_remove)) {

        message("No records were obtained from both the NCBI and BOLD")

        final_refdb <- buildBarcodeOre(refdbs)
        return(final_refdb)

      } else {

        final_refdb <- refdbs %>% dplyr::filter(!(.data$recordID %in% procs_remove)) %>%
          buildBarcodeOre()

        return(final_refdb)

      }}

  }

  # else, if resolve.conflicts is set to FALSE, simply merge all data from all Ores
  final_refdb <- buildBarcodeOre(refdbs)
  return(final_refdb)

}
