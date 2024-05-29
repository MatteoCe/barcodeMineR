#' Merge multiple refdb objects to create a single one.
#'
#' @param ... Any number of refdb objects.
#' @param resolve.conflicts (default = TRUE) If set to FALSE, the script will
#' merge refdb objects from different sources without any quality control steps.
#' If left to default settings, it is used to perform a basic quality control step
#' required if the refdb objects were obtained from a pipeline analysis that
#' searched the same species in both the NCBI and BOLD databases. Records from
#' both databases.
#'
#' @return A refdb object, including the records and sequences from all refdb
#' objects provided as arguments.
#'
#' @export
#' @importFrom rlang .data
#'
#' @examples \dontrun{
#' (totBO <- mergeBarcodeOres(ncbiBO, boldBO, customBO, resolve.conflicts = TRUE))
#'
#' }
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
