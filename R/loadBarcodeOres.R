#' Load sequences and corresponding records as a refdb object
#'
#' @param pathToRecords character string with path (or file name if in the working
#' directory) leading to strictly "tsv file" with records information.
#' @param pathToSequences character string with path (or file name if in the working
#' directory) leading to fasta file with sequences corresponding to each record.
#' @param prefix defaults to NULL. Character string that will be used to create
#' numbered custom ids for each record in ascending order.
#'
#' @return (refdb-class) a refdb object as those created with the download_*
#' functions.
#' @export
#'
#' @examples \dontrun{
#' myBO <- loadBarcodeOre("path/to/table.tsv", "path/to/sequences.fasta")
#' }
#'
loadBarcodeOre <- function(pathToRecords, pathToSequences, prefix = NULL) {

  # set error messages in case paths are incorrect and quit
  if (!(file.exists(pathToRecords))) {

    stop(paste(pathToRecords,"does not exists, is the file in another path?"))

  } else if (!(file.exists(pathToSequences))) {

    stop(paste(pathToSequences,"does not exists, is the file in another path?"))

  }

  # get sequences: simple command to load the sequences using Biostrings method
  sequences <- Biostrings::readDNAStringSet(pathToSequences, format = "fasta")

  # check that all headers have the correct format "recordID|markercode"
  badHeaders <- purrr::discard(sequences@ranges@NAMES, ~ stringr::str_detect(.x, "^[:graph:]+\\|[:graph:]+$"))

  if (length(badHeaders) > 0) {

    stop(paste("\nThese headers are not in the correct format 'recordID|markercode':\n", badHeaders))

  }

  # get records: simple command to load the records table
  records <- utils::read.delim(pathToRecords, na.strings = c("NA", ""))

  accepted_fields <- c(
    "sourceID",
    "source",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "lengthSource",
    "sampleID",
    "identified_by",
    "taxNotes",
    "db_xref",
    "NCBI_ID",
    "institutionStoring",
    "collected_by",
    "collection_date",
    "altitude",
    "depth",
    "country",
    "lat",
    "lon",
    "directionPrimers",
    "PCR_primers",
    "note"
  )

  # stop if the imported fields are not correct
  if (!(any(colnames(records) %in% accepted_fields))) {

    stop(paste(
      "\nThese fields are not accepted:\n",
      colnames(records)[!(colnames(records) %in% accepted_fields)]
    ))

  }

  # update missing unnecessary fields
  records$QueryName <- records$species

  # construct additional fields
  p <- progressr::progressor(steps = length(sequences@ranges@NAMES), on_exit = TRUE)

  records <- future.apply::future_lapply(sequences@ranges@NAMES, function(name) {

    p(message = sprintf("Final processing of records and sequences"))

    buildRecord(name, records = records, sequences = sequences)

  }) %>% purrr::compact() %>% do.call(rbind, .)

  # create refdb data.frame
  buildBarcodeOre(records, prefix)

}

