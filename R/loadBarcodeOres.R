#' Load sequences and corresponding records as a refdb object
#'
#' @param records `data.frame` or `character` A character string with path (or
#'   file name if in the working directory) leading to strictly "tsv file" with
#'   records information. It can also correspond to a data frame, provided that
#'   it has the fields included in the example data 'barcodeMineR::example_record'.
#' @param sequences `DNAStringSet` or `character` A character string with path
#'   (or file name if in the working directory) leading to fasta file with
#'   sequences corresponding to each record. It can also correspond to a
#'   "DNAStringSet" object, as the one in the example data
#'   `barcodeMineR::example_sequences`. The name of each sequence from the
#'   DNAStringSet object must correspond to the concatenation of the fields
#'   `sourceID` and `markerCode`, separated by a pipe `|`.
#' @param prefix `character` A character string that will be used to create
#'   numbered custom ids for each record in ascending order. The prefix will
#'   compose the recordID field in the final object. Default to `NULL`, using
#'   the information extracted from the field `sourceID`.
#'
#' @return `data.frame` A refdb data frame, including the DNA sequence as a
#'   field.
#'
#' @export
#'
#' @description
#' This function allows the user to load a custom refdb-formatted data frame
#' object with additional records and sequences obtained from private analyses.
#' The objects can be loaded from both files (tsv and fasta) or loaded objects
#' (data.frame and DNAStringSet).
#'
#' @examples \dontrun{
#' # load from tsv and fasta files
#' loadBarcodeOre("path/to/table.tsv", "path/to/sequences.fasta")
#' }
#'
#' # load from data.frame and DNAStringSet loaded objects
#' loadBarcodeOre(example_record, example_sequence)
#'
loadBarcodeOre <- function(records, sequences, prefix = NULL) {

  # check which type of data has been supplied to the function
  if (all(c(class(records), class(sequences)) %in% "character")) {

    # set error messages in case paths are incorrect and quit
    if (!(file.exists(records))) {

      stop(paste(records,"does not exists, is the file in another path?"))

    } else if (!(file.exists(sequences))) {

      stop(paste(sequences,"does not exists, is the file in another path?"))

    } else {

      # get sequences: simple command to load the sequences using Biostrings method
      sequences <- Biostrings::readDNAStringSet(sequences, format = "fasta")

      # get records: simple command to load the records table
      records <- utils::read.delim(records, na.strings = c("NA", ""))

    }

  } else if (all("data.frame" %in% class(records) & class(sequences) == "DNAStringSet")) {

  } else {

    stop("The function does not support object types other than 'character' or 'data.frame' and 'DNAStringSet' as inputs")

  }

  # check that all headers have the correct format "recordID|markercode"
  badHeaders <- purrr::discard(sequences@ranges@NAMES, ~ stringr::str_detect(.x, "^[:graph:]+\\|[:graph:]+$"))

  if (length(badHeaders) > 0) {

    stop(paste("\nThese headers are not in the correct format 'recordID|markercode':\n", badHeaders))

  }

  # set correct fields
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

