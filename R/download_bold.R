#' Download records from BOLD
#'
#' @param bold_tax a data.frame, as returned from the get_bold_taxonomy function.
#' @param rate `integer` The number of records to be downloaded at a time. It
#'   can be lowered for unstable internet connections. However, due to the
#'   structure of the bold package, it is not possible to download a specific
#'   number of records when a species is represented by more than `rate` records.
#'   This function only groups species whose sum of records is inferior to rate.
#'   Defaults to `100`.
#' @param api_rate `integer` The API rate with which to iterate each separate
#'   request. Must be a number between 3 and 10 which will translate in a rate
#'   of `1 / api_rate` seconds.
#' @param ask `logical` Should the function ask the user whether to filter the
#'   final output for taxonomic ranks. Default `TRUE`.
#' @param prefix `character` A character string that will be used to create
#'   numbered custom ids for each record in ascending order. The prefix will
#'   compose the recordID field in the final object. Default to `NULL`, using
#'   the internal recordID generator that will use the accession number for NCBI
#'   records and the processID for BOLD records, avoiding duplicates by adding
#'   `_1`, `_2` etc.
#'
#' @return `data.frame` A refdb data frame, including the DNA sequence as a
#'   field.
#'
#' @export
#'
#' @description
#' This function searches for BOLD records corresponding to the species found in
#' the argument `bold_tax`, i.e. the output from the function `get_bold_taxonomy`.
#'
#' For a thorough explanation of the function usage and capabilities, see the
#' 'Introduction to the barcodeMineR package' vignette:
#' \code{vignette("Introduction to the barcodeMineR package", package = "barcodeMineR")}
#'
#' @seealso [download_ncbi()]
#'
#' @examples
#' tax <- get_bold_taxonomy("Polymastia invaginata")
#'
#' download_bold(tax, ask = FALSE)
#'
download_bold <- function(bold_tax, rate = 100, api_rate = NULL, ask = TRUE, prefix = NULL) {

  if (!is(future::plan(), "sequential")) {
    stop("BOLD data retrieval currently do not support parallelization")
  }

  # set the api rate equal to the number of workers available
  if (is.null(api_rate)) {
    api_rate <- future::nbrOfWorkers()
  }

  # divide the taxa based on the number of records. The rate parameter will group
  # taxa if the cumulative sum of the corresponding records do not exceed "rate".
  # Taxa corresponding to more than "rate" records will be searched alone.
  ids_groups <- bold_record_grouper(bold_tax, rate)

  # download records
  records_tab <- ncbi_limit_handler(ids_groups, api_rate = api_rate, function(id) {

    bold_fetcher(ids_groups[[id]], bold_tax)

  }, message = "Downloading BOLD records", seed = NULL) %>% purrr::compact() %>% do.call(rbind, .)

  # format and remove unwanted records
  records_tab <- extractRecordsTabBOLD(records_tab, bold_tax)

  # stop if no records are found
  if (is.null(records_tab)) {
    stop(paste("No", paste(bold_tax$taxon, collapse = ", "),"records were found on the BOLD database.\nRead the vignette 'Searching taxonomy' at https://matteoce.github.io/barcodeMineR/index.html for more information."))
  }

  ### Select (step)
  selection_tab <- data.frame(BOLD_processID = records_tab$sourceID,
                              Gene_value = records_tab$markerCode)

  names <- select_accessions(selection_tab, ask = ask) %>% unique() %>% .[!is.na(.)]

  # change format to fasta and then to DNAStringSet (it will only be used in the
  # next function to check records mined from genbank multiple times with different
  # lenghts)
  dnaStringFormat <- textToDNAStringSet(records_tab)

  # first check if there are records "mined from genbank" and if so, run a cleaning
  # function that allows to remove duplicates from the record tables, records that
  # were mined multiple times from the NCBI in the BOLD database but correspond to
  # the same sequence/record. The cleaning function will check which duplicate
  # record has the correct length after fetching the original NCBI record using the
  # download_ncbi_records.R script, if this is not sufficient to discriminate the
  # records, the user will be asked to choose which BOLD duplicate record to keep.

  names_sel <- cleanMinedFromGenBankRecords(records_tab, dnaStringFormat, names)

  # apply buildRecords to all chosen accession numbers|gene_name, using the new
  # total sequences as they contain the correct region of the chosen CDS or rRNA
  p <- progressr::progressor(steps = length(names_sel), on_exit = TRUE)

  records_total <- future.apply::future_lapply(names_sel, function(name) {

    p(message = sprintf("Final processing of records and sequences"))

    buildRecord(name, records = records_tab, sequences = dnaStringFormat)

  }) %>% purrr::compact() %>% do.call(rbind, .)

  # create refdb data.frame
  buildBarcodeOre(records_total, prefix)

}
