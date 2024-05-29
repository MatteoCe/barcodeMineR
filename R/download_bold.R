#' Download records from BOLD
#'
#' @param bold_tax a data.frame, as returned from the get_bold_taxonomy function.
#' @param rate numeric, the number of taxonomic names for grouping. These will be
#' queried at the same time.
#' @param api_rate numeric, the timing rate at which to send requests. In this
#' function, it defaults to NULL and is always set to 1, one request per second,
#' due to the incompatibility of the bold functions with the future parallelizing
#' framework.
#' @param ask logical, defaults to TRUE. If disabled the selection process will
#' be skipped and all feature keys from each accession number will be retrieved.
#' @param prefix character, defaults to NULL. defaults to NULL. Character string
#' that will be used to create numbered custom ids for each record in ascending
#' order.
#'
#' @return a refdb tibble data.frame, including the DNA sequence as a field.
#'
#' @export
#'
#' @description
#' A short description...
#'
download_bold <- function(bold_tax, rate = 100, api_rate = NULL, ask = TRUE, prefix = NULL) {

  if (!is(future::plan(), "sequential")) {
    stop("BOLD data retrieval currently do not support parallelization")
  }

  # set the api rate equal to the number of workers available
  if (is.null(api_rate)) {
    api_rate <- future::nbrOfWorkers()
  }

  # use bold_stats to count the number of records corresponding to each taxon
  bold_count <- bold_record_counter(bold_tax, api_rate)

  # divide the taxa based on the number of records. The rate parameter will group
  # taxa if the cumulative sum of the corresponding records do not exceed "rate".
  # Taxa corresponding to more than "rate" records will be searched alone.
  ids_groups <- bold_record_grouper(bold_count[order(bold_count$records), ], bold_tax, rate)

  # download records
  records_tab <- ncbi_limit_handler(ids_groups, api_rate = api_rate, function(id) {

    bold_fetcher(ids_groups[[id]], bold_count, bold_tax)

  }, message = "Downloading BOLD records", seed = NULL) %>% purrr::compact() %>% do.call(rbind, .)

  # format and remove unwanted records
  records_tab <- extractRecordsTabBOLD(records_tab, bold_tax)

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
