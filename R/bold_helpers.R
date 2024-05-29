#' Select the taxa with records on bold
#'
#' @param bold_tax data.frame the output of get_bold_taxonomy
#' @param api_rate the api_rate as set in download_bold
#'
#' @return data.frame with taxa represented by at least 1 record on bold
#'
#' @keywords internal
#' @importFrom rlang .data
#'
#'
bold_record_counter <- function(bold_tax, api_rate) {

  purrr::map(unique(bold_tax$rank), function(rank_sel) {

    bold_tax_sel <- dplyr::filter(bold_tax, rank == rank_sel)

    taxa <- unique(sort(bold_tax_sel$taxon)) %>% split(., ceiling(seq_along(.) / 10))

    # use bold_stats to count the number of records corresponding to each taxon
    bold_count <- ncbi_limit_handler(taxa, api_rate = api_rate, function(id) {

      idStats <- bold::bold_stats(taxon = taxa[[id]])

      if (all(idStats[[rank_sel]]$drill_down$entity$records == 0)) {

        return(NULL)

      }

      idStats[[rank_sel]]$drill_down$entity %>% dplyr::filter(.data$records != 0)

    }, message = "Counting number of records per taxa on BOLD", seed = NULL) %>% purrr::compact() %>% do.call(rbind, .)

  }) %>% do.call(rbind, .)

}

#' Group taxonomic names by maximum cumulative value
#'
#' @param bold_count data.frame output of the bold_record_counter function
#' @param bold_tax data.frame output of the get_bold_taxonomy function
#' @param rate maximum number of records to group
#'
#' @return a list of character vectors
#'
#' @keywords internal
#' @importFrom rlang .data
#'
bold_record_grouper <- function(bold_count, bold_tax, rate) {

  # set visible binding to variable
  name <- NULL

  groups <- list()
  counter <- 0

  for (rank in unique(bold_tax$rank)) {

    bold_count_sel <- bold_count[bold_count$name %in% bold_tax[bold_tax$rank == rank, "taxon"], ]

    while (nrow(bold_count_sel) > 0) {

      counter <- counter + 1

      group <- bold_count_sel %>% dplyr::filter(cumsum(.data$records) < rate)

      # if no record is filtered, but bold_counts still has rows, the first
      # record will be filtered as it will correspond to records > rate
      if (nrow(group) == 0) {

        group <- bold_count_sel[1, ]

      }

      group <- group %>% dplyr::select(name) %>% c(., recursive = TRUE) %>% unname()

      groups[[counter]] <- group

      bold_count_sel <- bold_count_sel %>% dplyr::filter(!(.data$name %in% group))

    }

  }

  groups

}

#' Download records from the BOLD database
#'
#' @param taxon_group taxonomic names
#' @param bold_counts data.frame output of bold_record_grouper
#' @param bold_tax data.frame output of get_bold_taxonomy
#'
#' @return a data.frame
#'
#' @keywords internal
#' @importFrom rlang .data
#'
bold_fetcher <- function(taxon_group, bold_count, bold_tax) {

  records_number <- bold_count[bold_count$name %in% taxon_group, "records"] %>% sum()

  idRecords <- bold::bold_seqspec(taxon = taxon_group,
                                  seqFasta = FALSE,
                                  format = "tsv")

  if (records_number != length(unique(sort(idRecords$processid)))) {
    stop("The number of records from bold_stats do not correspond to those retrieved by bold_seqspec")
  }

  # some columns might have formats that do not allow to properly clean the
  # entire table (like POSIX) and thus must be changed to character. In any case,
  # most often these columns are not necessary for the following operations
  # (run_dates, for example)

  # get the class of each column
  classes <- lapply(idRecords, class) %>% do.call(c, .)

  # if some are not in the specified formats, change to character
  if (any(!(classes %in% c("character", "integer", "logical", "numeric")))) {

    idRecords <- dplyr::mutate(idRecords, dplyr::across(which(!(classes %in% c("character",
                                                                               "integer",
                                                                               "logical",
                                                                               "numeric"))), as.character))

  }

  # clean table from empty strings
  idRecords[idRecords == ""] <- NA

  # some records may not have sequences, "ASULB1341-21" may be an example.
  # Thus, if the nucleotide field is empty for some rows, remove it.
  if (any(is.na(idRecords$nucleotides))) {

    # check if all records do not have sequences
    if (nrow(idRecords[is.na(idRecords$nucleotides), ]) ==
        records_number) {

      return(NULL)

    } else {

      # clean the records table from any record with "NA" as DNA sequence
      idRecords <- idRecords[!is.na(idRecords$nucleotides), ]

    }

  }

  # filter records to include only those corresponding to the searched taxonomic
  # name, not downstream taxonomy records
  purrr::map(unique(bold_tax[bold_tax$taxon %in% taxon_group, "rank"]), function(rank) {

    lower_ranks <- NULL

    rank_sel <- rank

    while (rank_sel != "species") {

      rank_sel <- get_lower_tax_rank(rank_sel)
      lower_ranks <- c(lower_ranks, paste0(rank_sel, "_name"))

    }

    rank_sel <- paste0(rank, "_name")

    idRecords_sel <- idRecords %>% dplyr::filter(.data[[rank_sel]] %in% bold_tax[bold_tax$rank == rank, ]$taxon) %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(lower_ranks), is.na))

    if (nrow(idRecords_sel) == 0) {
      return(NULL)
    } else {
      idRecords_sel
    }

  }) %>% purrr::compact() %>% do.call(rbind, .)

}

#'  Extract records table from the output of bold::bold_seqspec
#'
#' @param bold_records data.frame the output of bold::bold_seqspec
#' @param bold_tax data.frame the output of get_bold_taxonomy
#'
#' @return data.frame a records_tab with selected columns
#'
#' @keywords internal
#'
#' @description
#' This function is very similar to extractRecordsTab, thus allowing to extract
#' records information, but from the output of a bold_fetcher function, which gives
#' a "tsv-formatted" records table, very easily processable to create the records
#' tab format used here. The taxonomic table, obtained with get_bold_taxonomy is
#' here inspected to extract solely the BOLD taxid for the specific species processed
#' at each iteration.
#'
extractRecordsTabBOLD <- function(bold_records, bold_tax) {

  # set progressbar
  p <- progressr::progressor(steps = nrow(bold_records))
  message <- "Formatting BOLD records..."

  # if some records have been retrieved searching a higher taxon, this will search
  # iteratively on all upstream taxonomies until finding the taxonomy in bold_tax.
  # In this way, the taxonomy table will have the correct queryName
  purrr::map(1:nrow(bold_records), function(row) {

    taxonomic_table <- data.frame()
    counter <- 0

    while (nrow(taxonomic_table) == 0) {

      counter <- counter + 1

      name <- bold_records[row,] %>% dplyr::select(tidyselect::matches("_name$")) %>%
        dplyr::select(dplyr::where(~!all(is.na(.x)))) %>%
        c(., recursive = TRUE) %>%
        utils::tail(., n = counter) %>% .[1] %>% unname()

      taxonomic_table <- bold_tax[bold_tax$taxon == name, ]

    }

    taxid <- bold_records[row,] %>% dplyr::select(tidyselect::matches("_name$")) %>%
      dplyr::select(dplyr::where(~!all(is.na(.x)))) %>%
      c(., recursive = TRUE) %>%
      utils::tail(., n = 1) %>%
      names(.) %>%
      stringr::str_replace(., "_name", "_taxID") %>%
      bold_records[row, .]

    # calculate entire and without gaps length of the sequence
    lengthSource <- stringr::str_count(bold_records[row, "nucleotides"], "[:graph:]")
    lengthGene <- stringr::str_count(bold_records[row, "nucleotides"], "[:alpha:]")

    # extract info from the institution_storing field, unless it is Mined from NCBI
    institution <- ifelse((bold_records[row, "institution_storing"] == "Mined from GenBank, NCBI" | is.na(bold_records[row, "institution_storing"])),
                          NA, bold_records[row, "institution_storing"])

    # concurrently, if the record is mined from NCBI, store that information in another field
    note <- ifelse(bold_records[row, "institution_storing"] == "Mined from GenBank, NCBI",
                   "Mined from GenBank, NCBI",
                   NA)

    # merge start date and end date
    date <- ifelse((is.na(bold_records[row, "collectiondate_start"]) | is.na(bold_records[row, "collectiondate_end"])),
                   NA, paste(bold_records[row, "collectiondate_start"], bold_records[row, "collectiondate_end"], sep = "|"))

    # send input to progressr
    if (!is.null(p)) {

      p(message = sprintf(message))

    }

    # create data.frame with the previously extracted information
    data.frame(
      source = "BOLD",
      sourceID = bold_records[row, "processid"],
      DNA_seq = bold_records[row, "nucleotides"],
      markerCode = bold_records[row, "markercode"],
      lengthGene = lengthGene,
      sampleID = bold_records[row, "sampleid"],
      QueryName = taxonomic_table$queryName,
      phylum = bold_records[row, "phylum_name"],
      class = bold_records[row, "class_name"],
      order = bold_records[row, "order_name"],
      family = bold_records[row, "family_name"],
      genus = bold_records[row, "genus_name"],
      species = bold_records[row, "species_name"],
      identified_by = bold_records[row, "identification_provided_by"],
      taxNotes = bold_records[row, "tax_note"],
      db_xref = paste0("taxon:", taxid),
      NCBI_ID = bold_records[row, "genbank_accession"],
      institutionStoring = institution,
      collected_by = bold_records[row, "collectors"],
      collection_date = date,
      lat = bold_records[row, "lat"],
      lon = bold_records[row, "lon"],
      altitude = bold_records[row, "elev"],
      depth = bold_records[row, "depth"],
      country = bold_records[row, "country"],
      directionPrimers = bold_records[row, "directions"],
      lengthSource = lengthSource,
      PCR_primers = bold_records[row, "seq_primers"],
      note = note
    )

  }) %>% purrr::compact() %>% do.call(rbind, .)

}
