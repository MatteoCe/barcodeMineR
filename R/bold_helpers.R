#' Select the taxa with records on bold
#'
#' @param bold_tax data.frame the output of get_bold_taxonomy
#' @param api_rate the api_rate as set in download_bold
#'
#' @return data.frame with taxa represented by at least 1 record on bold
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom rlang .data
#'
bold_record_counter <- function(bold_tax, api_rate) {

  # set visible binding to variable
  name <- records <- taxon <- taxid <- NULL

  taxa <- unique(sort(bold_tax$taxon)) %>% split(., ceiling(seq_along(.) / 1))

  # use bold_stats to count the number of records corresponding to each taxon
  stats <- ncbi_limit_handler(taxa, api_rate = api_rate, function(id) {

    idStats <- bold::bold_stats(taxa[[id]], dataTypes = "drill_down", simplify = TRUE)

    rank <- bold_tax[bold_tax$taxon == taxa[[id]], "rank"]

    if (length(rank) > 1) {

      taxids <- bold_tax[bold_tax$taxon == taxa[[id]], "taxid"]

      taxidStats <- bold::bold_tax_id2(taxids, dataTypes = "basic")

      stats <- lapply(1:nrow(taxidStats), function(row) {

        parent <- taxidStats[row, "parentname"]
        rank <- taxidStats[row, "tax_rank"]

        upper_rank <- get_lower_tax_rank(rank, upper = TRUE)

        if (upper_rank == "subfamily") {
          upper_rank <- "family"
        }

        rec_stats <- idStats$drill_down[[upper_rank]]

        if (!(parent %in% rec_stats$name)) {
          records <- rec_stats[rec_stats[, 1] == "Unspecified*", "records"]
        } else {
          records <- rec_stats[rec_stats[, 1] == parent, "records"]
        }

        data.frame(
          taxid = taxids[row],
          taxon = taxa[[id]],
          parentname = parent,
          records = records)

      }) %>% do.call(rbind, .)

    } else {

      upper_rank <- get_lower_tax_rank(rank, upper = TRUE)

      if (rank == "order") {
        stop("The function 'download_bold' does not support downloading taxa higher than the family level.\nProvide lower taxonomic levels with 'get_bold_taxonomy'")
      }

      if (rank != "subspecies" & is.null(idStats$drill_down[[rank]]$records)) {
        return(NULL)
      }

      if (rank == "subspecies" & is.null(idStats$drill_down[[rank]]$records)) {
        rank <- upper_rank <- "species"
      }

      if (upper_rank == "subfamily") {
        upper_rank <- "family"
      }

      data.frame(
        taxid = bold_tax[bold_tax$taxon == taxa[[id]], "taxid"],
        taxon = taxa[[id]],
        parentname = idStats$drill_down[[upper_rank]]$name[!(idStats$drill_down[[upper_rank]]$name %in% "Others")],
        records = sum(idStats$drill_down[[rank]]$records))

    }

  }, message = "Counting number of records per taxa on BOLD", seed = NULL) %>% purrr::compact() %>% do.call(rbind, .)

  final_stats <- stats %>% dplyr::filter(., !(.data$parentname %in% .data$taxon)) %>%
    dplyr::select(., taxon, taxid, records)

  if (all(final_stats$records == 0)) {

    return(NULL)

  }

  final_stats %>% dplyr::filter(.data$records != 0)

}

#' Group taxonomic names by maximum cumulative value
#'
#' @param bold_tax data.frame output of the get_bold_taxonomy function
#' @param rate maximum number of records to group
#'
#' @return a list of character vectors
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom rlang .data
#'
bold_record_grouper <- function(bold_tax, rate) {

  # set visible binding to variable
  name <- taxon <- records <- NULL

  # extract counts from bold_tax
  bold_count <- bold_tax[!is.na(bold_tax$records), c("taxon", "records")] %>%
    dplyr::group_by(taxon) %>%
    dplyr::summarise(records = sum(records))

  groups <- list()
  counter <- 0

  while (nrow(bold_count) > 0) {

    counter <- counter + 1

    group <- bold_count %>% dplyr::filter(cumsum(.data$records) < rate)

    # if no record is filtered, but bold_counts still has rows, the first
    # record will be filtered as it will correspond to records > rate
    if (nrow(group) == 0) {

      group <- bold_count[1, ]

    }

    group <- group %>% dplyr::select(taxon) %>% c(., recursive = TRUE) %>% unname()

    groups[[counter]] <- group

    bold_count <- bold_count %>% dplyr::filter(!(.data$taxon %in% group))

  }

  groups

}

#' Download records from the BOLD database
#'
#' @param taxon_group taxonomic names
#' @param bold_tax data.frame output of get_bold_taxonomy
#'
#' @return a data.frame
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom rlang .data
#'
bold_fetcher <- function(taxon_group, bold_tax) {

  records_number <- bold_tax[bold_tax$taxon %in% taxon_group, "records"] %>% sum()

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

  idRecords

}

#'  Extract records table from the output of bold::bold_seqspec
#'
#' @param bold_records data.frame the output of bold::bold_seqspec
#' @param bold_tax data.frame the output of get_bold_taxonomy
#' @param filter whether to filter the downloaded output for the exact
#'   taxonomies present in the bold_tax or to include children taxonomy records
#'
#' @return data.frame a records_tab with selected columns
#'
#' @keywords internal
#' @noRd
#'
#' @description
#' This function is very similar to extractRecordsTab, thus allowing to extract
#' records information, but from the output of a bold_fetcher function, which gives
#' a "tsv-formatted" records table, very easily processable to create the records
#' tab format used here. The taxonomic table, obtained with get_bold_taxonomy is
#' here inspected to extract solely the BOLD taxid for the specific species processed
#' at each iteration.
#'
extractRecordsTabBOLD <- function(bold_records, bold_tax, filter = TRUE) {

  # set progressbar
  p <- progressr::progressor(steps = nrow(bold_records))
  message <- "Formatting BOLD records..."

  # filter records to include only those corresponding to the searched taxonomic
  # name, not downstream taxonomy records
  if (filter) {
    bold_records <- purrr::map(seq(1, nrow(bold_tax)), function(row) {

      lower_ranks <- NULL

      rank_sel <- rank <- bold_tax[row, "rank"]

      while (rank_sel != "subspecies") {

        rank_sel <- get_lower_tax_rank(rank_sel)
        lower_ranks <- c(lower_ranks, paste0(rank_sel, "_name"))

      }

      rank_sel <- paste0(rank, "_name")

      idRecords_sel <- bold_records %>% dplyr::filter(.data[[rank_sel]] %in% bold_tax[row, "taxon"]) %>%
        dplyr::filter(dplyr::if_all(dplyr::all_of(lower_ranks), is.na))

      if (nrow(idRecords_sel) == 0) {
        return(NULL)
      } else {
        idRecords_sel
      }

    }) %>% purrr::compact() %>% do.call(rbind, .)
  }
  # if no records remain, stop and message it
  if (length(bold_records) == 0) {
    return(NULL)
  }

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
