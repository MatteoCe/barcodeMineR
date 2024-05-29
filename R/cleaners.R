#' Clean duplicate BOLD records that were mined from the NCBI multiple times
#'
#' @param records data.frame the output of bold fetcher format
#' @param sequences DNAStringSet the sequences corresponding to the records table, but in DNAStringSet format
#' @param names character a vector of processID and corresponding marker codes as selected prevously by the user
#'
#' @return character the same as parameter name, but including the correct mined record
#'
#' @keywords internal
#' @importFrom rlang .data
#'
#' @description
#' This is a convenient function that allows cleaning duplicate BOLD records that
#' were mined from the NCBI but multiple times and thus have a different processID.
#' In those cases, the accession numbers those records reports will be the same,
#' thus the original NCBI record is searched using different utils functions listed
#' in utils_helpers.R to search those records and check the length of the original
#' CDS/rRNA.
#'
#' Create, from the records object, a new table that includes two new columns
#' one corresponding to the combination of processID and gene name that is also
#' supplied to this function as accn, and another with the corresponding accession
#' number. The latter is created using the coalesce function of dplyr which allows
#' to merge two vectors, in this case columns of a table, replacing the NA items
#' with the corresponding value, if present, in the other column. Important to
#' remember that coalesce works giving preference to the first vector supplied,
#' thus if both columns present the value the one coming from the column
#' NCBI_ID will be used instead of the other. This allows to recover all
#' those accession numbers that are not present in the correct field of the record
#' mined from Genbank but instead in the sampleID field.
#' The lengthGene field is quickly calculated using the sequences supplied to the
#' function. Moreover, some records on BOLD have no sequence associated. Those
#' will be removed beforehand.
#'
cleanMinedFromGenBankRecords <- function(records, sequences, accn) {

  mined <- records %>% dplyr::filter(.data$note == "Mined from GenBank, NCBI")

  if (nrow(mined) == 0) {

    # no records were mined from genbank, exit and continue with all accn
    return(accn)

  }

  # get chosen markerCode names and processIDs stored in the accn argument to help selection
  markerCodes <- stringr::str_remove(string = accn, pattern = ".*\\|") %>% unique()

  filt_mined <- records %>% dplyr::filter(.data$note == "Mined from GenBank, NCBI",
                                          .data$sourceID %in% stringr::str_remove(accn, "\\|.*"),
                                          .data$markerCode %in% markerCodes) %>%
    dplyr::mutate(id = paste(.data$sourceID, .data$markerCode, sep = "|"),
                  accn = dplyr::coalesce(.data$NCBI_ID, .data$sampleID)) %>%
    dplyr::mutate(lengthGene = sequences[.data$id]@ranges@width)

  # get vector of accession numbers without the final ".1"
  ac_norm <- stringr::str_remove(filt_mined$accn, "\\.1")

  # select only those occurring multiple times
  duplicates <- unique(ac_norm[duplicated(ac_norm)])

  if (length(duplicates) == 0) {

    # No duplicates, continue normally with "parent" script
    #message("No duplicates in all records that were mined from the NCBI, continue...")
    return(accn)

  } else {

    # check if the duplicates are not simply multiple markers from the same big NCBI sequence
    areMultipleMarkers <- filt_mined[filt_mined$accn %in% duplicates, ] %>% dplyr::filter(duplicated(.data$markerCode) & duplicated(.data$accn))

    if (nrow(areMultipleMarkers) == 0) {
      return(accn)
    }

    # search the accession numbers on the NCBI using the download_ncbi_records function
    # This greatly increases speed and simplify this code as well.
    ncbi_records <- download_ncbi(ncbi_ids = paste(duplicates, "1", sep = "."))

    wrong_accn_list <- lapply(duplicates, function(ac) {

      # extract original length information corresponding to that accession number
      lengthOriginal <- ncbi_records %>% dplyr::filter(.data$NCBI_ID == paste0(ac, ".1")) %>%
        dplyr::select(.data$lengthGene) %>%
        c(.data, recursive = TRUE) %>%
        unname

      # select, from the duplicates found for that accession number in the BOLD
      # retrieved records, that record including the length of the gene that exactly
      # corresponds to the original
      filt_rec <- filt_mined %>% dplyr::filter(.data$accn == ac | .data$accn == paste0(ac, ".1")) %>%
        dplyr::filter(.data$lengthGene %in% .data$lengthOriginal)

      # check if the selected information is unique and thus unequivocal as correct
      # record to include. Return the other, incorrect, which at the end will be used
      # to update the "accn" parameter supplied to this function (cleanMinedFromGenBankRecords)
      # removing the wrong record. Otherwise, if both (or none) records show the
      # correct length, then ask the user to choose which record, of the two mined
      # from genbank, should be kept, showing the user the available information
      # for the two duplicated records.

      if (nrow(filt_rec) == 1) {

        wrong_accn <- filt_mined %>% dplyr::filter((.data$accn == ac | .data$accn == paste0(ac, ".1")),
                                                   .data$sourceID != filt_rec$sourceID) %>%
          dplyr::select(.data$id) %>%
          c(., recursive = TRUE) %>%
          unname()

        return(wrong_accn)

      } else {

        # get original NCBI complete record and transpose
        filt_original <- ncbi_records %>% dplyr::filter(.data$NCBI_ID == paste0(ac, ".1"))

        # transpose the duplicates table for visual clearness and bind the original
        # NCBI info
        filt_rec_t <- filt_mined %>% dplyr::filter(.data$accn == ac | .data$accn == paste0(ac, ".1")) %>%
          dplyr::select(-.data$accn, -.data$id) %>%
          dplyr::bind_rows(.data, filt_original) %>%
          dplyr::select(-.data$DNA_seq, -.data$note) %>%
          t()

        # change column names to 1, 2 and originalNCBI
        colnames(filt_rec_t) <- c(seq(1:(ncol(filt_rec_t) - nrow(filt_original))), rep("originalNCBI", nrow(filt_original)))

        # Ask the user to choose the correct record. Both records to choose from
        # and the original NCBI record here downloaded will be printed.
        repeat {

          print(filt_rec_t)

          message("Length confirmation method failed. Choose which 'Mined from the NCBI' BOLD record to keep.\nIndicate the column number of the processID to keep:")
          numb <- scan(n = 1, quiet = TRUE, what = "raw")

          if ((length(numb) != 0) & (numb %in% colnames(filt_rec_t)))
            break

        }

        procs_select <- filt_rec_t[2, which(colnames(filt_rec_t) != c(numb, "originalNCBI"))] %>%
          unname() %>%
          unique()

        wrong_accn <- filt_mined %>% dplyr::filter((.data$accn == ac | .data$accn == paste0(ac, ".1")),
                                                   .data$sourceID %in% procs_select) %>%
          dplyr::select(.data$id) %>%
          c(., recursive = TRUE) %>%
          unname()
        return(wrong_accn)

      }

    })

    wrongs <- do.call(c, wrong_accn_list)

    clean_accns <- accn[!(accn %in% wrongs)]
    return(clean_accns)

  }

}

#' Remove BOLD records if the corresponding, mined NCBI records are present
#'
#' @param refdb a refdb object
#'
#' @return a refdb object
#'
#' @keywords internal
#' @importFrom rlang .data
#'
boldDuplicateCheck <- function(refdb) {

  # set visible binding to variable
  sampleID <- accn <- recordID <- NULL

  numMined <- refdb %>% dplyr::filter(.data$note == "Mined from GenBank, NCBI") %>% nrow()
  message(paste0("'", numMined, "' records from BOLD were mined from the NCBI.\nIf they are already represented by the NCBI barcodeOre they will be removed to avoid duplicates."))

  # now find all accessionNumbers corresponding to those BOLD processID. These
  # can be found in both "NCBI_ID" and "sampleID" fields, the latter is
  # usually helpful when the former is missing.

  # select accession number in those records. In this code, the coalesce function
  # has an important detail: it substitutes any NA with a value from the other
  # column, but if both columns have a different value (not NA), then the value
  # from the first column will be used instead. The field "NCBI_ID" is thus
  # preferred
  filt_mined <- refdb %>% dplyr::filter(.data$source == "BOLD",
                                        .data$note == "Mined from GenBank, NCBI") %>%
    dplyr::mutate(accn = dplyr::coalesce(.data$NCBI_ID, .data$sampleID)) %>%
    dplyr::select(recordID, accn)

  procs <- lapply(filt_mined$accn, function(sing_accn) {

    if (length(grep(sing_accn, refdb$recordID)) == 1) {

      proc <- filt_mined[filt_mined$accn == sing_accn, ]$recordID
      return(proc)

    } else if (length(grep(sing_accn, refdb$recordID)) > 1) {

      # in this case the NCBI_ID value may be incomplete (for example
      # not have the final .1 and there may be an accession number identical
      # and differing only because ends with a .2) and correspond to multiple
      # recordID from the ncbi barcodeOre. The sampleID may be different
      # and more complete and thus help in this case

      other_accn <- refdb %>% dplyr::filter(.data$source == "BOLD", .data$note == "Mined from GenBank, NCBI") %>%
        dplyr::mutate(accn = dplyr::coalesce(.data$NCBI_ID, .data$sampleID)) %>%
        dplyr::filter(.data$accn == sing_accn) %>%
        dplyr::select(sampleID) %>%
        c(., recursive = TRUE) %>%
        unname()

      if (length(grep(other_accn, refdb$recordID)) == 1) {

        proc <- filt_mined[filt_mined$accn == other_accn, ]$recordID
        return(proc)

      } else if (length(grep(other_accn, refdb$recordID)) > 1) {

        # at this point, the most probable explanation is that the multiple matches
        # derive from multiple copies "or versions" of the same NCBI record and marker
        # This would mean that the recordID is to exclude anyways, as already represented
        # by one or multiple records on the NCBI ore
        proc <- filt_mined[filt_mined$accn == other_accn, ]$recordID
        return(proc)

      }

    }

  })

  # now eliminated the BOLD records that are already represented by NCBI ones
  # and create a new barcode Ore object
  procs_remove <- do.call(c, procs)
  return(procs_remove)

}

#' Remove NCBI records if the corresponding, mined BOLD records are present
#'
#' @param refdb a refdb object
#'
#' @return a refdb object
#'
#' @keywords internal
#' @importFrom rlang .data
#'
ncbiDuplicateCheck <- function(refdb) {

  # set visible binding to variable
  db_xref <- NULL

  numMined <- length(grep("BOLD:", refdb$db_xref))
  message(paste0("'", numMined, "' records from NCBI were mined from BOLD.\nIf they are already represented by the BOLD barcodeOre they will be removed to avoid duplicates.\nDuplicated records obtained from the BOLD will be kept."))

  procs <- refdb[grep("BOLD:", refdb$db_xref), ] %>%
    dplyr::select(db_xref) %>%
    c(., recursive = TRUE) %>% unname() %>%
    stringr::str_remove("\\|.*") %>%
    stringr::str_remove("\\..*") %>%
    stringr::str_remove("BOLD:")

  accns <- lapply(procs, function(pr) {

    if (length(grep(pr, refdb$recordID)) == 1) {

      accn <- refdb[grep(pr, refdb$db_xref), ]$recordID
      return(accn)

    }


  })

  accns_remove <- do.call(c, accns)

}
