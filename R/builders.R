#' Clean and process chosen sequences based on the NCBI location qualifiers.
#'
#' @param accn a string in the form "accession number|marked code".
#' @param DNAString the DNA sequences in the format "DNAStringSet"
#' @param selection_tab the selection tab created by "extractSelectionTab"
#' @param skip.unknown.pos skip the definition of start and end position of
#'   sequence to extract based on ">" and "<" symbols as first and last position
#'   of entire sequence
#'
#' @return a list of "DNAStringSet" objects corresponding to the DNA sequences trimmed as indicated in the relative location qualifiers in the selection tab.
#'
#' @keywords internal
#' @noRd
#'
#' @description
#' To be used in a vectorized loop (like lapply) for all "accession number|marker code" chosen in the selection process. This function will select the correct sequences from the total obtained in select_ncbi_genes.R or download_ncbi_genes.R and, if necessary, cut the correct region of the sequence corresponding to the chosen CDS or rRNA. It is also used inside the "select_bold_genes" function. In this case, this script will simply select and extract the "accn" previously selected from the whole sequences gathered at the beginning of the selection script and remove gaps, if present.
#'
buildSequences <- function(accn, DNAString, selection_tab, skip.unknown.pos = TRUE) {

  if (length(grep("BOLD", colnames(selection_tab)[1], value = FALSE)) == 1) {

    DNAString[[accn]] <- Biostrings::DNAString(gsub("-", "", DNAString[accn]))
    return(DNAString[accn])

  } else {

    # extract gene_name and accession number from input accn
    accession <- stringr::str_remove(string = accn, pattern = "\\|.*")
    gene_name <- stringr::str_remove(string = accn, pattern = paste0(stringr::str_replace(accession, "\\.", "\\\\."), "\\|"))

    # it is possible that a single records on the NCBI contains multiple copies
    # or "variants" of the same CDS (see https://doi.org/10.1093/mollus/eyw052)
    # as in the case of the NCBI accession number "NC_016423.1". it is thus
    # necessary to differentiate the two variants in two different sequences, thus
    # records. A for loop would suffice as in any case, the copies would be limited
    # to only one or possibly just a few, and, if there are no copies at all, it
    # would not decrease the speed of the function.

    feat <- selection_tab[(selection_tab$GBInterval_accession %in% accession) &
                            (selection_tab$Gene_value %in% gene_name | selection_tab$Product_value %in% gene_name), ]

    ncopies <- nrow(feat)

    versions <- list()

    for (copy in seq(1:ncopies)) {

      if (ncopies > 1) {

        # if multiple copies are present, add a "_copy number" to the accession name
        # to discriminate the two in the final renaming of the sequences
        name <- paste(paste0(accession, "_", copy), gene_name, sep = "|")

      } else {

        name <- accn

      }

      # extract the single table
      sing_feat <- feat[copy, ]

      # extract Location qualifier
      location <- sing_feat$Feature_coords

      if (

        # check if it corresponds to a single, specific base in the sequence or...
        stringr::str_like(location, "^[0-9]+$") ||

        # ... a single base of unknown position or ...
        stringr::str_detect(location, "[0-9]+\\.[0-9]+") ||

        stringr::str_detect(location, "^\\<[0-9]+$") ||

        stringr::str_detect(location, "^\\>[0-9]+$") ||

        # ... a site between two specific positions or ...
        stringr::str_detect(location, "[0-9]+\\^[0-9]+") ||

        # ... to a site, base or sequence span in a remote entry
        stringr::str_detect(location, "[0-9A-Z]+\\.[0-9]:")

      ) {

        # left commented in case it will be needed in future {
          # extract whole sequence
          # seq <- DNAString[accession]

          # update name of the sequence with accn (thus "accession number|gene name")
          # seq@ranges@NAMES <- name
          # versions[[copy]] <- seq

        # }

        # return nothing
        return(NULL)

      } else {

        # extracts ranges of bases to extract from whole sequence
        # old regex changed the 27 Jun 2024 ([<0-9]+\\.\\.[>0-9]+)+|([<>0-9]+)+
        bases_ranges <- unlist(stringr::str_extract_all(location, "(([<0-9]+\\.\\.[>0-9]+)+|([<>0-9]+)+(?=[\\)|,]))+"))

        # create list with subsequences from whole sequence based on the range of bases extracted before
        subseqs <- sapply(bases_ranges, function(range) {

          # get start coord
          start <- stringr::str_replace(range, "(.+)\\.\\.(.+)", "\\1")
          # get end coord
          end <- stringr::str_replace(range, "(.+)\\.\\.(.+)", "\\2")

          # based on commit of 25 June 2024, skip using start and end position of
          # entire sequence as substitutes of unknown positions. Keep the code
          # for possible future use
          if (!skip.unknown.pos) {

            if (stringr::str_detect(start, "<")) {

              # if the symbol "<" is present, get the first position of the whole sequence
              start <- 1

            } else {

              start <- as.integer(start)

            }

            if (stringr::str_detect(end, ">")) {

              # if the symbol ">" is present, get the last position of the whole sequence
              end <- DNAString[accession]@ranges@width

            } else {

              end <- as.integer(end)

            }

          } else {

            if (all(start == end)) {
              end <- start <- as.integer(stringr::str_remove(start, "\\<|\\>"))
            } else {
              start <- as.integer(stringr::str_remove(start, "\\<"))
              end <- as.integer(stringr::str_remove(end, "\\>"))
            }

          }

          # extract subsequence
          subseq <- XVector::subseq(DNAString[accession], start, end)

          return(subseq)

        }, USE.NAMES=TRUE)

        if (

          # if the location information is not a simple "XXX..XXX" condition, meaning
          # that it corresponds to a simple subsequence of the whole...
          !(all(bases_ranges == location))

        ) {

          # extracts groups to operate directly on (join or complement)
          # old regex changed the 27 Jun 2024 [a-z]+\\((([<0-9]+\\.\\.[>0-9]+)+,?)+([<0-9]+\\.\\.[>0-9]+)?\\)
          operations_on_groups <- unlist(stringr::str_extract_all(location, "[a-z]+\\(((([<0-9]+\\.\\.[>0-9]+)|([<>]?[0-9]+))+,?)+([<0-9]+\\.\\.[>0-9]+)?\\)"))

          processed_subseqs <- sapply(operations_on_groups, function(group_operation) {

            # extract exact operation to perform
            operation <- stringr::str_extract(group_operation, "^[a-z]+")

            # get ranges of bases previously extracted to operate on
            range <- unlist(stringr::str_extract_all(group_operation, "([\\<0-9]+\\.\\.[\\>0-9]+)|([<>]?[0-9]+)"))

            if (operation == "complement") {

              # extract sequence and perform operation
              subseq <- Biostrings::reverseComplement(subseqs[[range]])

              return(subseq)

            } else if (operation == "join") {

              # extract sequence and perform operation
              subseq <- Biostrings::DNAStringSet(
                do.call(c, sapply(unname(subseqs[range]),
                                  function(x) {
                                    return(x[[1]])
                                    })))

              return(subseq)

            }

          }, USE.NAMES = TRUE)

          # added filter for cases where a join operation is performed on a
          # subsequence that has been processed while another was not. The
          # example can be found for NC_082072.1 having two copies of "psbA"
          # one of which (the second) has this location:
          # join(83517..83856,complement(238..959))
          # This will add the base ranges not processed to the new "processed_subseqs"
          # object in order to detect the next operation to perform
          if (!all(sapply(names(subseqs), function(pat) {any(stringr::str_detect(names(processed_subseqs), pat))}, simplify = TRUE, USE.NAMES = FALSE))) {

            unoperated <- which(!sapply(names(subseqs), function(pat) {any(stringr::str_detect(names(processed_subseqs), pat))}, simplify = TRUE, USE.NAMES = TRUE)) %>% names()

            processed_subseqs <- c(subseqs[unoperated], processed_subseqs)

          }

          # now the fun part.

          re_processed_subseqs <- processed_subseqs

          repeat {

            if (all(names(re_processed_subseqs) == location)) {

              break

            }

            # set the conditional OR names of each element of the processed_subseqs list. This
            # will be used by stringr to search all patterns corresponding to those names. Here, all
            # alphanumeric and blank characters, except the comma, will be appended with a "\\" symbology
            # that allows stringr to consider them as literal and no special characters (like "(" for example)
            or_names_coll <- paste(stringr::str_replace_all(names(re_processed_subseqs), "([^[:alnum:][:space:],])", "\\\\\\1"), collapse = "|")

            # get group operations including the sequences previously processed to operate on again
            re_operations_on_groups <- unlist(stringr::str_extract_all(location, paste0("[a-z]+\\(", paste0("((", or_names_coll, ")+,?)+(", or_names_coll, ")?"),"\\)")))

            re_processed_subseqs_repeat <- sapply(re_operations_on_groups, function(group_operation) {

              # extract exact operation to perform
              operation <- stringr::str_extract(group_operation, "^[a-z]+")

              # get range of bases previously extracted to operate on
              # use a command similar to the previous one creating the object or_names to
              # create a query for stringr to search which elements from the list re_processed_subseqs
              # should be extracted for the operation to perform
              or_names <- stringr::str_replace_all(names(re_processed_subseqs), "([^[:alnum:][:space:],])", "\\\\\\1")

              range <- names(re_processed_subseqs[stringr::str_detect(group_operation, or_names)])

              if (operation == "complement") {

                # extract sequence and perform operation
                subseq <- Biostrings::reverseComplement(re_processed_subseqs[[range]])

                return(subseq)

              } else if (operation == "join") {

                # extract sequence and perform operation
                subseq <- Biostrings::DNAStringSet(
                  do.call(c, sapply(unname(re_processed_subseqs[range]),
                                    function(x) {
                                      return(x[[1]])
                                      })))

                return(subseq)

              }

            })

            re_processed_subseqs <- re_processed_subseqs_repeat

          }

          # extract sequence
          seq <- unname(re_processed_subseqs[location])[[1]]

          # update name of the sequence with accn (thus "accession number|gene name")
          seq@ranges@NAMES <- name
          versions[[copy]] <- seq

        } else {

          # if the location information is a simple "XXX..XXX" condition
          # extract whole sequence
          seq <- unname(subseqs[bases_ranges])[[1]]

          # update name of the sequence with accn (thus "accession number|gene name")
          seq@ranges@NAMES <- name
          versions[[copy]] <- seq

        }

      }

    }

    # merge all sequences in the versions list, in case multiples copies are present
    # also, remove duplicates, as sometimes the different copies might actually be
    # the same sequence, or the initial and final positions of the copies are
    # unknown (for example copy1: <1..382>, copy2: <383..745>), thus each copy will
    # actually retain the whole sequence
    seqs <- do.call(c, versions) %>% unique()

    # if the final sequence, originated from multiple "copies" of duplicate sequences)
    # correspond to one sequence only, confirm original name as accession/gene_name
    if (all(length(versions) > 1 & length(seqs) == 1)) {
      seqs@ranges@NAMES <- accn
    }

    return(seqs)

  }}

#' Update record information with DNA sequence data coming from the buildSequences script.
#'
#' @param accn a string in the form "accession number|marked code".
#' @param records the record table.
#' @param sequences the DNA sequences in the format "DNAStringSet".
#'
#' @return a data.frame, with information updated based on the additional information from the "built" sequences.
#'
#' @keywords internal
#' @noRd
#'
#' @description
#' To be used after buildSequences as it calculates the final length of each marker. It only extracts the info of each chose accession number and marker combination and adds the gene_name and final length of each sequence. There is a difference depending on the origin of the records, whether they come from a search on BOLD records or NCBI records. In the former, the information obtained from the fetching of the sample data may not be unequivocally selected using the sourceID, especially for those cases where the sourceID is the same for different sequences (16S and COI CDS coming from same source sequence on the NCBI). In this case, records must be selected by sourceID and markerCode altogether as it is already available in the records table. Instead, for the NCBI the records can only be selected by sourceID, and, considering the rest of the information is unique, possible duplicate rows of the selected data.frame can be removed to only have a dataframe with a single row.
#'
buildRecord <- function(accn, records, sequences) {

  # extract gene_name and record name from input accn
  gene_name <- stringr::str_remove(string = accn, pattern = ".*\\|")
  accession <- stringr::str_remove(string = accn, pattern = "\\|.*")

  # extract specific record
  if (unique(records$source) == "BOLD") {
    record <- records[((records$sourceID %in% accession) & (records$markerCode %in% gene_name)), ]

    # update "markerCode" and "lengthGene" fields
    record$markerCode <- gene_name
    record$lengthGene <- sequences[accn]@ranges@width
    record$DNA_seq <- as.character(sequences[[accn]])

  } else {

    # extract single record from records_tab, which is always single (corresponds
    # to the entire feature)
    record <- unique(records[(records$sourceID %in% accession), ])

    if (accn %in% sequences@ranges@NAMES) {

      # update "processID", "markerCode" and "lengthGene" fields
      record$markerCode <- gene_name
      record$lengthGene <- sequences[accn]@ranges@width
      record$DNA_seq <- as.character(sequences[[accn]])

    } else {

      # in case of multiple copies of the same gene, run it in a for loop
      versions <- list()
      pos <- 0

      # xxx
      copies <- purrr::keep(stringr::str_split(sequences@ranges@NAMES, "\\|"), function(name) {

        stringr::str_detect(name[1],  accession) && (name[2] == gene_name)

      }) %>%

        purrr::map_chr(., function(name) {

          stringr::str_c(name[1], name[2], sep = "|")

        })

      for (copy in copies) {

        pos <- pos + 1

        # update "sourceID", "markerCode" and "lengthGene" fields
        record$sourceID <- stringr::str_remove(string = copy, pattern = "\\|.*")
        record$markerCode <- gene_name
        record$lengthGene <- sequences[copy]@ranges@width
        record$DNA_seq <- as.character(sequences[[copy]])

        versions[[pos]] <- record
      }

      record <- do.call(rbind, versions)

    }

  }

  record

}

#' Create a refdb from the final records data.frame
#'
#' @param records the records data.frame, as obtained by the "download" functions.
#' @param prefix the prefix to use in case custom IDs are chosen by the user instead of the original processID/accession numbers
#'
#' @return a refdb data.frame
#'
#' @keywords internal
#' @noRd
#'
buildBarcodeOre <- function(records, prefix = NULL) {

  # create tibble data.frame
  records <- tibble::as_tibble(records)

  # group taxonomy fields
  taxonomy <- c(phylum = "phylum",
                class = "class",
                order = "order",
                family = "family",
                genus = "genus",
                species = "species")

  # create custom IDs prefix if chosen or not by the user
  if (is.null(prefix)) {

    records$recordID <- records$sourceID

  } else {

    records$recordID <- stringr::str_c(prefix, seq_along(records$sourceID))

  }

  # set correct class for latitude and longitude values
  records$lat <- as.numeric(records$lat)
  records$lon <- as.numeric(records$lon)

  # define additional fields as a character vector
  additional_data <- c("lengthGene",
                       "sampleID",
                       "QueryName",
                       "identified_by",
                       "taxNotes",
                       "db_xref",
                       "sourceID",
                       "NCBI_ID",
                       "institutionStoring",
                       "collected_by",
                       "collection_date",
                       "altitude",
                       "depth",
                       "country",
                       "directionPrimers",
                       "lengthSource",
                       "PCR_primers",
                       "note")

  # create refdb database
  db <- refdb::refdb_set_fields(records,
                                id = "recordID",
                                marker = "markerCode",
                                sequence = "DNA_seq",
                                taxonomy = taxonomy,
                                source = "source",
                                latitude = "lat",
                                longitude = "lon",
                                reference = additional_data)

  # update order of columns
  db %>% dplyr::relocate(dplyr::any_of(c("recordID",
                                         "markerCode",
                                         "DNA_seq",
                                         unname(taxonomy),
                                         "source",
                                         "lat",
                                         "lon",
                                         additional_data
                                         )))

}
