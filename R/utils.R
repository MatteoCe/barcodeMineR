#' Ask the user to interactively chose which values of a tables' column to retain from a table.
#'
#' @param df The data.frame that will be filtered.
#' @param field The column name to extract the values from. It must be a character string.
#'
#' @return A filtered data.frame.
#'
#' @keywords internal
#' @importFrom rlang .data
#'
#' @description
#' This function is used to let a user filter a data.frame based on the values of a column
#'
ask_user <- function(df, field, ask = TRUE) {

  table <- dplyr::filter(df, !is.null(df[[field]]) & !is.na(df[[field]])) %>%
    dplyr::mutate(dplyr::across({{field}}, as.factor)) %>%
    dplyr::count(dplyr::across({{field}}), sort = TRUE) %>%
    tibble::as_tibble()

  if (ask) {

    if (field == "rank") {

      message("\n")
      print(table, n = nrow(table))

      message(
        "\nThe names searched retrieved info for multiple taxonomic levels. Indicate which levels to keep.
        (0 will return none, other inputs will return all of them)\n"
      )

    } else {

      message("\n")
      table <- table %>% dplyr::mutate(dplyr::across({{field}}, as.character)) %>%
        dplyr::arrange(., {{field}}, .locale = "en")

      print(table, n = nrow(table))

      message(
        "\nIndicate the row numbers corresponding to the CDS and rRNA of interest.
        (0 will return none, other inputs will return all of them)\n"
      )
    }

    progressr::without_progress({

      numb <- scan(n = nrow(table),
                   quiet = TRUE,
                   what = 'raw')

    })

  } else {

    numb <- NULL

  }

  if (!(0 %in% numb) & !is.null(numb) & length(numb) %in% seq(1:nrow(table)) &
      all(numb %in% c(1:nrow(table)))) {

    filt_table <- dplyr::filter(df, df[[field]] %in% as.character(table[numb, ][, 1, drop = TRUE]))
    return(filt_table)

  } else if (0 %in% numb) {
    return(NULL)

  } else {
    return(df)

  }

}

#' Filter a selection table interactively to return a vector of accession number concatenated with a gene or product name
#'
#' @param selection_tab A selection table, as produced by the ::: function
#'
#' @return A vector of the type: "accession number|gene or product name"
#'
#' @keywords internal
#'
#' @description
#' It will take the tab with the basin info from each accession number (CDS and/or rRNA name, and coordinates of the sequence) and output a character vector with the accNum and the feature name (gene or product) separated by a pipe.
#'
#'
select_accessions <- function(selection_tab, ask = TRUE) {

  # ask user to select the Gene names
  filt_tab <- ask_user(selection_tab, "Gene_value", ask = ask)

  # check the correct path:
  # (1) the user may have chosen a group or a single gene name that wants to keep
  # (2) the user do not want any of the shown gene names

  if (!is.null(filt_tab)) {
    # Case (1) no 0 values have been set by the user, the input is not null, the
    # number of input values does not exceed the possible choices and the input
    # values themselves are included in the choices available

    names_gene <- stringr::str_c(filt_tab[, 1],
                                 filt_tab$Gene_value,
                                 sep = "|") %>% .[!is.na(.)]

  } else {
    # Case (2) "0" has been chosen, thus none of the available choices have been
    # selected. Returning no accession numbers from this selection

    # set no accession number as output
    names_gene <- NULL

    # moreover, the original filt_tab will be kept as a data.frame with no rows
    # this will allow to keep the information on how many nomenclatures for the gene
    # (gene/product qualifier) are available

    filt_tab <- selection_tab[NULL,]

  }

  # Are there other nomenclatures to filter from? If filt_tab has only 2, then
  # it means there is only one nomenclature to choose from (the "gene" qualifier
  # but not the "product", for example)
  if (ncol(filt_tab) > 2) {

    # filter out the gene names previously selected from the tax_tab and count
    # the number of occurrences of each product name. This will be prompt to the user
    # which will have to select the product of interest
    selection_tab_2 <- selection_tab[!(selection_tab$Gene_value %in% unique(stringr::str_remove_all(names_gene, ".*\\|"))), ]

    if (nrow(selection_tab_2) > 0) {

      filt_tab_2 <- ask_user(selection_tab_2, "Product_value", ask = ask)

      # check the correct path, as before
      if (!is.null(filt_tab_2)) {

        names_product <- stringr::str_c(filt_tab_2$GBInterval_accession, filt_tab_2$Product_value, sep = "|")

      } else {

        names_product <- NULL

      }

      names <- c(names_gene, names_product)

    } else {

      names <- names_gene

    }

  } else {

    names <- names_gene

  }

  if (is.null(names)) {

    stop("No gene or product qualifier has been chosen, no records will be returned")

  } else {

  return(names)

  }

}

#' Calculate the right retstart and retmax for the NCBI fetcher functions.
#'
#' @param split The starting position.
#' @param hits The total number of ids from the web_history object that will be used to download data from the NCBI.
#' @param hits_x_iteration The number of ids to fetch each execution time.
#'
#' @keywords internal
#'
#' @return A list with two values: the retstart parameter and the retmax. The last one is calculated based on the total number of ids in the web_history object and the starting position.
#'
web_history_parameter <- function(split, hits, hits_x_iteration) {

  parameters <- list()

  if (hits < 2) {

    parameters[1] <- NULL
    parameters[2] <- 1

  } else if ((hits - split) < hits_x_iteration) {

    parameters[1] <- split
    parameters[2] <- (hits - split)

  } else {

    parameters[1] <- split
    parameters[2] <- hits_x_iteration

  }

  parameters

}

#' Get all the retstart positions based on total number of counts and xml rate of a web_history object.
#'
#' @param counts total number of records corresponding to a web_istory object
#' @param rate the number of record to download at eahc iteration of the ncbi_xml_fetcher function.
#'
#' @return a numeric vector with each retstart position for a particular
#'
#' @keywords internal
#'
web_history_splitter <- function(counts, rate) {

  splits <- seq(0, counts, rate)

  # if one of the starting numbers coincide with the total number of hits, remove that number from the sequence
  if ((counts %in% splits) & (counts != 1)) {

    splits <- splits[-length(splits)]

  }

  splits

}

#' Choose the optimal http requests limit based on the presence or absence of an NCBI API key.
#'
#' @param api_rate the api rate limit, as set or not by the user
#'
#' @return A number, defining the correct api rate.
#'
#' @keywords internal
#'
set_ncbi_rate <- function(api_rate) {

  # define default values if api_rate has not been manually set
  if (is.null(api_rate)) {

    # get entrez key from the environment variables
    key <- Sys.getenv("ENTREZ_KEY")

    if (nchar(key) == 0) {

      rate_ncbi <- 3

    } else {

      rate_ncbi <- 9

    }

  # if the rate has been manually defined, check whether it has correct values
  } else if (!(api_rate %in% seq(1, 9, 0.1))) {

    stop("The api rate parameter must be a number between 1.0 and 9.0")

  } else {

    rate_ncbi <- api_rate

  }

  return(rate_ncbi)

}

#' Is the vector a species or a taxid.
#'
#' @param ids The character wtring with a taxonomic name or a taxid.
#'
#' @keywords internal
#'
#' @return A condition, either TRUE or FALSE.
#'
set_id_type <- function(ids) {

  if (any(base::grepl("\\D", ids))) {

    areSpecies <- TRUE

  } else {

    areSpecies <- FALSE

  }

  areSpecies

}

#' Get lower taxonomic level from input
#'
#' @param tax
#'
#' @keywords internal
#'
#' @return character vector
#'
get_lower_tax_rank <- function(tax) {

  if (tax == "species") {

    return(NULL)

  }

  tax_ranks <- c('superkingdom',
                 'kingdom',
                 'subkingdom',
                 'infrakingdom',
                 'phylum',
                 'division',
                 'subphylum',
                 'subdivision',
                 'infradivision',
                 'superclass',
                 'class',
                 'subclass',
                 'infraclass',
                 'superorder',
                 'order',
                 'suborder',
                 'infraorder',
                 'superfamily',
                 'family',
                 'subfamily',
                 'tribe',
                 'subtribe',
                 'genus',
                 'subgenus',
                 'section',
                 'subsection',
                 'species group',
                 'species',
                 'subspecies',
                 'variety',
                 'form',
                 'subvariety',
                 'race',
                 'stirp',
                 'morph',
                 'aberration',
                 'subform',
                 'unspecified',
                 'no rank')

  accepted_tax_ranks <- c('kingdom',
                          'phylum',
                          'class',
                          'order',
                          'family',
                          'subfamily',
                          'genus',
                          'species')

  pos <- grep(paste0("\\b", tax,"\\b"), tax_ranks)
  val_tax <- tax_ranks[pos + 1]

  while (!(val_tax %in% accepted_tax_ranks)) {

    pos <- pos + 1

    val_tax <- tax_ranks[pos]

  }

  val_tax

}

#' Turn a data.frame or fasta text object of sequences into DNAStringSet
#'
#' @param sequences Object containing sequences as data.frame or fasta text object
#'
#' @return DNAStringSet
#'
#' @keywords internal
#'
textToDNAStringSet <- function(sequences) {

  if (inherits(sequences, "data.frame")) {

    sequences <- paste(paste(
      paste0(">",
             paste0(sequences$sourceID,
                    "|",
                    sequences$markerCode))),
      sequences$DNA_seq,
      sep = "\n") %>% stringr::str_flatten(., collapse = "\n\n")

  }

  fasta_list <- seqinr::read.fasta(file = textConnection(sequences),
                                   as.string = TRUE)

  # update as DNAStringSet object from Biostrings
  vapply(fasta_list, function(single_fasta) {

    paste(seqinr::getSequence(single_fasta), collapse = "")

    }, character(1)) %>% as(., "DNAStringSet")

}
