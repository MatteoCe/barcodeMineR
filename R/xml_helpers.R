#' Append root of XML to new xml doc.
#'
#' @param xml an XML object as fetched using ncbi_xml_fetcher
#'
#' @return a new xml object without the original root.
#'
#' @keywords internal
#' @noRd
#'
XML_root <- function(xml) {

  new_xml <- XML::newXMLDoc()

  node_xml <- XML::newXMLNode("Merging", doc = new_xml)

  root_xml <- XML::xmlRoot(xml)

  XML::addChildren(node_xml, root_xml)

  node_xml

}

#' Extract nodes from an NCBI nucleotide xml object.
#'
#' @param featureKey What feature key to extract. Currently, "source", "CDS" and
#' "rRNA" are supported.
#' @param xml An "XMLInternalDocument", "XMLAbstractDocument" class object.
#'
#' @return An xml of 'XMLNodeSet' class.
#'
#' @keywords internal
#' @noRd
#'
#' @description
#' Supported feature keys are "source", "CDS" and "rRNA".
#'
#'
XML_extract_nodes <- function(featureKey, xml) {

  if (!all(featureKey %in% c("source", "CDS", "rRNA"))) {

    stop("Currently the function supports 'source', 'CDS' and 'rRNA' GBFeature keys only.")

  }

  feature <- XML::getNodeSet(
    xml,
    paste0(
      "//GBFeature[GBFeature_key[contains(., '",
      paste(featureKey, collapse = "') or contains(., '"),
      "')]]"
    )
  )

  feature

}

# Very similar to the previous, but works on xml objects from the taxonomy
# database of the NCBI. It is used for the function get_ncbi_taxonomy,
# that might come from a search of species names or a search of taxid coming from
# the download_ncbi_records function. Thus, it has a conditional crossroad, if
# the id provided is a taxid, than that is used to filter the potential multiple
# taxonomies, if it was a species name, than scinetificName is used to filter.

#' Extract nodes from an NCBI taxonomy xml object.
#'
#' @param xml An "XMLInternalDocument", "XMLAbstractDocument" class object.
#'
#' @return An xml of 'XMLNodeSet' class.
#'
#' @keywords internal
#' @noRd
#'
#' @description
#' Similar function to XML_extract_nodes, but works on an NCBI taxonomy xml.
#'
XML_Taxon_extract_nodes <- function(xml) {

  path <- "//TaxaSet/Taxon"

  feature <- XML::getNodeSet(
    xml,
    path
  )

  feature

}

#' Extract GBFeature data from an NCBI xml node object.
#'
#' @param xml An xml of 'XMLNodeSet' class obtained with ncbi_xml_fetcher (db = 'nucleotide').
#' @param path One of 'accession', 'location' or 'qualifier'. This parameter will determine what type of output will be generated.
#' @param quals With path set as 'qualifier' or 'location', determines which qualifiers will be extracted from the xml object. If chosen with path 'location' it must be one of 'from' or 'to'.
#'
#' @return Either a string text with the accession number (path = 'accession') or the location definition (path = 'location'), or a named list with the qualifiers extracted from GBFeature of the NCBI xml object.
#'
#' @keywords internal
#' @noRd
#'
#' @description
#' This function will extract different data from an xml file obtained from the NCBI.
#' An xml corresponding to a particular feature of the entire xml obtained from the
#' NCBI for a specific accession number must be supplied, thus it may be of a GBFeature
#' corresponding to Source, CDS or rRNA. This function accepts three cases (here
#' chosen with the argument "path") "accession" (for the extraction of the accession
#' number), "location" (for the extraction of the position of the sequence in the
#' entire sequence) and "qualifier" (used to extract information from the source key
#' feature). The quals argument is used to specify different metadata the user want
#' to extract from the source feature key of the xml file supplied like "organism",
#' collected_by".
#'
XML_extract <- function(xml, path, quals = NULL) {

  # set cases for this function to work properly
  paths <- c("accession", "location", "qualifier")

  if (!(path %in% paths)) {

    stop("wrong choice for XML_extract path structure")

  } else if (path == "accession") {

    interval_accession <- XML::xpathSApply(
      xml,
      "./GBFeature_intervals/GBInterval/GBInterval_accession",
      XML::xmlValue
    )

    return(interval_accession)

  } else if (path == "location") {

    if (is.null(quals)) {

      pos <- XML::xpathSApply(
        xml,
        "./GBFeature_location",
        XML::xmlValue
      )

    } else if (!(quals %in% c("from", "to"))) {

      stop("a quals in the form of character string 'from' or 'to' should be supplied")

    } else {

      pos <- XML::xpathSApply(
        xml,
        paste0(
          "./GBFeature_intervals/GBInterval//GBInterval_", quals),
        XML::xmlValue
      )

    }

    return(pos)

  } else if (path == "qualifier") {

    values <- sapply(quals, function(x) {
      XML::xpathSApply(
        xml,
        paste0(
          "./GBFeature_quals/GBQualifier[GBQualifier_name='",
          x,
          "']/GBQualifier_value"
        ),
        XML::xmlValue
      )
    }, USE.NAMES=TRUE)

    values_coll <- lapply(values, function(x) {paste(x, collapse = "|")})

    values_coll[nchar(values_coll) == 0] <- NA

    return(values_coll)
  }
}

#' Extract taxonomic data from an NCBI xml node object.
#'
#' @param feature_nodes An xml of 'XMLNodeSet' class obtained with ncbi_xml_fetcher (db = 'taxonomy').
#'
#' @return A named list
#'
#' @keywords internal
#' @noRd
#'
#' @description
#' Works similarly to XML_extract. A named list with taxid, taxonomic rank, scientific
#' name, phylum, class, order, family, genus and species will be returned.
#'
extractTaxonomyTab <- function(feature_nodes) {

  listed_tab <- lapply(feature_nodes, function(feat) {

    taxid <- XML::xpathSApply(
      feat,
      "./TaxId",
      XML::xmlValue)  %>% stats::setNames("taxid")

    scientificName <- XML::xpathSApply(
      feat,
      "./ScientificName",
      XML::xmlValue)  %>% stats::setNames("scientificName")

    rank <- XML::xpathSApply(
      feat,
      "./Rank",
      XML::xmlValue) %>% stats::setNames("rank")

    quals <- c("phylum", "class", "order", "family", "genus", "species")

    values <- sapply(quals, function(x) {
      XML::xpathSApply(
        feat,
        paste0(
          "./LineageEx/Taxon[Rank = '",
          x,
          "']/ScientificName"
        ),
        XML::xmlValue
      )
    }, USE.NAMES=TRUE)

    # as some elements might by empty, turn them to NA or errors might occur
    values_coll <- lapply(values, function(x) {paste(x, collapse = "|")})

    values_coll[nchar(values_coll) == 0] <- NA

    if (rank %in% quals) {

      values_coll[names(values_coll) == rank] <- scientificName

    }

    as.data.frame(c(taxid, rank, scientificName, values_coll))

  })

  taxonomy_tab <- do.call(base::rbind, listed_tab)

  taxonomy_tab
}

#' Extract selection table from a list of xml CDS and rRNA features.
#'
#' @param feature_nodes A list of xml CDS and rRNA features.
#' @param accn a character vector of accession number, those recovered using extractRecordsTab
#'
#' @return A data.frame counting the number of CDS and rRNA features.
#'
#' @keywords internal
#' @noRd
#'
#' @description
#' This function process a list of feature XML objects in order to create the selection tab, the table that will later allow to choose which CDS and rRNA to keep from all of the accession numbers gathered. It works on a lapply function to iteratively extract data.frame-like objects for each CDS and rRNA.
#'
#'
extractSelectionTab <- function(feature_nodes, accn) {

  listed_tab <- purrr::map(feature_nodes, function(feat) {

    # Extract the GBInterval_accession values
    interval_accession <- XML_extract(feat, "accession")

    # for "join" cases, check if with unique is singular, else it may be a sequence
    # at a "remote" location. In that case, choose the first accession and for now
    # 10 november 2023 hope for the better
    if (length(unique(interval_accession)) > 1) {

      interval_accession <- interval_accession[1]

    } else {

      interval_accession <- unique(interval_accession)

    }

    # if the accession number retrieved is not in those obtained from
    # extractRecordsTab, exclude it from the selection table
    if (!(interval_accession %in% accn)) {

      return(NULL)

    }

    # Extract the GBQualifier_value for Qualifier_name 'gene'
    gene_product_values <- XML_extract(feat, "qualifier", c("gene", "product"))

    # extract position coordinates of the corresponding sequence
    location <- XML_extract(feat, "location")

    # create data.frame with the previously extracted informations
    sing_tab <- data.frame(
      GBInterval_accession = interval_accession,
      Gene_value = ifelse(is.null(gene_product_values[["gene"]]),
                          NA,
                          gene_product_values[["gene"]]),
      Product_value = ifelse(is.null(gene_product_values[["product"]]),
                             NA,
                             gene_product_values[["product"]]),
      Feature_coords = ifelse(is.null(location),
                              NA,
                              location
      )
    )

    return(sing_tab)

  }) %>% purrr::compact()

  selection_tab <- do.call(rbind, listed_tab)

  return(selection_tab)

}

#' Extract records table from a list of xml source features.
#'
#' @param feature_nodes A list of xml source features.
#' @param taxonomic_table The taxonomic table obtained from the get_ncbi_taxonomy function.
#'
#' @return A data.frame with source data from each xml source features, thus the metadata reported for each corresponding accession number.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom rlang .data
#'
#' @description
#' Similar to "extractSelectionTab", this function allows to extract record information from a list of XML source features. The taxonomic table, obtained with get_ncbi_taxonomy is here inspected for taxonomic names.
#'
extractRecordsTab <- function(feature_nodes, taxonomic_table) {

  listed_tab <- purrr::map(feature_nodes, function(feat) {

    # Extract the GBInterval_accession values
    interval_accession <- XML_extract(feat, "accession")

    # Extract the GBQualifier_value for the useful "source" qualifier
    pos_from <- XML_extract(feat, "location", "from")
    pos_to <- XML_extract(feat, "location", "to")

    # in case multiple GBInterval qualifiers are available for the same source
    # see HQ615872.1 for example, then multiple positions will be retrieved.
    # remove duplicates and get entire sequence position:
    if (any(duplicated(interval_accession))) {
      interval_accession <- unique(interval_accession)
      pos_from <- pos_from[1]
      pos_to <- pos_to[length(pos_to)]
    }

    # set the qualifiers to extract from the XML to get the records infos
    qualifiers <- c(
      "organism",
      "altitude",
      "identified_by",
      "note",
      "PCR_primers",
      "mol_type",
      "specimen_voucher",
      "collected_by",
      "collection_date",
      "country",
      "lat_lon",
      "db_xref"
    )

    # extract those infos as a list
    source_data <- XML_extract(feat, "qualifier", quals = qualifiers)

    # this gets the taxId corresponding to the record here processed takes the
    # "taxon:xyz" from the source table extracted before
    id <- grep("taxon:", stringr::str_split_1(source_data$db_xref, "\\|"), value = TRUE) %>%
      stringr::str_remove(., "taxon:")

    # if no taxonomic table is provided (download_ncbi with accession numbers)
    # then create a dummy table
    if (is.null(taxonomic_table)) {
      taxonomic_table <- get_ncbi_taxonomy(id, ask = FALSE)
    }

    # if the taxonomic table was obtained using the download_ncbi function,
    # then it may be corresponding to multiple species and thus not work in the
    # following if path. The taxonomic table must then be filtered
    taxonomic_table <- taxonomic_table %>% dplyr::filter(., .data$taxid == id)

    # when the record recovered from the NCBI corresponds to a lower taxonomic
    # rank than those retrieved with get_ncbi_taxonomy (e.g. subspecies excluded
    # by the selection process), skip the recovery of this record
    if (nrow(taxonomic_table) == 0) {
        return(NULL)
    }

    # Complex strings processing:

    ### Geographical coordinates

    # process lat_lon qualifier to extract latitude coordinates
    source_data$lat <- ifelse(is.na(source_data$lat_lon), NA,
                              ifelse(stringr::str_extract(source_data$lat_lon,
                                                          " [S,N] ") == " S ",
                                     stringr::str_extract(source_data$lat_lon,
                                                          "^[0-9]*.[0-9]* [S,N]") %>%
                                       stringr::str_remove(" [S,N]") %>%
                                       stringr::str_replace("^", "-"),
                                     stringr::str_extract(source_data$lat_lon,
                                                          "^[0-9]*.[0-9]* [S,N]") %>%
                                       stringr::str_remove(" [S,N]")))

    # process lat_lon qualifier to extract longitude coordinates
    source_data$lon <- ifelse(is.na(source_data$lat_lon), NA,
                              ifelse(stringr::str_extract(source_data$lat_lon,
                                                          " [W,E]$") == " W",
                                     stringr::str_extract(source_data$lat_lon,
                                                          "[0-9]*[\\.]?[0-9]* [W,E]*$") %>%
                                       stringr::str_remove(" [W,E]$") %>%
                                       stringr::str_replace("^", "-"),
                                     stringr::str_extract(source_data$lat_lon,
                                                          "[0-9]*[\\.]?[0-9]* [W,E]*$") %>%
                                       stringr::str_remove(" [W,E]$")))

    ### Depth or altitude

    # decode depth information. If a negative number is shown, included it as
    # as a new object of the source data list
    source_data$depth <- ifelse(is.na(source_data$altitude), NA,
                                ifelse(length(grep("-", source_data$altitude)) == 0, NA,
                                       paste0(stringr::str_remove_all(source_data$altitude, " m"))))

    # decode alrtitude information. If a positive number is shown, included it
    # as a new object of the source data list. The altitude must be done after
    # depth as it would overwrite the original info
    source_data$altitude <- ifelse(is.na(source_data$altitude), NA,
                                   ifelse(length(grep("-", source_data$altitude)) > 0, NA,
                                          paste0(stringr::str_remove_all(source_data$altitude, " m"))))

    ### Primers information

    # decode primer's direction
    source_data$directionPrimers <- ifelse(length(grep("name", source_data$PCR_primers, value = TRUE)) == 0, NA, "F|R")

    # primer's names will be obtained from the specific naming included in the INSDC
    source_data$PCR_primers <- ifelse(length(grep("name", source_data$PCR_primers, value = TRUE)) == 0, NA,
                                      paste(paste(stringr::str_extract_all(source_data$PCR_primers, "fwd_name: [[:graph:]]*,")[[1]] %>%
                                                    stringr::str_remove_all("fwd_name: ") %>%
                                                    stringr::str_remove_all(","), collapse = ","),
                                            paste(stringr::str_extract_all(source_data$PCR_primers, "rev_name: [[:graph:]]*,")[[1]] %>%
                                                    stringr::str_remove_all("rev_name: ") %>%
                                                    stringr::str_remove_all(","), collapse = ","), sep = "|"))

    ### Length of source

    # Calculate the length of the entire source, the original fasta sequence
    # containing all the CDS and rRNA of each accession
    source_data$sourceLength <- ifelse((is.null(pos_from) | is.null(pos_to)),
                                       NA,
                                       length(seq(pos_from, pos_to, 1)))

    # create data.frame with the previously extracted information
    data.frame(
      source = "NCBI",
      sourceID = interval_accession,
      DNA_seq = NA,
      markerCode = NA,
      lengthGene = NA,
      sampleID = source_data$specimen_voucher,
      QueryName = taxonomic_table$queryName,
      phylum = taxonomic_table$phylum,
      class = taxonomic_table$class,
      order = taxonomic_table$order,
      family = taxonomic_table$family,
      genus = taxonomic_table$genus,
      species = taxonomic_table$species,
      identified_by = source_data$identified_by,
      taxNotes = NA,
      db_xref = source_data$db_xref,
      NCBI_ID = interval_accession,
      institutionStoring = NA,
      collected_by = source_data$collected_by,
      collection_date = source_data$collection_date,
      lat = source_data$lat,
      lon = source_data$lon,
      altitude = source_data$altitude,
      depth = source_data$depth,
      country = source_data$country,
      directionPrimers = source_data$directionPrimers,
      lengthSource = source_data$sourceLength,
      PCR_primers = source_data$PCR_primers,
      note = source_data$note
    )

  }) %>% purrr::compact()

  # merge all data.frames into a single records tab
  records_tab <- do.call(base::rbind, listed_tab)

  records_tab

}
