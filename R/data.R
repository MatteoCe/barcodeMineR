#' Example records data frame
#'
#' An example data frame, including meta-data for a fake sequence used as example
#' dataset for the function loadBarcodeOre or for the vignettes.
#'
#' @format ## `example_record`
#' A data frame with 1 row and 25 columns:
#' \describe{
#'   \item{sourceID}{Record code}
#'   \item{source}{Acronym/name of the database/department from which the record/sequence originates}
#'   \item{phylum, class, order, family, genus, species}{Taxonomic classification of the record}
#'   \item{lengthSource}{Length of the sequence, in base pairs}
#'   \item{sampleID}{Sample code, corresponding to the specimen from which the sequence originates}
#'   \item{identified_by, collected_by, taxNotes, collection_date, note}{Character class additional informations on the sample's origin}
#'   \item{db_xref}{Generally unused for user's records, include the taxid of the NCBI record (taxid:"00001") separated by a pipe for the corresponding BOLD processID}
#'   \item{NCBI_ID}{Accession number of the record on the NCBI database}
#'   \item{institutionStoring}{Acronym/name of the institution storing the sample}
#'   \item{altitude, depth}{Elevation data}
#'   \item{country}{Name of the locality from which the sample originates}
#'   \item{lat, lon}{Geographical coordinates, in decimal degrees}
#'   \item{directionPrimers}{Gives the type of primers that are shown in the foeld "PCR_primers": "F|R" indicates that the primers' names in the "PCR_primers" separated by a pipe correspond to the forward and revers primers in this order}
#'   \item{PCR_primers}{Names of the primers used to amplify the sequence of the record}
#'   \item{note}{Additional notes}
#'   ...
#' }
#' @source generated manually, consult data_raw directory for details
"example_record"

#' Example DNAStringSet sequence
#'
#' An example sequence, corresponding to the data in example_record. Used as
#' example dataset for the function loadBarcodeOre or for the vignettes.
#'
#' @format ## `example_sequence`
#' A DNAStringSet object (from Biostrings) of 658 base pairs of "width". The
#' item NAMES includes the sourceID and the markerCode separated by a pipe.
#'
#' @source generated manually, consult data_raw directory for details
"example_sequence"
