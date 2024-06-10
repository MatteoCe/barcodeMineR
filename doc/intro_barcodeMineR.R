## -----------------------------------------------------------------------------
library(barcodeMineR)

# search taxonomic information for a species on the NCBI
tax <- get_ncbi_taxonomy("Dissostichus mawsoni")
tax

## -----------------------------------------------------------------------------
rec_NCBI <- download_ncbi(tax, ask = FALSE)


## ----echo = FALSE, eval = TRUE------------------------------------------------
rec_NCBI

## ----eval = TRUE--------------------------------------------------------------
# the default filters exclude whole genome shotgun sequences and transcribed shotgun assembly products
rentrez::entrez_search(db="nucleotide", term="(((txid36200[ORGN] NOT wgs[Keyword]) NOT tsa[Keyword]) AND biomol_genomic[PROP]) AND (cds[Feature key] OR rrna[Feature key])")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  rec_NCBI <- download_ncbi(tax)

## -----------------------------------------------------------------------------
# search taxonomic information for a species on BOLD
tax <- get_bold_taxonomy("Dissostichus mawsoni")
tax


## -----------------------------------------------------------------------------
rec_BOLD <- download_bold(tax, ask = FALSE)
rec_BOLD

## -----------------------------------------------------------------------------
total <- mergeBarcodeOres(rec_NCBI, rec_BOLD)
total

## -----------------------------------------------------------------------------
library(refdb)

refdb_check_tax_conflict(total)

## ----message = FALSE, warning = FALSE-----------------------------------------
refdb_plot_map(total)

## -----------------------------------------------------------------------------


