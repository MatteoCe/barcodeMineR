## -----------------------------------------------------------------------------
library(barcodeMineR)

# search taxonomic information for a species on the NCBI
tax <- get_ncbi_taxonomy("Dissostichus mawsoni")
tax

## -----------------------------------------------------------------------------
rec_NCBI <- download_ncbi(tax, ask = FALSE)


## ----echo = FALSE-------------------------------------------------------------
rec_NCBI

## -----------------------------------------------------------------------------
# the default filters exclude whole genome shotgun sequences and transcribed 
# shotgun assembly products
rentrez::entrez_search(db="nucleotide", term="(((txid36200[ORGN] NOT wgs[Keyword]) NOT tsa[Keyword]) AND biomol_genomic[PROP]) AND (cds[Feature key] OR rrna[Feature key])")

## -----------------------------------------------------------------------------
rec_NCBI <- download_ncbi(tax)

## -----------------------------------------------------------------------------
# search taxonomic information for a species on BOLD
tax <- get_bold_taxonomy("Dissostichus mawsoni")
tax


## -----------------------------------------------------------------------------
rec_BOLD <- download_bold(tax, ask = FALSE)
rec_BOLD

## -----------------------------------------------------------------------------
total <- mergeBarcodeOres(rec_NCBI, rec_BOLD)

## -----------------------------------------------------------------------------
total

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  # load example datasets from barcodeMineR...
#  rec <- barcodeMineR::example_record
#  seq <- barcodeMineR::example_sequence
#  
#  #... or include the path to your files
#  rec <- "~/home/my/path/data.tsv"
#  seq <- "~/home/my/path/data.fasta"
#  
#  new_records <- loadBarcodeOre(rec, seq)
#  new_records

## ----echo = FALSE, eval = TRUE------------------------------------------------
rec <- barcodeMineR::example_record
seq <- barcodeMineR::example_sequence
new_records <- loadBarcodeOre(rec, seq)
new_records

## -----------------------------------------------------------------------------
total <- mergeBarcodeOres(total, new_records)
total

## -----------------------------------------------------------------------------
library(refdb)

refdb_check_tax_conflict(total)

## -----------------------------------------------------------------------------
total_filt <- refdb_filter_seq_length(total, max_len = 658)

## ----echo = TRUE, eval = TRUE, message = FALSE--------------------------------
refdb_plot_seqlen_hist(total_filt)

## ----echo = TRUE, eval = TRUE-------------------------------------------------
plot_length(total_filt)

## ----echo = TRUE, eval = TRUE-------------------------------------------------
plot_primers(total_filt)

## ----message = FALSE, warning = FALSE-----------------------------------------
refdb_plot_map(total)

