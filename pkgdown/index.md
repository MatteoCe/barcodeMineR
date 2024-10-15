
<!-- README.md is generated from README.Rmd. Please edit that file -->

**NOTE: This package is still in development and awaiting for CRAN
approval.**

# barcodeMineR <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/MatteoCe/barcodeMineR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/MatteoCe/barcodeMineR/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/MatteoCe/barcodeMineR/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/MatteoCe/barcodeMineR/actions/workflows/test-coverage.yaml)
[![codecov](https://codecov.io/gh/MatteoCe/barcodeMineR/graph/badge.svg?token=62OUVOL8MP)](https://app.codecov.io/gh/MatteoCe/barcodeMineR)
<!-- badges: end -->

## The barcodeMineR package

This package allows to query multiple taxonomic names on the NCBI and
BOLD repositories and retrieve DNA barcodes and associated metadata for
any wanted marker. It heavily relies on the
[bold](https://github.com/ropensci/bold) and
[rentrez](https://github.com/ropensci) packages from the
[rOpenSci](https://ropensci.org/), and takes advantage of the
asynchronous framework from the
[future](https://github.com/HenrikBengtsson/future) package to speed up
the retrieval of data from the NCBI, respecting its API requests rate
limit.

The final output is a data frame object modified following the
formatting requirements of the [refdb](https://github.com/fkeck/refdb)
package, and it is cleaned from mining duplicates, differences in
formats of BOLD and NCBI metadata (*e.g.* geographic coordinates, dates)
and commonly occurring issues that derive from downloading and merging
data from those two repositories.

In synthesis, it provides a unified framework for mining DNA Barcodes
from the main online genomic repositories, providing clean,
metadata-rich sequences in a programmatic way or interactively,
depending on the wanted usage.

## Installation

The development version of this package can be installed directly from
this GitHub repository, using the
[Bioconductor](https://www.bioconductor.org/install/) package:

``` r
BiocManager::install("MatteoCe/barcodeMineR")
```

The last stable version is available through CRAN:

``` r
install.packages("barcodeMineR")
```

## Basic usage

The basic functioning of the package includes only two steps, which must
be run consequently:

``` r
library(barcodeMineR)

# check if a species is on the NCBI nucleotide database:
tax <- get_ncbi_taxonomy("Dissostichus mawsoni")

# download all records of this species from the database:
rec <- download_ncbi(tax, ask = FALSE)

# display output:
rec
#> # A tibble: 192 × 30
#>    recordID   markerCode DNA_seq  phylum class order family genus species source
#>    <chr>      <chr>      <DNA>    <chr>  <chr> <chr> <chr>  <chr> <chr>   <chr> 
#>  1 HM422302.1 COI        CTCTACT… Chord… Acti… Perc… Notot… Diss… Dissos… NCBI  
#>  2 ON000293.1 COX1       GCGCCTG… Chord… Acti… Perc… Notot… Diss… Dissos… NCBI  
#>  3 MK843765.1 COI        TCTCTAC… Chord… Acti… Perc… Notot… Diss… Dissos… NCBI  
#>  4 DQ498816.1 Cytb       GCCACCC… Chord… Acti… Perc… Notot… Diss… Dissos… NCBI  
#>  5 DQ498794.1 Rhod       GCCTACA… Chord… Acti… Perc… Notot… Diss… Dissos… NCBI  
#>  6 MK500763.1 enc1       TCTGACG… Chord… Acti… Perc… Notot… Diss… Dissos… NCBI  
#>  7 MG729451.1 COI        GAACTTA… Chord… Acti… Perc… Notot… Diss… Dissos… NCBI  
#>  8 KY656477.1 COI        GCCGGAA… Chord… Acti… Perc… Notot… Diss… Dissos… NCBI  
#>  9 LC138011.1 ND1        ATGCTTT… Chord… Acti… Perc… Notot… Diss… Dissos… NCBI  
#> 10 LC138011.1 ND2        ATGAGCC… Chord… Acti… Perc… Notot… Diss… Dissos… NCBI  
#> # ℹ 182 more rows
#> # ℹ 20 more variables: lat <dbl>, lon <dbl>, lengthGene <int>, sampleID <chr>,
#> #   QueryName <chr>, identified_by <chr>, taxNotes <lgl>, db_xref <chr>,
#> #   sourceID <chr>, NCBI_ID <chr>, institutionStoring <lgl>,
#> #   collected_by <chr>, collection_date <chr>, altitude <lgl>, depth <lgl>,
#> #   country <lgl>, directionPrimers <chr>, lengthSource <int>,
#> #   PCR_primers <chr>, note <chr>
```

Take a look at the [introductory
vignette](https://matteoce.github.io/barcodeMineR/articles/intro_barcodeMineR.html)
for a more complete tutorial!

After that, proceed with the [‘Speeding up’
vignette](https://matteoce.github.io/barcodeMineR/articles/api_rate_barcodeMineR.html)
to learn how to increase speed and reliability of the package functions.
