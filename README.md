
<!-- README.md is generated from README.Rmd. Please edit that file -->

**NOTE: This is a temporary GitHub repository created for testing
purposes.**

# barcodeMineR <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
<!-- badges: end -->

The goal of barcodeMineR is to facilitate the download of DNA sequences
from the NCBI and BOLD repositories, taking advantage of the future
asynchronous framework to speed the operations while respecting the API
requests limits. It outputs a
[refdb](https://github.com/fkeck/refdb?tab=readme-ov-file#refdb-a-reference-database-manager-for-r)
object.

## Installation

As the package is still in development, you can install this version of
barcodeMineR, directly from the GitHub repository, using the
[Bioconductor](https://www.bioconductor.org/install/) package:

``` r
BiocManager::install("MatteoCe/barcodeMineR")
```

The installation might require the update of various packages, which may
fail. It is not necessarily an issue, thus you an proceed nonetheless.

## Example

The first step consists in loading the package:

``` r
library(barcodeMineR)
```

Take a look at the [vignette](doc/intro_barcodeMineR.html) for more
information!

## TO DO

- Check sequence extract method in buildSequences. If a [location
  descriptor](https://www.insdc.org/submitting-standards/feature-table/#3.4.3)
  of the INSDC include a start or end location of unknown position (for
  example, “\<345..500”), should the function extract the sequence using
  the first or end position of the whole sequence instead of the unknown
  position or include the last/first known position only? (the 345
  position as mentioned in the example). As example of a problem
  arising, try downloading all CDS of the records corresponding to
  *Ophthalmolycus amberensis* (specifically the acc num ON417737.1 and
  the CDS CYTB), or the acc num KX362346.1 for the muts CDS of the
  species *Notisis elongata*.

Another issue related to this problem can be seen with the accession
number “JN034583.1” of the species *Dissostichus mawsoni*, where
different portions of “immunoglobulin M heavy chain” are separated by
exons.

- 
