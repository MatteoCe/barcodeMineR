
<!-- README.md is generated from README.Rmd. Please edit that file -->

**NOTE: This is a temporary GitHub repository created for testing
purposes.**

# barcodeMineR <img src='man/logo/barcodeMineR_logo.png' align="right" height="220" />

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
