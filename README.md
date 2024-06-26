
<!-- README.md is generated from README.Rmd. Please edit that file -->

**NOTE: This is a temporary GitHub repository created for testing
purposes.**

# barcodeMineR <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/MatteoCe/barcodeMineR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/MatteoCe/barcodeMineR/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/MatteoCe/barcodeMineR/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/MatteoCe/barcodeMineR/actions/workflows/test-coverage.yaml)
[![codecov](https://codecov.io/gh/MatteoCe/barcodeMineR/graph/badge.svg?token=62OUVOL8MP)](https://codecov.io/gh/MatteoCe/barcodeMineR)
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
fail. It is not necessarily an issue, thus you may proceed nonetheless.

## Getting started

Take a look at the [barcodeMineR
website](https://matteoce.github.io/barcodeMineR/) for more information!

## Issues and ‘to dos’

- Remove field *lengthGene*? It might have no purpose, considering the
  sequence itself is in the refdb object and can be modified using the
  refdb package functions.

- Change default *api_rate* in ncbi functions from 3 to 2.8, to slow
  down requests that might build up to 4 per seconds with slow internet
  connections.

- Change default *ask* argument from TRUE to FALSE in both bold and ncbi
  functions?

- Add stringr code in taxonomy functions to search for both genus and
  species names in the output in order to get also those records
  identified as *genus* cf. *species*, for example.
