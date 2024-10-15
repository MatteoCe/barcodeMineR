
<!-- README.md is generated from README.Rmd. Please edit that file -->

**NOTE: This package is still in development and awaiting for CRAN
approval.**

# barcodeMineR <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/MatteoCe/barcodeMineR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/MatteoCe/barcodeMineR/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/MatteoCe/barcodeMineR/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/MatteoCe/barcodeMineR/actions/workflows/test-coverage.yaml)
[![codecov](https://codecov.io/gh/MatteoCe/barcodeMineR/graph/badge.svg?token=62OUVOL8MP)](https://app.codecov.io/gh/MatteoCe/barcodeMineR)
<!-- badges: end -->

The goal of barcodeMineR is to facilitate the download of DNA sequences
from the NCBI and BOLD repositories, taking advantage of the future
asynchronous framework of the [future
package](https://github.com/HenrikBengtsson/future) to speed up the
operations while respecting the API requests limits. It outputs a
[refdb](https://github.com/fkeck/refdb?tab=readme-ov-file#refdb-a-reference-database-manager-for-r)
object, including all CDS and/or rRNA feature corresponding to the
searched taxa.

## Installation

The development version
[![](https://img.shields.io/badge/devel%20version-0.1.0-blue.svg)](https://github.com/MatteoCe/barcodeMineR)
of this package can be installed directly from this GitHub repository,
using the [Bioconductor](https://www.bioconductor.org/install/) package:

``` r
BiocManager::install("MatteoCe/barcodeMineR")
```

The last stable version \#put next line up to hashtag in single quotes
when available# r badger::badge_cran_release(NULL, “orange”)# is
available through CRAN:

``` r
install.packages("barcodeMineR")
```

## Getting started

Take a look at the [barcodeMineR
website](https://matteoce.github.io/barcodeMineR/) for more information!

## Citation

If you use this package or its code in your work, please cite:

Cecchetto et al. XXX

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
  identified as *genus* cf. *species*, for example?

## Acknowledgments

This package was written following the instructions and suggestions
described in the fantastic [R Packages book](https://r-pkgs.org) by
[Hadley Wickham](https://hadley.nz/) and [Jennifer
Bryan](https://jennybryan.org/).

The barcodeMineR package was built using the following R packages: \*
[devtools](https://github.com/r-lib/devtools) \*
[usethis](https://github.com/r-lib/usethis) \*
[pkgdown](https://github.com/r-lib/pkgdown) \*
[testthat](https://github.com/r-lib/testthat) \*
[roxygen2](https://github.com/r-lib/roxygen2)
