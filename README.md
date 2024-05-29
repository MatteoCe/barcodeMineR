
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

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
