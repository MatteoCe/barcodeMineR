## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(barcodeMineR)

# search taxonomic information for a species on the NCBI
t_tax <- get_ncbi_taxonomy("Thouarella variabilis")
t_tax

## -----------------------------------------------------------------------------
t_rec <- download_ncbi(t_tax, ask = FALSE)


## ----echo = FALSE, eval=TRUE--------------------------------------------------
t_rec

## ----echo = T, eval=FALSE-----------------------------------------------------
#  t_rec <- download_ncbi(t_tax)

