
<!-- README.md is generated from README.Rmd. Please edit that file -->

# aof package

[![Build
Status](https://travis-ci.org/frareb/osfi.svg?branch=master)](https://travis-ci.org/frareb/osfi)
[![CRAN
version](https://www.r-pkg.org/badges/version/aof)](https://CRAN.R-project.org/package=aof)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/grand-total/aof)](https://CRAN.R-project.org/package=aof)
[![Coverage
Status](https://img.shields.io/codecov/c/gh/frareb/osfi/master.svg)](https://codecov.io/gh/frareb/osfi?branch=master)

## Purpose of the package

A breakpoint-based method to detect ontogenetic shifts in univariate
time-activity budget series of central-place foraging insects. The
method finds a single breakpoint according to the likelihood function.
The method was developed with honey bees in order to detect the Age at
Onset of Foraging (AOF), but can be used for the detection of other
ontogenetic shifts in other central-place foraging insects.

## Installation instructions

``` r
# from CRAN
install.packages("aof")

# Or the development version from GitHub:
# install.packages("devtools")
devtools::install_github("frareb/aof")
```

## Authors’ affiliation

UMR EGCE (IRD, CNRS, Univ. Paris-Saclay), Gif-sur-Yvette, France

To cite this package please use the article describing the methods:

Requier, F., Henry, M., Decourtye, A., Brun, F., Aupinel, P., Rebaudo,
F., Bretagnolle, V. (2020) Measuring ontogenetic shifts in central-place
foragers: a case study with honey bees. Journal of Animal Ecology
<https://doi.org/10.1111/1365-2656.13248>

and/or the package itself:

``` r
citation("aof") 
#> 
#> Requier F, Rebaudo F (2020). _aof: Ontogenetic Shifts in Central-Place
#> Foraging Insects_. R package version 0.1.3, <URL:
#> https://cran.r-project.org/package=aof>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {aof: Ontogenetic Shifts in Central-Place Foraging Insects},
#>     author = {Fabrice Requier and François Rebaudo},
#>     year = {2020},
#>     note = {R package version 0.1.3},
#>     url = {https://cran.r-project.org/package=aof},
#>   }
```
