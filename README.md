
<!-- README.md is generated from README.Rmd. Please edit that file -->

# horf

<!-- badges: start -->

[![R-CMD-check](https://github.com/robingenuer/horf/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/robingenuer/horf/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The `horf` package contains necessary functions to reproduce the
simulation experiments of the article Arlot and Genuer (XXXX). In this
article the rate of convergence of Hold-Out Random Forests bias is
studied. Thus, the package gives an implementation of this Random
Forests variant and a plot function to visualize the bias convergence
rate of trees v.s. forests.

## Installation

You can install the development version of horf from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("robingenuer/horf")
```

## Example

``` r
library(horf)
## 
```

A basic example code already included in the help pages of the package
(quite quick to run, but not so interesting because the different
parameters are set to quite small values):

``` r
rawForestBias <- comput_bias(nbleaves = c(2^5, 2^6), freg.name = "sinus",
   xdim = 1, nbobs = 640, nbobs_test = 60, nfor = 10, var_estim = TRUE,
   mc.cores = 2)
#> [1] "bias forest k = 32 xdim = 1 mtry =  1 done"
#> Time difference of 2.44238 secs
#> [1] "bias forest k = 64 xdim = 1 mtry =  1 done"
#> Time difference of 5.287689 secs
#> [1] "Total time = 7.73211169242859 secs"
rawForestConfint <- comput_confint(rawForestBias, nbBootstrap = 250)
rawTreeBias <- comput_bias(nbleaves = c(2^5, 2^6), freg.name = "sinus",
   xdim = 1, nbobs = 640, nbobs_test = 60, nfor = 10, tree = TRUE,
   var_estim = TRUE, mc.cores = 2)
#> [1] "bias tree k = 32 xdim = 1 mtry =  1 done"
#> Time difference of 0.1073971 secs
#> [1] "bias tree k = 64 xdim = 1 mtry =  1 done"
#> Time difference of 0.1904933 secs
#> [1] "Total time = 0.299024820327759 secs"
rawTreeConfint <- comput_confint(rawTreeBias, nbBootstrap = 250)
plot_bias(rawForestBias, rawForestConfint, rawTreeBias, rawTreeConfint)
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

Another example (that takes several minutes to run on a 3 cores laptop)
that allows to get a graph closer to the ones included in the paper
would be:

``` r
forestBias <- comput_bias(nbleaves = 2^(5:8), freg.name = "sum",
   xdim = 2, nbobs = 1280, nbobs_test = 100, nfor = 25, 
   var_estim = TRUE, mc.cores = 3)
#> [1] "bias forest k = 32 xdim = 2 mtry =  1 done"
#> Time difference of 6.181252 secs
#> [1] "bias forest k = 64 xdim = 2 mtry =  1 done"
#> Time difference of 13.36967 secs
#> [1] "bias forest k = 128 xdim = 2 mtry =  1 done"
#> Time difference of 39.89349 secs
#> [1] "bias forest k = 256 xdim = 2 mtry =  1 done"
#> Time difference of 3.656062 mins
#> [1] "Total time = 4.64681521256765 mins"
forestConfint <- comput_confint(forestBias, nbBootstrap = 500)
treeBias <- comput_bias(nbleaves = 2^(5:8), freg.name = "sum",
   xdim = 2, nbobs = 1280, nbobs_test = 100, nfor = 100, 
   tree = TRUE, var_estim = TRUE, mc.cores = 3)
#> [1] "bias tree k = 32 xdim = 2 mtry =  1 done"
#> Time difference of 5.561715 secs
#> [1] "bias tree k = 64 xdim = 2 mtry =  1 done"
#> Time difference of 5.772569 secs
#> [1] "bias tree k = 128 xdim = 2 mtry =  1 done"
#> Time difference of 6.633436 secs
#> [1] "bias tree k = 256 xdim = 2 mtry =  1 done"
#> Time difference of 8.674646 secs
#> [1] "Total time = 26.6431660652161 secs"
treeConfint <- comput_confint(treeBias, nbBootstrap = 500)
plot_bias(forestBias, forestConfint, treeBias, treeConfint)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />
