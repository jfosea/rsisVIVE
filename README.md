
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rsisVIVE

![](man/figures/logo.png)

The goal of rsisVIVE (robustified some invalid some valid instrumental
variable estimator) is to estimate the causal relationship between an
outcome and an exposure in the presence of invalid instruments while
tolerating large proportions of contamination in the exposure and/or
outcome values. The algorithm follows closely with the sisVIVE method of
Kang et. al. (2016). The rsisVIVE replaces the two step algorithm of the
sisVIVE with robust counterparts.

More details of the method can be seen in the master’s thesis submission
by Jana Osea at the University of British Columbia as a partial
fulfillment of the requirements for Master’s of Science, Department of
Statistics under the supervision of [Dr. Gabriela
Cohen-Freue](https://gcohenfr.github.io/).

## Installation

You can install the development version of rsisVIVE like so:

``` r
devtools::install_github("jfosea/rsisVIVE")
library(rsisVIVE)
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(rsisVIVE)
#> Warning: replacing previous import 'MASS::select' by 'dplyr::select' when
#> loading 'rsisVIVE'

n <- 100
L <- 10
s <- 3
m <- floor(0.2*n)
Sigma <- diag(1, nrow = L)

a0 <- 1
a1 <- rep(0.05, L)

b0 <- 1
b1 <- c(rep(1,s), rep(0, L-s)) # alpha
b2 <- 1                        # beta (parameter of interest)

error <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, 0.9, 0.9, 1), 2, 2))

Z <- MASS::mvrnorm(n, mu = rep(0, L), Sigma)
D <- a0 + Z %*% a1 + error[, 1]
Y <- b0 + Z %*% b1 + D * b2 + error[, 2]

Y[sample(1:n, m)] <- rnorm(m, 70, 1) # add contamination
rsisVIVE(Y, D, Z, method = 'PE_SE', ncores = 1)
#> $alpha
#>          x1          x2          x3          x4          x5          x6 
#>  0.00000000  0.00000000 -0.09379886  0.00000000  0.00000000 -0.25057890 
#>          x7          x8          x9         x10 
#> -1.94291405  0.32676347  0.00000000  1.31882800 
#> 
#> $beta
#> [1] 1.920595
```

We compare with the results of `sisVIVE`.

``` r
# compare with sisVIVE
sisVIVE::cv.sisVIVE(Y, D, Z)
#> $lambda
#> [1] 160.7302
#> 
#> $estCVError
#> [1] 7155.194
#> 
#> $alpha
#>  [1] 0 0 0 0 0 0 0 0 0 0
#> 
#> $beta
#> [1] 8.768325
#> 
#> $whichInvalid
#> [1] ""
```

## References

Kang, H., Zhang, A., Cai, T. T., and Small, D. S. (2016). [Instrumental
Variables EStiamtion with Some Invalid Instruments and its Application
to Mendelian
Randomization](https://www.tandfonline.com/doi/full/10.1080/01621459.2014.994705).
Journal of the American Statistical Association , 111, 132-144.
