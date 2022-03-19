
<!-- README.md is generated from README.Rmd. Please edit that file -->

# elastes

<!-- badges: start -->

[![R-CMD-check](https://github.com/mpff/elastes/workflows/R-CMD-check/badge.svg)](https://github.com/mpff/elastes/actions)
[![codecov](https://codecov.io/gh/mpff/elastes/branch/main/graph/badge.svg?token=p0xOHfDpnk)](https://codecov.io/gh/mpff/elastes)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

elastes is a R package to estimate elastic shape means from sparse and
irregular planar curves.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mpff/elastes")
```

## Example

Calculate a smooth elastic shape mean for sparse spirals.

``` r
library(elastes)
#> Loading required package: elasdics

# define spiral curve
curve <- function(t){
  rbind(t*cos(13*t), t*sin(13*t))
}

# randomly draw sparse spirals with noise
set.seed(18)
data_curves <- lapply(1:10, function(i){
  m <- sample(10:15, 1)
  delta <- abs(rnorm(m, mean = 1, sd = 0.05))
  t <- cumsum(delta)/sum(delta)
  data.frame(t(curve(t)) + 0.07*t*matrix(cumsum(rnorm(2*length(delta))),
             ncol = 2))
})

# randomly rotate and scale curves
rand_scale <- function(curve){ ( 0.5 + runif(1) ) * curve }
rand_rotate <- function(curve){
  names <- colnames(curve)
  theta <- 2*pi*runif(1)
  mat <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2, ncol = 2)
  curve.rot <- as.matrix(curve) %*% t(mat)
  curve.rot <- as.data.frame(curve.rot)
  colnames(curve.rot) <- names
  return(curve.rot)
}
data_curves <- lapply(data_curves, rand_scale)
data_curves <- lapply(data_curves, rand_rotate)

# compute smooth procrustes mean with 2nd order penalty
knots <- seq(0,1, length = 11)
elastic_shape_mean <- compute_elastic_shape_mean(
  data_curves,
  knots = knots,
  type = "smooth",
  penalty = 2
)
plot(elastic_shape_mean)
```

![](man/figures/README-example-1.png)<!-- -->
