# Define data curves
curve <- function(t){
  rbind(t*cos(13*t), t*sin(13*t))
}
set.seed(18)
data_curves <- lapply(1:4, function(i){
  m <- sample(10:15, 1)
  delta <- abs(rnorm(m, mean = 1, sd = 0.05))
  t <- cumsum(delta)/sum(delta)
  data.frame(t(curve(t)) + 0.07*t*matrix(cumsum(rnorm(2*length(delta))),
                                         ncol = 2))
})
#compute elastic means
knots <- seq(0,1, length = 11)
esm <- compute_elastic_shape_mean(data_curves, knots = knots)

test_that("Test plot runs without warning",{
  expect_warning(plot(esm), regexp = NA)
})

test_that("Test can plot srv curve",{
  expect_warning(plot(esm, srv = TRUE), regexp = NA)
})

test_that("Test can plot centered curve",{
  expect_warning(plot(esm, centering = TRUE), regexp = NA)
})
