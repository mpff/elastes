# See https://github.com/cran/elasdics/blob/master/tests/testthat/test_get_evals.R

test_that("get evals default",{
  curve <- function(t){
    c(sin(2*pi*t), cos(2*pi*t))
  }
  expect_equal(as.numeric(get_evals(curve, seq(0,1,0.1))[11,]), c(0,1))
  expect_equal(get_evals(curve, seq(0,1,0.01)), get_evals(curve))
})


test_that("get evals data.frame",{
  data_frame <- data.frame(s = c(0, 0.3, 0.6, 1),
                           x1 = c(1, 0.5, -1, -1),
                           x2 = c(1, -0.5, -1, 1))
  expect_error(get_evals(data_frame), "Parametrisation t must be in the first column!")
  names(data_frame)[1] <- "t"
  expect_equal(nrow(get_evals(data_frame)), 101)
})


data_curves <- list(data.frame(x1 = c(1, 0.5, -1, -1), x2 = c(1, -0.5, -1, 1)),
                    data.frame(x1 = c(0.1,0.7)*sin(1:6), x2 = cos(1:6)))
knots1 <- seq(0, 1, length=5)
knots2 <- seq(0, 1, length=8)
esm1 <- compute_elastic_shape_mean(data_curves, knots = knots1, max_iter = 2)
esm2 <- compute_elastic_shape_mean(data_curves, knots = knots2, type = "polygon")


test_that("Test can get evals for elastic shape mean",{
  expect_equal(colnames(get_evals(esm1)), c("x1","x2"))
  expect_equal(nrow(get_evals(esm2)), length(knots2))
})


test_that("Test can get evals for elastic shape mean with given t_grid",{
  expect_equal(nrow(get_evals(esm1, t_grid = 0:20/20)), 21)
  expect_equal(nrow(get_evals(esm2, t_grid = 0:10/10)), 11)
})
