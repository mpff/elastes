test_that("Test variance negative/inf in issue #8 example.", {
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+~~~
  # Manuel's negative variance example   15.02.2022
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+~~~~~~~~
  # https://github.com/mpff/elasdicsproc2d/issues/8

  # define curve
  curve <- function(t) {
    rbind(t * cos(2 * t), t * sin(2 * t))
  }

  # randomly draw sparse curve with noise
  data_curves <- lapply(1:4, function(i) {
    m <- sample(10:15, 1)
    delta <- abs(rnorm(m, mean = 1, sd = 0.05))
    t <- cumsum(delta) / sum(delta)
    data.frame(t(curve(t)) + 0.07 * t * matrix(cumsum(rnorm(2 * length(delta))),
                                               ncol = 2
    ))
  })

  expect_error(compute_elastic_shape_mean(data_curves, var_type = "zero"), NA)
  expect_error(compute_elastic_shape_mean(data_curves, type = "polygon"), NA)
})


test_that("Test no dropped points in digit3 polygon example.",{
  data_curves <- list(
    data.frame(
      X1 = c(-14.30, -11.30, -6.30, 2.69, 10.69, 12.69, 14.69, 11.69, 6.69, -2.30, -2.30, -7.30, -15.30),
      X2 = c(-4.15, -8.15, -13.15, -16.15, -14.15, -10.15, -4.15, 3.84, 7.84, 8.84, 14.84, 16.84, 17.84)
    ),
    data.frame(
      X1 = c(-5.61, -1.61, 3.38, 4.38, 2.38, -0.61, -3.61, 1.38, 3.38, 5.38, 3.38, -4.61, -7.61),
      X2 = c(-14.23, -12.23, -10.23, -6.23, -2.23, -1.23, -3.23, 0.76, 5.76, 9.76, 12.76, 11.76, 8.76)
    ),
    data.frame(
      X1 = c(-6.92, -1.92, 3.07, 4.07, 1.07, -4.92, -8.92, 1.07, 4.07, 5.07, 5.07, 1.07, -1.92),
      X2 = c(-13.07, -13.07, -8.07, -4.07, 0.92, -0.07, -1.07, 0.92, 2.92, 5.92, 8.92, 9.92, 9.92)
    ),
    data.frame(
      X1 = c(-9.30, -3.30, 5.69, 10.69, 5.69, 1.69, -6.30, 1.69, 5.69, 2.69, -0.30, -5.30, -9.30),
      X2 = c(-14.46, -17.46, -15.46, -10.46, -4.46, -0.46, 3.53, 3.53, 5.53, 9.53, 11.53, 13.53, 15.53)
    )
  )
  expect_error(compute_elastic_shape_mean(data_curves, type = "polygon"), NA)

})
