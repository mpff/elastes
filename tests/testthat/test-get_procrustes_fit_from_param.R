test_that("proc2d from with param and proc2d fit from inverse param are inverse",{
  data_curve <- data.frame(x1 = 1:6*sin(1:6), x2 = cos(1:6))
  data_curve_new <- get_procrustes_fit_from_param(get_procrustes_fit_from_param(data_curve, pi, 2, 1, 1, 2), -pi, 0.5, 1, -1, 0.5)
  expect_equal(apply(data_curve_new, 2, diff),
               apply(data_curve, 2, diff))
  data_curve$t <- 0:5/5
  data_curve_new2 <- get_procrustes_fit_from_param(get_procrustes_fit_from_param(data_curve, pi, 2, 1, 1, 2), -pi, 0.5, 1, -1, 0.5)
  expect_equal(apply(data_curve_new2, 2, diff),
               apply(data_curve, 2, diff))
})
