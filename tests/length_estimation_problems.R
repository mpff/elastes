
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+~~~
# Manuel's length estimation example   25.01.2022
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+~~~~~~~~
# https://github.com/mpff/elasdicsproc2d/issues/2#issue-1114129028

library(elasdicsproc2d)

# define spiral curve
curve <- function(t) {
  rbind(t * cos(13 * t), t * sin(13 * t))
}

# randomly draw sparse spirals with noise
set.seed(18)
data_curves <- lapply(1:10, function(i) {
  m <- sample(10:15, 1)
  delta <- abs(rnorm(m, mean = 1, sd = 0.05))
  t <- cumsum(delta) / sum(delta)
  data.frame(t(curve(t)) + 0.07 * t * matrix(cumsum(rnorm(2 * length(delta))),
                                             ncol = 2
  ))
})

# Use smoothing in procrustes fit calculation and re-normalisation.
elastic_proc2d_mean_smooth <- compute_elastic_proc2d_mean(
  data_curves,
  knots = seq(0, 1, length = 13),
  type = "smooth",
  penalty = 2,
  var_type = "smooth",
  pfit_method = "smooth"
)

# Plot the estimated mean shape with procrustes fits
plot(elastic_proc2d_mean_smooth, main = "Mean (pfit_method = 'smooth')")
#> Error in fit_alignment_proc2d(q, type, knots, var_type, coefs.compl, method = pfit_method,  :
#>  In calculation of curve length: Length is negative or imaginary. l = -0.781801207724315+0i
