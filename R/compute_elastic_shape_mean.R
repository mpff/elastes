#' Compute an elastic full Procrustes mean for a collection of curves
#' @name compute_elastic_shape_mean
#' @description Computes an elastic full Procrustes mean for curves stored in \code{data_curves}.
#' Constructor function for class \code{elastic_shape_mean}.
#' @param data_curves list of \code{data.frame}s with observed points in each row. Each
#' variable is one coordinate direction. If there is a variable \code{t},
#' it is treated as the time parametrisation, not as an additional coordinate.
#' @param knots set of knots for the mean spline curve
#' @param type if "smooth" linear srv-splines are used which results in a differentiable mean curve
#' if "polygon" the mean will be piecewise linear.
#' @param penalty the penalty to use in the covariance smoothing step. use '-1' for no penalty.
#' @param pfit_method (experimental) "smooth" or "polygon"
#' @param var_type (experimental) assume "smooth", "constant" or "zero" measurement-error variance along t
#' @param smooth_warp (experimental) controls the weighting of original and smoothed observations
#' over the iterations, if pfit_method == "smooth".
#' @param eps the algorithm stops if L2 norm of coefficients changes by less than \code{eps}
#' @param max_iter maximal number of iterations
#' @param verbose print iterations
#' @param cluster (experimental) use the parallel package for faster computation
#' @return an object of class \code{elastic_shape_mean}, which is a \code{list}
#' with entries
#'   \item{type}{"smooth" if mean was modeled using linear srv-splines, "polygon" if constant srv-splines}
#'   \item{coefs}{spline coeffiecients}
#'   \item{knots}{spline knots}
#'   \item{variance}{sample elastic shape variance}
#'   \item{data_curves}{list of \code{data.frame}s with observed points in each row.
#'   First variable \code{t} gives the initial parametrisation, second variable \code{t_optim}
#'   the optimal parametrisation when the curve is aligned to the mean. Has the
#'   attributes 'rotation', 'scaling', 'translation' and 'dist_to_mean'. Use
#'   \code{\link{get_procrustes_fit}} to get the elastic full Procrustes fit.}
#'  \item{fit}{see \code{fit_mean}}
#' @export
#' @examples
#' curve <- function(t){
#'   rbind(t*cos(13*t), t*sin(13*t))
#' }
#' set.seed(18)
#' data_curves <- lapply(1:4, function(i){
#'   m <- sample(10:15, 1)
#'   delta <- abs(rnorm(m, mean = 1, sd = 0.05))
#'   t <- cumsum(delta)/sum(delta)
#'   data.frame(t(curve(t)) + 0.07*t*matrix(cumsum(rnorm(2*length(delta))),
#'              ncol = 2))
#' })
#'
#' #randomly rotate and scale curves
#' rand_scale <- function(curve){ ( 0.5 + runif(1) ) * curve }
#' rand_rotate <- function(curve){
#'   names <- colnames(curve)
#'   theta <- 2*pi*runif(1)
#'   mat <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2, ncol = 2)
#'   curve.rot <- as.matrix(curve) %*% t(mat)
#'   curve.rot <- as.data.frame(curve.rot)
#'   colnames(curve.rot) <- names
#'   return(curve.rot)
#' }
#' data_curves <- lapply(data_curves, rand_scale)
#' data_curves <- lapply(data_curves, rand_rotate)
#'
#' #compute smooth procrustes mean with 2 order penalty
#' knots <- seq(0,1, length = 11)
#' elastic_shape_mean <- compute_elastic_shape_mean(
#'     data_curves,
#'     knots = knots,
#'     type = "smooth",
#'     penalty = 2
#'     )
#' plot(elastic_shape_mean)


compute_elastic_shape_mean <- function(data_curves, knots = seq(0, 1, len = 13),
                                  type = c("smooth", "polygon"), penalty = 2, var_type = c("smooth", "constant", "zero"),
                                  pfit_method = c("smooth", "polygon"), smooth_warp = function(i) 0,
                                  eps = 0.05, max_iter = 50, verbose = F, cluster = NULL) {

  if(verbose) message("Initializing...")
  # Input checks
  stopifnot(all(sapply(data_curves, is.data.frame)))

  # Match arguments
  type <- match.arg(type, c("smooth", "polygon"))
  var_type <- match.arg(var_type, c("smooth", "constant", "zero"))
  pfit_method <- match.arg(pfit_method, c("smooth", "polygon"))

  # remove duplicated points
  data_curves <- lapply(data_curves, remove_duplicate, closed = F)
  if(sum(sapply(data_curves, function(curve){attributes(curve)$points_rm}) > 0)){
    warning("Duplicated points in data curves have been removed!")
  }
  data_curves <- lapply(data_curves, function(curve){
    attr(curve, "points_rm") <- NULL
    curve
  })

  # parametrisation with respect to arc length if not given,
  # after this, parametrisation is always in the first column
  data_curves <- lapply(data_curves, function(data_curve){
    if(!("t" %in% colnames(data_curve))){
      data.frame("t" = elasdics::get_arc_length_param(data_curve), data_curve)
    } else {
      param <- data_curve$t
      data_curve$t <- NULL
      data.frame("t" = param, data_curve)
    }
  })

  # Adjust scaling and centering of data_curves
  translations <- lapply(data_curves, get_center)
  lengths <- lapply(data_curves, get_polygon_length)
  data_curves_adj <- lapply(data_curves, elasdics::center_curve)
  # "Prescale" data curves to polygon unit-length.
  data_curves_adj <- lapply(1:length(data_curves_adj), function(i){
    data.frame("t" = data_curves_adj[[i]]$t, data_curves_adj[[i]][,-1]/lengths[[i]])
  })

  if(length(unique(sapply(data_curves, ncol))) != 1) stop("All curves must have same number of dimensions!")

  # Calculate normalized SRV data curves. (Note: Normalization comes from unit-length of the data_curves_adj!)
  srv_data_curves <- lapply(data_curves_adj, elasdics::get_srv_from_points)

  if(ncol(srv_data_curves[[1]]) != 3) stop("This package was designed to analyse only planar curve data!")

  # Calculate elastic full Procrustes mean
  elastic_shape_mean <- fit_mean(srv_data_curves = srv_data_curves,
    knots = knots, type = type, penalty = penalty, var_type = var_type, pfit_method = pfit_method,
    max_iter = max_iter, eps = eps, verbose = verbose, cluster = cluster, smooth_warp = smooth_warp)

  # Add scaling, rotation, translation and distance attributes to the original data curves.
  data_curves <- lapply(1:length(data_curves), function(j){
    data_curves[[j]]$t_optim <- elastic_shape_mean$t_optims[[j]]
    attributes(data_curves[[j]]$t_optim) <- NULL
    data_curve <- data_curves[[j]][, c(1, 4, 2, 3)]
    attr(data_curve, "dist_to_mean") <- elastic_shape_mean$distances[[j]]
    attr(data_curve, "rotation") <- elastic_shape_mean$fit$G_optims[[j]]
    attr(data_curve, "scale") <- 1/elastic_shape_mean$fit$b_optims[[j]]^2
    attr(data_curve, "polygon_length") <- lengths[[j]]
    attr(data_curve, "norm_factor") <- elastic_shape_mean$fit$l_optims[[j]]
    attr(data_curve, "translation") <- translations[[j]]
    data_curve
  })

  # Calculate variance.
  eigenvals <- elastic_shape_mean$fit$cov_pca$values
  eigenvals <- eigenvals[eigenvals > 0]
  variance <- 1 - eigenvals[1]/sum(eigenvals)

  # Return
  elastic_shape_mean$data_curves <- data_curves
  elastic_shape_mean$shift_idxs <- NULL
  elastic_shape_mean$t_optims <- NULL
  elastic_shape_mean$variance <- variance
  class(elastic_shape_mean) <- "elastic_shape_mean"

  if(verbose) message("Finished.")
  elastic_shape_mean
}
