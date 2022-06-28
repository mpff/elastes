#' elastes: Elastic Full Procrustes Means for Sparse and Irregular Planar Curves
#'
#' Provides functions for the computation of functional elastic shape
#' means over sets of open planar curves. The package is particularly suitable for
#' settings where these curves are only sparsely and irregularly observed. It uses
#' a novel approach for elastic shape mean estimation, where planar curves are
#' treated as complex functions and a full Procrustes mean is estimated from the
#' corresponding smoothed hermitian covariance surface, which is combined with the
#' methods for elastic mean estimation proposed in Steyer, Stöcker, Greven	(2022).
#' See Stöcker et. al. (2022) for details on the method.
#'
#' Compute a mean for a set of observed curves: \code{\link{compute_elastic_shape_mean}}
#'
#' @docType package
#' @name elastes
NULL
