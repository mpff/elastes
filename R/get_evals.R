#' @title Evaluate a curve on a grid
#' @param curve a one parameter function which is to be evaluated on a grid
#' @param t_grid the curve is evaluated at the values in t_grid, first value needs
#' to be 0, last value needs to be 1.
#' If t_grid = NULL, a default regular grid with grid length 0.01 is chosen
#' @param centering TRUE if curves shall be centered
#' @param srv TRUE if SRV curve shall be evaluated
#' @param ... other arguments
#' @return a \code{data.frame} with evaluations of the curve
#' at the values in \code{t_grid} in its rows.
#' @export
#' @examples
#' curve <- function(t){c(t*sin(10*t), t*cos(10*t))}
#' plot(get_evals(curve), type = "b")
#' @seealso See \code{\link[elasdics]{get_evals}} for the original code.

get_evals <- function(curve, t_grid = NULL, ...){
  UseMethod("get_evals")
}

#' @export
get_evals.default <- function(curve, t_grid = NULL,  ...){
  if(is.null(t_grid)) t_grid <- seq(0,1, by = 0.01)
  data.frame(t(sapply(t_grid, curve)))
}


#' @rdname get_evals
#' @export

get_evals.data.frame <- function(curve, t_grid = NULL, ...){
  if("t" != names(curve)[1]) stop("Parametrisation t must be in the first column!")
  if(is.null(t_grid)) t_grid <- seq(0,1, by = 0.01)

  points <- sapply(t_grid, function(t){
    idx <- findInterval(t, curve$t, rightmost.closed = T)
    weights <- c(1,-1)*(curve$t[c(idx +1, idx)] - t)/(curve$t[idx + 1] - curve$t[idx])
    t(weights)%*%as.matrix(curve[c(idx, idx + 1), -1, drop = FALSE])
  })

  if(!is.matrix(points)) points <- t(matrix(points))
  points <- as.data.frame(t(points))
  names(points) <- names(curve)[-1]
  points
}


#' @rdname get_evals
#' @export

get_evals.elastic_shape_mean <- function(curve, t_grid = NULL,
                                   centering = TRUE, srv = FALSE, ...){
  if(curve$type == "smooth"){
    if(is.null(t_grid)) t_grid <- seq(0,1, by = 0.01)
    srv_mean_curve <- function(t){
      t(make_design(t, knots = curve$knots, type = curve$type) %*% curve$coefs)
    }
    if(srv == FALSE){
      mean_data <- as.data.frame(t(srvf_to_curve(t_grid, srv_mean_curve)))
      if(centering) mean_data <- elasdics::center_curve(mean_data)
    } else {
      mean_data <- data.frame(t(srv_mean_curve(t_grid)))
    }
  } else if(curve$type == "polygon"){
    if(!is.null(t_grid)){
      t <- t_grid[-t_grid]
      idx <- findInterval(t, curve$knots, rightmost.closed = T)
      srv_data <- data.frame(t, curve$coefs[idx,])
    } else {
      srv_data <- data.frame("t" = curve$knots[-length(curve$knots)], curve$coefs)
    }
    if(srv == FALSE) {
      mean_data <- elasdics::get_points_from_srv(srv_data)
      if(centering) mean_data <- elasdics::center_curve(mean_data)
    } else {
      mean_data <- srv_data[,-1]
    }
  }
  colnames(mean_data) <- colnames(curve$coefs)
  mean_data
}
