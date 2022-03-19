#' Plot method for planar elastic procrustes mean curves
#' @description Plots objects of class \code{elastic_shape_mean}.
#' @param x object of class \code{elastic_shaped_mean},
#' usually a result of a call to \code{\link{compute_elastic_shape_mean}}
#' @param asp numeric, giving the aspect ratio of the two coordinates,
#' see \code{\link{plot.window}} for details.
#' @param col color of the mean curve.
#' @param srv TRUE if the SRV curve shall be plotted
#' @param centering TRUE if mean and pfits shalle be centered
#' @param ... further plotting parameters.
#' @importFrom graphics plot lines
#' @export
#'
#' @seealso For examples see documentation of \code{\link{compute_elastic_shape_mean}}.

plot.elastic_shape_mean <- function(x, srv = FALSE, centering = TRUE, asp = 1, col = "red", ...){
  if(ncol(x$coefs) != 2){
    stop("Plotting option only for planar curves")
  }
  data_curves <- lapply(x$data_curves, function(data_curve){
    data_curve <- get_procrustes_fit(data_curve)
    if(srv == TRUE) {
      if("t_optim" %in% colnames(data_curve)) {
        data_curve$t <- data_curve$t_optim
        data_curve <- data_curve[, names(data_curve) != "t_optim"]
      }
      data_curve <- elasdics::get_srv_from_points(data_curve) }
    else {
      if(centering == FALSE){
        data_curve[,3] <- data_curve[,3] - data_curve[1,3]
        data_curve[,4] <- data_curve[,4] - data_curve[1,4]
      }
    }
    data_curve
    })
  data_curves <- lapply(data_curves, function(data) data[,colnames(x$coefs)])
  data_all <- do.call("rbind", data_curves)

  #empty plot
  plot(NULL, xlim = range(data_all[,1]), ylim = range(data_all[,2]), xlab = colnames(x$coefs)[1],
       ylab = colnames(x$coefs)[2], asp = 1, ...)

  #plot data
  lapply(data_curves, lines, col = "gray")

  #plot mean
  lines(get_evals(x, srv = srv, centering = centering), col = col, lwd = 2)
}
