#' Align two curves measured at discrete points
#' @name align_curves_proc2d
#' @description Finds the optimal reparametrisation and procrustes fit of the second
#' curve (stored in \code{pfit_curve2}) onto the normalized unit-length
#' representaton of the first one (stored in \code{pfit_curve1}) with respect
#' to the elastic full Procrustes distance. Constructor function for class
#' \code{aligned_curves_proc2d}.
#' @param data_curve1 \code{data.frame} with observed points in each row. Each
#' variable is one coordinate direction. If there is a variable \code{t},
#' it is treated as the time parametrisation, not as an additional coordinate.
#' @param data_curve2 same as \code{data_curve1}
#' @param eps convergence tolerance
#' @return an object of class \code{aligned_curves_proc2d}, which is a \code{list}
#' with entries
#'   \item{pfit_curve1}{\code{pfit_curve1} with parametrisation variable \code{t}}
#'   \item{pfit_curve2_aligned}{\code{pfit_curve2} with initial parametrisation
#'   variable \code{t} and optimal parametrisation \code{t_optim}}
#'   \item{elastic_full_procrustes_dist}{elasic full procrustes distance between curve1 and curve2}
#' @export
#' @exportClass aligned_curves_proc2d
#' @examples
#' data_curve1 <- data.frame(x1 = c(1, 0.5, -1, -1), x2 = c(1, -0.5, -1, 1))
#' data_curve2 <- data.frame(x1 = c(0.1,0.7)*sin(1:6), x2 = cos(1:6))
#' aligned_curves <- align_curves(data_curve1, data_curve2)
#'
#' #different parametrisation of the first curve
#' data_curve1$t <- 0:3/3
#' align_curves(data_curve1, data_curve2)


align_curves_proc2d <- function(data_curve1, data_curve2, eps = 0.01){
  closed = FALSE  # for now only open curves
  #remove duplicated points
  data_curve1 <- remove_duplicate(data_curve1)
  data_curve2 <- remove_duplicate(data_curve2)
  if(attributes(data_curve1)$points_rm | attributes(data_curve2)$points_rm){
    warning("Duplicated points in data curves have been removed!")
  }
  # input checking given parametrisation t
  if("t" %in% names(data_curve1)) check_param(data_curve1, closed)
  if("t" %in% names(data_curve2)) check_param(data_curve2, closed)

  # input checking for closed curves
  if(closed){
    data_curve1 <- check_closed(data_curve1)
    data_curve2 <- check_closed(data_curve2)
  }

  # Prepare pfit curves
  length1 <- get_polygon_length(data_curve1)
  length1 <- get_polygon_length(data_curve1)
  pfit_curve_1 <- elasdics::center_curve(data_curve1)/length1
  pfit_curve_2 <- elasdics::center_curve(data_curve2)/length2
  # Prepare normalized SRV curves
  srv_data_1 <- get_srv_from_points(pfit_curve1)
  srv_data_2 <- get_srv_from_points(pfit_curve2)

  #remove parametrisation
  if("t" %in% names(pfit_curve1)) pfit_curve1$t  <- NULL
  if("t" %in% names(pfit_curve2)) pfit_curve2$t  <- NULL
  # after computing the srv transformations the parametrisation t is definitly
  # in the first column
  if(closed){
    #pre alignment
    #t <- get_arc_length_param(data_curve2)
    #loss <- sapply(1:length(srv_data_2$t), function(i){
    #  t_new <- t - t[i]
    #  compute_distance(srv_data_1, srv_data_2, t_new, closed)
    #})
    #initial_t <- t - t[which.min(loss)]
    #optimal_t <- find_optimal_t_discrete_closed(c(srv_data_1$t,1), t(srv_data_1[,-1]),
    #                                     c(srv_data_2$t,1), t(srv_data_2[,-1]),
    #                                     initial_t = initial_t, eps = eps)
    #t_optim <- optimal_t + (optimal_t < 0)
    #t_optim <- ifelse(t_optim > 1, t_optim - 1, t_optim)
  } else {
    initial_t <- get_arc_length_param(pfit_curve2)
    pfit_srv_data_2 <- srv_data_2

    ##adjust convergence criterium for number of points
    #eps <- eps/length(t)
    #delta <- eps
    #while(delta >= eps){

    #stop_crit <- sum((coefs - coefs_old)^2)/sum(coefs^2)
    #if(stop_crit < eps |  == 0){
    #    }

    t_optim <- find_optimal_t_discrete(c(srv_data_1$t,1), t(srv_data_1[,-1]),
                                       c(srv_data_2$t,1), t(srv_data_2[,-1]),
                                       initial_t = initial_t, eps = eps)
  }

  elastic_fp_dist <- compute_distance(srv_data_1, pfit_srv_data_2,
                                                   t_optim, closed)

  pfit_curve1 <- cbind(t = c(srv_data_1$t, 1), pfit_curve1)
  pfit_curve2_aligned <- cbind(t = c(srv_data_2$t, 1),
                               t_optim = t_optim, pfit_curve2)

  aligned_curves_proc2d <- list("pfit_curve1" = pfit_curve1,
                         "pfit_curve2_aligned" = pfit_curve2_aligned,
                         "elastic_fp_dist" = elastic_fp_dist,
                         "closed" = closed)
  class(aligned_curves_proc2d) <- "aligned_curves_proc2d"
  return(aligned_curves_proc2d)

}

#' Input checking for given parametrisation
#' @inheritParams align_curves
#' @param data_curve data of curve like in \code{align_curves}
#' @noRd

check_param <- function(data_curve, closed = closed){
  if(!(data_curve$t[1] >= 0 & data_curve$t[nrow(data_curve)] <= 1 &
       all(diff(data_curve$t) >= 0))){
    stop("Parametrisation t needs to be within 0 and 1 and increasing!")
    }
  if(data_curve$t[1] != 0){
    stop("Parametrisation t needs to start at 0!")
    }
  if(!closed & data_curve$t[nrow(data_curve)] != 1){
    stop("Last value of parametrisation t needs to be 1!")
  }
}



#' Input checking for closed curves
#' @param data_curve data of curve like in \code{align_curves}
#' @noRd

check_closed <- function(data_curve){
  if("t" %in% names(data_curve)){
    if(data_curve$t[nrow(data_curve)] == 1){
      if(!all(data_curve[1, names(data_curve) != "t"] ==
         data_curve[nrow(data_curve), names(data_curve) != "t"])){
        stop("Curve is not closed")
      }
    } else {
      data_curve <- rbind(data_curve, data_curve[1,])
      data_curve$t[nrow(data_curve)] <- 1
    }
  } else {
    err_non_closed <- sum((data_curve[1,] - data_curve[nrow(data_curve),])^2)/sum(data_curve[1,]^2)
    if(err_non_closed > sqrt(.Machine$double.eps)) {
      data_curve <- rbind(data_curve, data_curve[1,])
    } else {
      data_curve[nrow(data_curve),] <- data_curve[1,]
    }
  }
  return(data_curve)
}



#' Input checking for planar curves
#' @param data_curve data of curve like in \code{align_curves}
#' @noRd

check_planar <- function(data_curve){
  if("t" %in% names(data_curve)){
    if(ncol(data_curve) - 1 != 2){
      stop("Curve is not planar. Option 'proc2d' is only valid for planar curves.")
    }
  } else {
    if(ncol(data_curve) != 2){
      stop("Curve is not planar. Option 'proc2d' is only valid for planar curves.")
    }
  }
  return(data_curve)
}



#' Remove duplicated points in data curves
#' @param data_curve data of curve like in \code{align_curves}
#' @noRd

remove_duplicate <- function(data_curve){
  points <- as.data.frame(data_curve)
  try(points$t <- NULL, silent = TRUE)
  moves <- c(TRUE, rowSums(apply(points, 2, diff)^2) != 0)
  data_curve <- data_curve[moves,]
  attr(data_curve, "points_rm") <- !all(moves)
  data_curve
}

#' Computes elastic distance
#' @noRd
#' @param srv_data_1 srv of curve1
#' @param srv_data_2 srv of curve2
#' @param t_optim optimal parametrisation of curve2

compute_distance<- function(srv_data_1, srv_data_2, t_optim, closed){
  norm_1 <- sum(t(srv_data_1[,-1]^2)%*%diff(c(srv_data_1$t,1)))
  norm_2 <- sum(t(srv_data_2[,-1]^2)%*%diff(c(srv_data_2$t,1)))
  if(closed){
    srv_data_1_extended <- rbind(srv_data_1, srv_data_1, srv_data_1)
    srv_data_1_extended$t <- c(srv_data_1$t - 1, srv_data_1$t, srv_data_1$t + 1)
    cross_prod <- get_loss_discrete(t = t_optim, srv_data_1 = srv_data_1_extended, srv_data_2)
  } else {
    cross_prod <- get_loss_discrete(t = t_optim, srv_data_1, srv_data_2)
  }
  sqrt(norm_1 + norm_2 - 2*cross_prod)
}


