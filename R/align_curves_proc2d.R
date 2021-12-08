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


align_curves_elastic_proc2d <- function(data_curve1, data_curve2, delta = 0.01, eps = 0.01, h = 0.01){
  warning("This function is not fully implemented yet. At the moment only the alignment inside mean calculation works.")
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

  # Prepare centered, unit-length curves
  length1 <- get_polygon_length(data_curve1)
  length2 <- get_polygon_length(data_curve1)
  pfit_curve1 <- elasdics::center_curve(data_curve1)/length1
  pfit_curve2 <- elasdics::center_curve(data_curve2)/length2
  # Prepare normalized SRV curves
  srv_data_1 <- get_srv_from_points(pfit_curve1)
  srv_data_2 <- get_srv_from_points(pfit_curve2)

  #remove parametrisation
  if("t" %in% names(pfit_curve1)) pfit_curve1$t  <- NULL
  if("t" %in% names(pfit_curve2)) pfit_curve2$t  <- NULL
  # after computing the srv transformations the parametrisation t is definitly
  # in the first column

  initial_t <- get_arc_length_param(pfit_curve2)
  t_optim <- initial_t
  pfit_srv_data_2 <- srv_data_2

  diff <- Inf
  while(diff >= delta){

    pfit_srv_data_2_old <- pfit_srv_data_2

    w_optim <- find_optimal_w_discrete(c(srv_data_1$t,1), t(srv_data_1[,-1]),
                                    t_optim, t(srv_data_2[,-1]), h = h)

    pfit_srv_data_2.compl <- c(w_optim) * complex(re = srv_data_2[,1], im = srv_data_2[,2])
    pfit_srv_data_2 <- data.frame(t = t_optim[-length(t_optim)], X1 = Re(pfit_srv_data_2.compl), X2 = Im(pfit_srv_data_2.compl))

    # TODO: Warping alignment gets all wonky.
    t_optim <- find_optimal_t_discrete(c(srv_data_1$t,1), t(srv_data_1[,-1]),
                                       c(pfit_srv_data_2$t,1), t(pfit_srv_data_2[,-1]),
                                       initial_t = t_optim, eps = eps)

    diff <- compute_distance(pfit_srv_data_2_old, pfit_srv_data_2, t_optim, closed)
    print(paste0('Diff : ', diff))
  }

  elastic_fp_dist <- compute_distance(srv_data_1, pfit_srv_data_2, t_optim, closed)

  pfit_curve1 <- cbind(t = c(srv_data_1$t, 1), pfit_curve1)
  pfit_curve2.compl <- c(w_optim * abs(w_optim)) * complex(re = pfit_curve2[,1], im = pfit_curve2[,2])
  pfit_curve2_aligned = data.frame(t = initial_t, t_optim = t_optim, X1 = Re(pfit_curve2.compl), X2 = Im(pfit_curve2.compl))

  aligned_curves_proc2d <- list("pfit_curve1" = pfit_curve1,
                         "pfit_curve2_aligned" = pfit_curve2_aligned,
                         "elastic_fp_dist" = elastic_fp_dist,
                         "w_optim" = w_optim,
                         "closed" = closed)

  class(aligned_curves_proc2d) <- "aligned_curves_proc2d"
  return(aligned_curves_proc2d)

}



# ' Test
# '
# '

align_curves_test <- function(data_curve1, data_curve2, delta = 0.01, eps = 0.01, h = 0.01){
  warning("This function is not fully implemented yet. At the moment only the alignment inside mean calculation works.")
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

  # Start
  efpfit_curve1 <- data_curve1
  efpfit_curve2 <- data_curve2
  efpfit_srv2 <- elasdics::get_srv_from_points(efpfit_curve2)
  t_optim <- c(efpfit_srv2$t, 1)

  diff <- Inf
  while(diff >= delta){
    efpfit_srv2_old <- efpfit_srv2

    # Procrustes alignment step
    aligned_curves_proc2d <- align_curves_proc2d(data_curve1, data_curve2, t_optim=t_optim, h=h)
    efpfit_curve2 <- aligned_curves_proc2d$pfit_curve2

    # Parametrizations alignment step
    aligned_curves <- elasdics::align_curves(efpfit_curve1, efpfit_curve2, eps=eps)
    efpfit_curve2 <- aligned_curves$data_curve2_aligned[,-1]
    efpfit_curve2 <- data.frame(t = efpfit_curve2$t_optim, efpfit_curve2[,-1])

    # Build SRV for distance calculation
    #t_optim <- efpfit_curve2$t_optim
    #efpfit_curve2$t <- efpfit_curve2$t_optim
    #efpfit_curve2$t_optim <- NULL

    print(efpfit_curve2)
    efpfit_srv2 <- get_srv_from_points(efpfit_curve2)
    diff <- compute_distance(efpfit_srv2_old, efpfit_srv2, t_optim, closed)
    print(diff)
  }

  efpfit_curve1 <- aligned_curves_proc2d$pfit_curve1
  efpfit_srv1 <- elasdics::get_srv_from_points(efpfit_curve1)
  efp_dist <- compute_distance(efpfit_srv1, efpfit_srv2, t_optim, closed)

  w_optim <- aligned_curves_proc2d$w_optim

  aligned_curves_proc2d <- list("efpfit_curve1" = efpfit_curve1,
                                "efpfit_curve2" = efpfit_curve2,
                                "elastic_fp_dist" = efp_dist,
                                "w_optim" = w_optim)

  class(aligned_curves_proc2d) <- "aligned_curves_proc2d"
  return(aligned_curves_proc2d)

}




align_curves_proc2d <- function(data_curve1, data_curve2, t_optim=NULL, h = 0.01){
  data_curve1 <- remove_duplicate(data_curve1)
  data_curve2 <- remove_duplicate(data_curve2)
  if(attributes(data_curve1)$points_rm | attributes(data_curve2)$points_rm){
    warning("Duplicated points in data curves have been removed!")
  }
  # input checking given parametrisation t
  if("t" %in% names(data_curve1)) check_param(data_curve1, closed=FALSE)
  if("t" %in% names(data_curve2)) check_param(data_curve2, closed=FALSE)

  # Prepare centered, unit-length curves
  length1 <- get_polygon_length(data_curve1)
  length2 <- get_polygon_length(data_curve1)
  pfit_curve1 <- elasdics::center_curve(data_curve1)/length1
  pfit_curve2 <- elasdics::center_curve(data_curve2)/length2
  # Prepare normalized SRV curves
  pfit_srv_data1 <- get_srv_from_points(pfit_curve1)
  pfit_srv_data2 <- get_srv_from_points(pfit_curve2)
  # Save parametrization
  t1 <- pfit_srv_data1$t
  t2 <- pfit_srv_data2$t
  if(!is.null(t_optim)){t2 <- t_optim[-length(t_optim)]}
  # To complex
  pfit_srv_data1.compl <- complex (re = pfit_srv_data1[,2], im = pfit_srv_data1[,3])
  pfit_srv_data2.compl <- complex (re = pfit_srv_data2[,2], im = pfit_srv_data2[,3])
  # TODO --- Approx. scalar product from piecewise constant srv curves. TODO: Exact solution
  pfit1.func <- sapply(seq(0,1,by=h), function(step){
    knot <- findInterval(step, t1, rightmost.closed=T)
    pfit_srv_data1.compl[knot]
  })
  pfit2.func <- sapply(seq(0,1,by=h), function(step){
    knot <- findInterval(step, t2, rightmost.closed=T)
    pfit_srv_data2.compl[knot]
  })
  w_optim <- h * t(Conj(pfit2.func)) %*% pfit1.func
  # Calculate pfit of curve2 onto curve1
  pfit_srv_data2.compl <- c(w_optim) * pfit_srv_data2.compl
  pfit_srv_data2 <- data.frame(t=t2, X1 = Re(pfit_srv_data2.compl), X2 = Im(pfit_srv_data2.compl))
  pfit_curve2 <- data.frame(t = c(t2,1), elasdics::get_points_from_srv(pfit_srv_data2))
  pfit_curve2 <- elasdics::center_curve(pfit_curve2)
  # Calculate procrustes distance
  full_proc2d_dist <- compute_distance(pfit_srv_data1, pfit_srv_data2, c(t2,1), closed=FALSE)
  # Build output
  pfit_curve1 <- cbind(t = c(pfit_srv_data1$t, 1), pfit_curve1)

  aligned_curves_proc2d <- list("pfit_curve1" = pfit_curve1,
                                "pfit_curve2" = pfit_curve2,
                                "full_proc2d_dist" = full_proc2d_dist,
                                "w_optim" = w_optim,
                                "closed" = FALSE)

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


