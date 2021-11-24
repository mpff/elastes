#' Fitting function for open curves
#' @name fit_mean_closed
#' @description Fits an elastic mean for open curves. Is usually called from
#' \code{\link{compute_elastic_mean}}.
#' @param srv_data_curves list of \code{data.frame}s with srv vectors in each row.
#' Usually a result of a call to \code{\link{get_srv_from_points}}
#' @param knots set of knots for the mean spline curve
#' @param max_iter maximal number of iterations
#' @param type if "smooth" linear srv-splines are used which results in a differentiable mean curve
#' @param eps the algorithm stops if L2 norm of coefficients changes less
#' @param pen_factor penalty factor forcing the mean to be closed
#' if "polygon" the mean will be piecewise linear.
#' @return a \code{list}
#' with entries
#'   \item{type}{"smooth" or "polygon"}
#'   \item{coefs}{\code{coefs} srv spline coefficients of the estimated mean}
#'   \item{knots}{spline knots}
#'   \item{t_optims}{optimal parametrisation}
#'   \item{shift_idxs}{index of the starting point of the closed curve after alignment}

fit_mean_closed <- function(srv_data_curves, knots, max_iter, type, weights,
                            eps, pen_factor){
  #pre align srv curves
  pre_aligned_curves <- lapply(srv_data_curves, function(data_curve_1){
    aligned_curves <- lapply(srv_data_curves, function(data_curve_2){
      pre_align_closed_srv_curves(data_curve_1, data_curve_2)
    })
    sum_dist <- sum(sapply(aligned_curves, function(x) x$dist)^2)
    shift_idxs <-  sapply(aligned_curves, function(x) x$shift_idx)
    aligned_curves <- lapply(aligned_curves, function(x) x[[1]])
    list(aligned_curves, "sum_dist" = sum_dist, "shift_idxs" = shift_idxs)
  })
  idx <- which.min(sapply(pre_aligned_curves, function(x) x$sum_dist))
  srv_data_curves <- pre_aligned_curves[idx][[1]][[1]]
  shift_idxs <- pre_aligned_curves[idx][[1]]$shift_idxs

  #initial param and coeffs
  t_optims <- lapply(srv_data_curves, function(srv_data_curve){
    c(srv_data_curve$t, 1)
  })
  model_data <- get_model_data(t_optims, srv_data_curves, knots, type, weights)
  coefs <- apply(model_data[,-(1:2)], 2, function(q_m_x_long){
    q_m_x_long[!is.finite(q_m_x_long)] <- NA
    coef(lm(q_m_x_long ~ -1 + make_design(model_data[,1], knots = knots,
                                          closed = TRUE, type = type),
            weights = model_data[,2]))
  })
  if(max_iter == 0){
    return(list("type" = type, "coefs" = coefs, "knots" = knots, "t_optims" = t_optims,
                "shift_idxs" = shift_idxs))
  }

  #iterate L_2 mean fit and warping
  for (i in 1:max_iter){
    ## warping fit, update t_optims
    if(type == "smooth"){
      pfun <- function(t){
        t(make_design(t, knots = knots,
                      closed = TRUE, type = type) %*% coefs)
      }

      t_optims <- lapply(1:length(srv_data_curves), function(j){
        t_optim <- find_optimal_t(srv_curve = pfun,
                                  s = c(srv_data_curves[[j]]$t, 1),
                                  q = t(srv_data_curves[[j]][,-1]),
                                  initial_t = t_optims[[j]],
                                  eps = eps*100/i)
        attr(t_optim, "dist_to_mean") <- attr(t_optim, "dist")
        attr(t_optim, "dist") <- NULL
        t_optim
      })
    } else {
      t_optims <- lapply(1:length(srv_data_curves), function(j){
        t_optim <- find_optimal_t_discrete(r = knots,
                                           p = t(coefs),
                                           s = c(srv_data_curves[[j]]$t, 1),
                                           q = t(srv_data_curves[[j]][,-1]),
                                           initial_t = t_optims[[j]])
        attr(t_optim, "dist_to_mean") <- attr(t_optim, "dist")
        attr(t_optim, "dist") <- NULL
        t_optim
      })
    }

    ############################################################################
    #L2 fit
    coefs_old <- coefs
    model_data <- get_model_data(t_optims, srv_data_curves, knots, type, weights)
    coefs <- apply(model_data[,-(1:2)], 2, function(q_m_x_long){
      q_m_x_long[!is.finite(q_m_x_long)] <- NA
      coef(lm(q_m_x_long ~ -1 + make_design(model_data[,1], knots = knots,
                                            closed = TRUE, type = type),
              weights = model_data[,2]))
    })

    #stop if coefficient don't change much anymore
    stop_crit <- sum((coefs - coefs_old)^2)/sum(coefs^2)

    #one step of penalised optimisation
    coefs <- one_step_grad(coefs, pen_factor = pen_factor*i, model_data = model_data,
                           knots = knots, type = type)

    if(stop_crit < eps){
      #one step of penalised optimisation
      coefs <- one_step_grad(coefs, pen_factor = Inf, model_data = model_data,
                             knots = knots, type = type)
      rownames(coefs) <- NULL
      colnames(coefs) <- colnames(srv_data_curves[[1]][,-1])
      return(list("type" = type, "coefs" = coefs, "knots" = knots, "t_optims" = t_optims,
                  "shift_idxs" = shift_idxs))
    }
  }
  warning("Stopping criteria eps has not been reached! Consider more iterations max_iter")
  rownames(coefs) <- NULL
  colnames(coefs) <- colnames(srv_data_curves[[1]][,-1])
  return(list("type" = type, "coefs" = coefs, "knots" = knots, "t_optims" = t_optims,
              "shift_idxs" = shift_idxs))
}

pre_align_closed_srv_curves <- function(srv_data_1, srv_data_2){
  if(identical(srv_data_1, srv_data_2)){
    return(list(srv_data_2, "dist" = 0, "shift_idx" = 1))
  }
  #pre alignment
  t <- c(srv_data_2$t, 1)
  dist <- sapply(1:length(srv_data_2$t), function(i){
    t_new <- t - t[i]
    compute_distance(srv_data_1, srv_data_2, t_new, closed = TRUE)
  })
  offset <- srv_data_2$t[which.min(dist)]
  new_t <- srv_data_2$t - offset
  srv_data_2$t <- ifelse(new_t < 0, new_t + 1, new_t)
  srv_data_2 <- srv_data_2[order(srv_data_2$t),]

  return(list(srv_data_2, "dist" = min(dist), "shift_idx" = which.min(dist)))
}


one_step_grad <- function(coefs, pen_factor, model_data, knots, type){
  coefs_long <- as.vector(coefs)
  grad_pen <- get_grad_penalty_fun(coefs_long, pen_factor, model_data, knots, type)
  get_loss_pen_1 <- function(x){
    get_loss_pen(coefs_long + x*grad_pen, pen_factor = pen_factor,
                 model_data = model_data, knots = knots, type = type)
  }
  x_optim <- optimise(get_loss_pen_1, interval = c(-1,1))$minimum
  coefs <- matrix(coefs_long + x_optim*grad_pen, ncol = 2)
  coefs
}


get_loss_pen <- function(coefs_long, pen_factor, model_data, knots, type){
  model_data[!is.finite(as.matrix(model_data))] <- NA
  model_data <- na.omit(model_data)
  m <- model_data[,1]
  coefs <- matrix(coefs_long,ncol = 2)
  pfun <- function(t){t(make_design(t, knots, closed = TRUE, type = type) %*% coefs)}
  L2_loss <- sum((t(pfun(m)) - model_data[, -(1:2)])^2)
  if(type == "smooth"){
    penalty <- sum((srvf_to_curve(c(0,1), pfun)[,2])^2)
  } else if (type == "polygon"){
    srv_data <- data.frame("t" = knots[-length(knots)], coefs)
    points <- get_points_from_srv(srv_data)
    penalty <- sum((points[nrow(points),])^2)
  }

  L2_loss/pen_factor + penalty
}

get_grad_penalty_fun <- function(coefs_long, pen_factor, model_data, knots, type){
  model_data[!is.finite(as.matrix(model_data))] <- NA
  model_data <- na.omit(model_data)
  m <- model_data[,1]
  coefs <- matrix(coefs_long, ncol = 2)
  design_mat <- make_design(m, knots, closed = TRUE, type = type)
  left_part <- t(design_mat)%*%design_mat%*%coefs
  right_part <- t(design_mat)%*%as.matrix(model_data[,-(1:2)])
  grad_part_linear <- as.vector(left_part - right_part)


  penalty_fun <- function(coefs_long){
    coefs <- matrix(coefs_long, ncol = 2)
    if(type == "smooth"){
      pfun <- function(t){t(make_design(t, knots, closed = TRUE, type = type) %*% coefs)}
      penalty <- sum((srvf_to_curve(c(0,1), pfun)[,2])^2)
    } else if (type == "polygon"){
      srv_data <- data.frame("t" = knots[-length(knots)], coefs)
      points <- get_points_from_srv(srv_data)
      penalty <- sum((points[nrow(points),])^2)
    }
    penalty
  }

  grad_part_linear/pen_factor + numDeriv::grad(penalty_fun, coefs_long, method.args=list(r = 2))
}
