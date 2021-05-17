#' @title Finds optimal alignment for discrete open curves
#' @description Finds optimal aligned time points for srv curve q to  srv curve p using
#' coordinate wise optimisation.
#' @param r time points for p, first has to be 0, last has to be 1
#' @param p square root velocity vectors, one less than time points in r
#' @param s time points for q, first has to be 0, last has to be 1
#' @param q square root velocity vectors, one less than time points in s
#' @param initial_t starting value for the optimisation algorithm
#' @param eps convergence tolerance
#' @return optimal time points for q, without first value 0 and last value 1
#' optimal time points have the distance of the observation to the srv_curve as an attribute

find_optimal_t_discrete <- function(r, p, s, q, initial_t = s, eps = 10^-3){
  t <- initial_t
  #adjust convergence criterium for number of points
  eps <- eps/length(t)
  delta <- eps
  while(delta >= eps){
    t_old <- t
    # optimise all even indices
    idx_even <- 2*(1:ceiling((length(t) - 2)/2))
    t[idx_even] <- sapply(idx_even, function(i){
      optimise_one_coord_analytic(t, i, r, p, s, q)[i]
    })
    #optimise all odd indices
    if(length(t) > 3){
      idx_odd <- 2*(1:ceiling((length(t) - 3)/2)) + 1
      t[idx_odd] <- sapply(idx_odd, function(i){
        optimise_one_coord_analytic(t, i, r, p, s, q)[i]
      })
    }
    delta <- max(abs(t_old - t))
  }

  t_optim <- t
  dist <- compute_distance(data.frame("t" = r[-length(r)], t(p)),
                           data.frame("t" = s[-length(s)], t(q)),
                           t_optim = t_optim, closed = FALSE)
  attr(t_optim, "dist") <- dist
  t_optim
}

#' Does optimisation in one parameter direction
#' @inheritParams find_optimal_t_discrete
#' @param t current time points, first has to be 0, last has to be 1
#' @param i index of t that should be updated
#' @return optimal time points for q with respect to optimisation only in the i-th
#' coordinate direction

optimise_one_coord_analytic <-function(t, i, r, p, s, q){
  # Find time points in the current interval
  idx_r <- which(r > t[i - 1] & r < t[i + 1])
  r_restr <- c(t[i-1], r[idx_r], t[i + 1])

  # Find corresponding srv vectors
  idx_p <- c(max(which(r <= t[i - 1])), idx_r)
  p_restr <- p[,idx_p, drop = FALSE]

  cross_plus <- function(p, q){
    max(0, crossprod(p,q))
  }

  # Optimise analytically
  t_optim_k <- sapply(1:(length(r_restr) - 1), function(j){

    left_k <- which(1:length(r_restr) < j)
    right_k <- which(1:(length(r_restr) - 1) > j)
    H_left_k <- as.numeric( sapply(left_k, function(k) cross_plus(p_restr[,k], q[, i -1])^2*(r_restr[k+1] - r_restr[k])) )
    H_right_k <- as.numeric( sapply(right_k, function(k) cross_plus(p_restr[,k], q[, i])^2*(r_restr[k+1] - r_restr[k])) )

    # determine parts of loss function
    A <- (s[i] - s[i-1])*cross_plus(p_restr[,j], q[, i-1])^2
    B <- (s[i] - s[i-1])*sum(H_left_k) - A * r_restr[j]
    C <- - (s[i + 1] - s[i])*cross_plus(p_restr[,j], q[, i])^2
    D <- (s[i + 1] - s[i])*sum(H_right_k) - C * r_restr[j + 1]

    t <- (C^2*B - A^2*D) / (A*C*(A-C))

    ##treat special cases seperatly
    if(A == 0) t <- ifelse(C <= 0, r_restr[j], r_restr[j + 1])
    if(C == 0) t <- ifelse(A <= 0, r_restr[j], r_restr[j + 1])
    if( t < r_restr[j]  ) t <- r_restr[j]
    if( t > r_restr[j+1]  ) t <- r_restr[j+1]

    value <- sqrt(A*t + B) + sqrt(C*t + D)

    c("value" = value, "t_i" = t)
  })

  t[i] <- t_optim_k[2, which.max(t_optim_k[1,])]
  t
}

