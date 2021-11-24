#' @title Finds optimal rotation and scaling alignment for discrete open curves
#' @description Finds optimal aligned rotation and scaling for srv curve q to srv curve p using
#' coordinate wise optimisation.
#' @param r time points for p, first has to be 0, last has to be 1
#' @param p square root velocity vectors, one less than time points in r
#' @param s time points for q, first has to be 0, last has to be 1
#' @param q square root velocity vectors, one less than time points in s
#' @return optimal rotation and scaling alignment for q
#' has the distance of the observation to the srv_curve as an attribute

find_optimal_w_discrete <- function(r, p, s, q){

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
