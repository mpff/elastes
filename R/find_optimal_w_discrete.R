#' @title Finds optimal rotation and scaling alignment for discrete open curves
#' @description Finds optimal aligned rotation and scaling for srv curve q to srv curve p using
#' coordinate wise optimisation.
#' @param r time points for p, first has to be 0, last has to be 1
#' @param p square root velocity vectors, one less than time points in r
#' @param s time points for q, first has to be 0, last has to be 1
#' @param q square root velocity vectors, one less than time points in s
#' @param h stepsize for pw constant integral approximation (TODO exact solution)
#' @return optimal rotation and scaling alignment for q
#' has the distance of the observation to the srv_curve as an attribute

find_optimal_w_discrete <- function(r, p, s, q, h = 0.01){

  p.compl <- complex(re = p[1,], im = p[2,])
  q.compl <- complex(re = q[1,], im = q[2,])

  # Ultra hacky for now <-- REPLACE!
  p.func <- sapply(seq(0,1,by=h), function(step){
    knot <- findInterval(step, r, rightmost.closed=T)
    p.compl[knot]
  })

  q.func <- sapply(seq(0,1,by=h), function(step){
    knot <- findInterval(step, s, rightmost.closed=T)
    q.compl[knot]
  })

  w_optim <- h * t(Conj(q.func)) %*% p.func

  w_optim
}
