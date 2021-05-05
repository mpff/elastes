#' Compute the loss function for discrete input
#' @inheritParams align_curves
#' @param srv_data_1 srv of curve1
#' @param srv_data_2 srv of curve2
#' @param t parametrisation
#' @noRd

get_loss_discrete <- function(t, srv_data_1, srv_data_2){
  stopifnot(length(t) == nrow(srv_data_2) + 1)
  #unfold srv_data
  s <- srv_data_2$t
  q <- t(srv_data_2[,-1])
  r <- srv_data_1$t
  p <- t(srv_data_1[,-1])

  delta_s <- diff(c(s,1))
  breaks <- sort(unique(c(t,r)))
  idx <- sapply(breaks, function(i) max(which(r <= i)))
  p_breaks <- sapply(idx, function(i) p[,i])

  integrals <- sapply(1:(length(t)-1), function(i){
    idx_i <- breaks >= t[i] & breaks <= t[i+1]
    delta_r <- diff(breaks[idx_i])
    if(length(delta_r) == 0){ value <- 0 } else {
      p_i <- p_breaks[,idx_i]
      p_i[,-ncol(p_i)]
      integrand <- as.vector(t(q[,i])%*%p_i[,-ncol(p_i)])
      value <- delta_s[i]*t(delta_r)%*%(ifelse(integrand >= 0, integrand, 0)^2)
    }
    value
  })

  sum(sqrt(integrals))
}
