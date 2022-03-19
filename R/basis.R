# Get inner and outer knots.
get_knots <- function(knots, type){
  deg <- ifelse(type == "smooth", 1, 0)
  knotl = 1 / ( length(knots) - 1 )  # mean length of a knot
  c(rep(-knotl, deg), knots, rep(1+knotl, deg))
}


# Get b-spline Gram matrix.
get_gram_matrix <- function(knots, type){
  knots = get_knots(knots, type)
  if( type == "polygon"){
    1/(length(knots) - 1) * diag(length(knots) - 1)  # Note: ONLY CORRECT FOR EQUIDISTANT KNOTS!!!
  } else {
    osb_smooth = orthogonalsplinebasis::SplineBasis(knots, order = 2)  # degree = order - 1
    orthogonalsplinebasis::GramMatrix(osb_smooth)
  }
}


# Design matrix of mean basis.
make_design <- function(t, knots, type, closed = FALSE) {
  deg <- ifelse(type == "smooth", 1, 0)
  knots <- get_knots(knots, type)
  design_mat <- splines::splineDesign(knots = knots, x = t, ord = deg + 1)
  if(closed == TRUE & type == "smooth"){
    design_mat[,1] <- design_mat[,1] + design_mat[, ncol(design_mat)]
    design_mat <- design_mat[,-ncol(design_mat)]
  }
  design_mat
}
