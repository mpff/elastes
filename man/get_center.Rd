% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_center}
\alias{get_center}
\title{Calculate the center of a curve}
\usage{
get_center(curve)
}
\arguments{
\item{curve}{a \code{data.frame} with observed points in each row. Each
variable is one coordinate direction. If there is a variable \code{t},
\code{t_optim} or \code{id}, it is treated as the time parametrization, not
as an additional coordinate.}
}
\value{
The average of observed points in \code{curve}.
}
\description{
Calculate the center of a curve
}
