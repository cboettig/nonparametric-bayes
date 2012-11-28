\name{f_transition_matrix}
\alias{f_transition_matrix}
\title{calculate transition matrix based on parameterized function}
\usage{
  f_transition_matrix(f, p, x_grid, h_grid = NULL)
}
\arguments{
  \item{f}{the transtion function, a function of state x,
  harvest h, and parameters p}

  \item{p}{the parameters p of the transition function}

  \item{x_grid}{the grid over which to solve}

  \item{h_grid}{the harvest grid over which to compute.
  Assumed to match the stock grid gp$X if not given.}
}
\value{
  a list of matrices corresponding to the state transition
  matrix at each harvest level h
}
\description{
  calculate transition matrix based on parameterized
  function
}

