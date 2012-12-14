\name{par_F}
\alias{par_F}
\title{calculate transition matrix based on parameterized function}
\usage{
  par_F(h, f, p, x_grid, sigma_g)
}
\arguments{
  \item{f}{the transtion function, a function of state x,
  harvest h, and parameters p}

  \item{p}{the parameters p of the transition function}

  \item{x_grid}{the grid over which to solve}

  \item{h_grid}{the harvest grid over which to compute.
  Assumed to match the stock grid gp$X if not given.}

  \item{sigma_g}{the scale of the noise}
}
\value{
  a list of matrices corresponding to the state transition
  matrix at each harvest level h
}
\description{
  calculate transition matrix based on parameterized
  function
}

