\name{gp_transition_matrix}
\alias{gp_transition_matrix}
\title{Determine the transition matrix given a Guassian process}
\usage{
  gp_transition_matrix(gp, h_grid = NULL)
}
\arguments{
  \item{gp}{a gpfit type object, a list with the GP mean
  gp$Ef, Covariance gp$Cf on grid gp$X}

  \item{h_grid}{the harvest grid over which to compute.
  Assumed to match the stock grid gp$X if not given.}
}
\value{
  a list of matrices corresponding to the state transition
  matrix at each harvest level h
}
\description{
  Determine the transition matrix given a Guassian process
}

