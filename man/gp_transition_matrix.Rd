\name{gp_transition_matrix}
\alias{gp_transition_matrix}
\title{Determine the transition matrix given a Guassian process}
\usage{
  gp_transition_matrix(Ef, Sf, x_grid, h_grid = NULL)
}
\arguments{
  \item{Ef}{GP mean}

  \item{Sf}{GP standard deviations (sqrt(diag(Cf))}

  \item{x_grid}{x grid over which to define transition
  matrices}

  \item{h_grid}{the harvest grid over which to compute.
  Assumed to match the stock grid x_grid if not given.}
}
\value{
  a list of matrices corresponding to the state transition
  matrix at each harvest level h
}
\description{
  Determine the transition matrix given a Guassian process
}

