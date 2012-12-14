\name{gp_F}
\alias{gp_F}
\title{Determine the transition matrix given a Guassian process}
\usage{
  gp_F(h, Ef, V, x_grid)
}
\arguments{
  \item{Ef}{GP mean}

  \item{V}{GP vector of variances}

  \item{x_grid}{x grid over which to define transition
  matrices}
}
\value{
  a list of matrices corresponding to the state transition
  matrix at each harvest level h
}
\description{
  Determine the transition matrix given a Guassian process
}

