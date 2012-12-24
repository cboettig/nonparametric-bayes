\name{sim_obs}
\alias{sim_obs}
\title{Simulate observed data for training the GP}
\usage{
  sim_obs(Xo, z_g, f, p, Tobs = 40, seed = 1, nz = 10,
    harvest = sort(rep(seq(1, 2, length = 5), 8)))
}
\arguments{
  \item{Xo}{the initial population size}

  \item{z_g}{a function of no arguments providing the
  random multiplicative growth noise}

  \item{f}{a function of (x, h, p) providing the growth in
  population x under harvest h given pars p}

  \item{p}{the parameters of function f}

  \item{Tobs}{the number of timesteps on which to observe
  the system}

  \item{seed}{the random number seed for the simulation}

  \item{nz}{the number of observations at 0 (a hack, should
  just condiition on 0,0)}

  \item{harvest}{a sequence of harvest values that occur
  durning the simulation}
}
\description{
  Simulate observed data for training the GP
}
\details{
  varying harvest values allow the system to explore the
  state space, making for better training data.
}

