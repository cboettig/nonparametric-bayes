\name{gp_predict}
\alias{gp_predict}
\title{predict the expected values and posterior distributions of the Gaussian Process}
\usage{
  gp_predict(gp, x_predict, burnin = 0, thin = 1,
    covs = FALSE)
}
\arguments{
  \item{x_predict}{the values at which we desire
  predictions}

  \item{covs}{if TRUE, return covariances instead of
  variances (set high thinning as this is memory
  intensive)}
}
\description{
  predict the expected values and posterior distributions
  of the Gaussian Process
}

