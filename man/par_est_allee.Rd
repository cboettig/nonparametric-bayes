\name{par_est_allee}
\alias{par_est_allee}
\title{MLE estimate of parameters given by function f}
\usage{
  par_est_allee(obs, f, p)
}
\arguments{
  \item{obs}{the observed data, two columns giving x_t,
  x_t+1 respectively}
}
\value{
  a list with the f given, the MLE estimated parameters,
  and estimated noise level
}
\description{
  MLE estimate of parameters given by function f
}

