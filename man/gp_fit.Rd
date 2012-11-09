\name{gp_fit}
\alias{gp_fit}
\title{gp_fit fits a Gaussian process to the observed data}
\usage{
  gp_fit(obs, X, pars = c(sigma_n = 1, l = 1))
}
\arguments{
  \item{obs}{a data frame of observations with columns
  obs$x and obs$y}

  \item{X}{the desired points over which to predict}

  \item{pars}{a named numeric specifying "sigma_n" for the
  (additive) and "l" for the covariance length scale}
}
\value{
  a list with items "mu", the expected Y values at X (mean
  of the posterior Gaussian process), Sigma, the covariance
  matrix for the posterior Gaussian process, and "loglik",
  the log likelihood of observering the given data under
  the process, marginalized over the prior
}
\description{
  API not yet stable! Arguments to this function will
  probably be updated to be much more generic
}
\details{
  so far treates the prior as mean 0 and covariance given
  by cov(X)
}
\examples{
X <- seq(-5,5,len=50)
obs <- data.frame(x = c(-4, -3, -1,  0,  2),
                  y = c(-2,  0,  1,  2, -1))
gp <- gp_fit(obs, X, c(sigma_n=0.8, l=1))

df <- data.frame(x=X, y=gp$mu, ymin=(gp$mu-2*sqrt((gp$var))), ymax=(gp$mu+2*sqrt((gp$var))))
require(ggplot2)
ggplot(df)  + geom_ribbon(aes(x,y,ymin=ymin,ymax=ymax), fill="gray80") + geom_line(aes(x,y)) + geom_point(data=obs, aes(x,y))
}

