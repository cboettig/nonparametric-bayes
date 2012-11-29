#' gp_fit fits a Gaussian process to the observed data
#' 
#' API not yet stable! Arguments to this function will probably be updated to be much more generic
#' @param obs a data frame of observations with columns obs$x and obs$y
#' @param X the desired points over which to predict
#' @param pars a named numeric specifying "sigma_n" for the (additive) and "l" for the covariance length scale
#' @param method select the method to use
#' @param out_var optionally force kernlab method to generate Cf using a different variance than sigma_n^2 assumed in the fit
#' @return a list with items "mu", the expected Y values at X (mean of the posterior Gaussian process), 
#' Sigma, the covariance matrix for the posterior Gaussian process, and 
#' "loglik", the log likelihood of observering the given data under the process, marginalized over the prior
#' @details so far treates the prior as mean 0 and covariance given by cov(X)
#' @export
#' @examples 
#' X <- seq(-5,5,len=50)
#' obs <- data.frame(x = c(-4, -3, -1,  0,  2),
#'                   y = c(-2,  0,  1,  2, -1))
#' gp <- gp_fit(obs, X, c(sigma_n=0.8, l=1))
#' 
#' df <- data.frame(x=X, y=gp$mu, ymin=(gp$mu-2*sqrt((gp$var))), ymax=(gp$mu+2*sqrt((gp$var))))
#' require(ggplot2)
#' ggplot(df)  + geom_ribbon(aes(x,y,ymin=ymin,ymax=ymax), fill="gray80") + geom_line(aes(x,y)) + geom_point(data=obs, aes(x,y))

gp_fit <- function(obs, X, pars=c(sigma_n=1, tau=1, l=1), method=c("direct", "sequential", "kernlab", "cholesky", "conditional"), fit = TRUE, out_var=NULL, ...){
  
  method <- match.arg(method)
  sigma_n <- pars["sigma_n"]
  l <- pars["l"]
  tau <- pars["tau"]
  llik <- NA 
  
  ## Parameterization-specific
  SE <- function(Xi,Xj, l=l) tau * exp(-0.5 * (Xi - Xj) ^ 2 / l ^ 2)
  cov <- function(X, Y) outer(X, Y, SE, l) 
  n <- length(obs$x)
  K <- cov(obs$x, obs$x)
  I <-  diag(1, n)
  
  
  ## Cholesky simultaneous method.  
  if(method=="cholesky"){
    L <- chol(K + sigma_n ^ 2 * I)
    alpha <- solve(t(L), solve(L, obs$y))
    loglik <- -.5 * t(obs$y) %*% alpha - sum(log(diag(L))) - n * log(2 * pi) / 2
    k_star <- cov(obs$x, X)
    mu <- t(k_star) %*% alpha
    v <- solve(L, k_star)
    var <- diag( cov(X,X) - t(v) %*% v )
  }
  
  ## Direct method 
  else if(method=="direct"){
    cov_xx_inv <- solve(K + sigma_n ^ 2 * I)
    Ef <- cov(X, obs$x) %*% cov_xx_inv %*% obs$y
    Cf <- cov(X, X) - cov(X, obs$x)  %*% cov_xx_inv %*% cov(obs$x, X)
    llik <- -.5 * t(obs$y) %*% cov_xx_inv %*% obs$y - 0.5 * log(det(cov_xx_inv)) - n * log(2 * pi) / 2
  }
  
  
  
  ## Direct method, conditional on passing through 0,0 (quick hack)
  else if(method=="conditional"){
    obs <- data.frame(x=c(0, obs$x), y = c(0, obs$y))
    n <- length(obs$x)
    K <- cov(obs$x, obs$x)
    sigmaI <-  sigma_n ^ 2 * diag(1, n)
    sigmaI[1,1] <- 1e-5 # no noise about zero.  
    cov_xx_inv <- solve(K + sigmaI)
    Ef <- cov(X, obs$x) %*% cov_xx_inv %*% obs$y
    Cf <- cov(X, X) - cov(X, obs$x)  %*% cov_xx_inv %*% cov(obs$x, X)
    llik <- -.5 * t(obs$y) %*% cov_xx_inv %*% obs$y - 0.5 * log(det(cov_xx_inv)) - n * log(2 * pi) / 2
  }
  
  ### Sequential method, using function recursion
  else if(method=="sequential"){
    mmult <- function(x,y){
      if(length(x) == 1){
        x <- as.numeric(x) 
        x * y
      } else if(length(y) == 1){ 
        y <- as.numeric(y)
        x * y
      }  else 
        x %*% y
    }
    mu <- numeric(length(obs$y))
    y <- obs$y
    x <- obs$x
    C_seq <- function(X, X_prime, i){
      if(i <= 1)
        cov(X, X_prime) - mmult(cov(X,x[i]), cov(x[i], X_prime)) / as.numeric( cov(x[i], x[i]) + sigma_n^2)
      else
        C_seq(X, X_prime,   i-1) - mmult(C_seq(X,x[i],   i-1), C_seq(x[i], X_prime,   i-1)) / as.numeric( C_seq(x[i], x[i],   i-1)  + sigma_n^2  )
    }
    mu_seq <- function(X, i){
      if(i <= 1)
        cov(x[i], X) * (y[i]-mu[i]) / as.numeric( cov(x[i], x[i]) + sigma_n^2)
      else
        mu_seq(X, i-1) + C_seq(x[i], X,  i-1) * (y[i]-mu[i]) / as.numeric( C_seq(x[i], x[i], i-1)  + sigma_n^2 )
    }
    Ef <- t(mu_seq(X, length(obs$x)))
    Cf <- C_seq(X, X, length(obs$x))
  }
  
  
  
  ## Kernlab method
  else if(method=="kernlab"){
  var <- sigma_n^2
  if(!fit)
    gp <- gausspr(obs$x, obs$y, kernel="rbfdot", kpar=list(sigma=1/(2*l^2)), fit=FALSE, scaled=FALSE, var=var) # fixed length-scale
  else
    gp <- gausspr(obs$x, obs$y, kernel="rbfdot", scaled=FALSE, var = var)
  Ef <- predict(gp, X)
  
  if(!is.null(out_var)) # permit an external variance
    var <- out_var
  
  ##Manually get the covariance, $K(x_*, x_*) - K(x_*, x)(K(x, x) - \sigma_n^2\mathbb{I})^{-1}K(x,x_*)$. 
  Inv <- solve(kernelMatrix(kernelf(gp), xmatrix(gp)) + diag(rep(var, length = length(xmatrix(gp)))))
  Cf <- kernelMatrix(kernelf(gp), as.matrix(X)) - 
    kernelMult(kernelf(gp), as.matrix(X), xmatrix(gp),  Inv) %*% 
    kernelMatrix(kernelf(gp), xmatrix(gp), as.matrix(X))
  llik <- -.5 * t(obs$y) %*% Inv %*% obs$y - 0.5 * log(det(Inv)) - n * log(2 * pi) / 2
  
  }
    
  out <- list(Ef = Ef, Cf = Cf, llik = llik, obs = obs, X = X, K = K, pars = pars)
  class(out) = "gpfit"
  out
  }



update.gpfit <- function(gp, obs, sigma){
  
}


