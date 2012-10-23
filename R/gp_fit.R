

gp_fit <- function(){
  
  ## Parameterization-specific
  SE <- function(Xi,Xj, l) exp(-0.5 * (Xi - Xj) ^ 2 / l ^ 2)
  cov <- function(X, Y) outer(X, Y, SE, l) 
   
  
  n <- length(obs$x)
  K <- cov(obs$x, obs$x)
  I <-  diag(1, n)
  
  
  ## Cholksy method
  L <- chol(K + sigma_n^2 * I)
  alpha <- solve(t(L), solve(L, obs$y))
  k_star <- cov(obs$x, X)
  Y <- t(k_star) %*% alpha
  v <- solve(L, k_star)
  Var <- cov(X,X) - t(v) %*% v
  loglik <- -.5 * t(obs$y) %*% alpha - sum(log(diag(L))) - n * log(2 * pi) / 2
  
  ## Direct method 
  cov_xx_inv <- solve(K + sigma_n^2*I)
  Ef <- cov(X, obs$x) %*% cov_xx_inv %*% obs$y
  Cf <- cov(X, X) - cov(X, obs$x)  %*% cov_xx_inv %*% cov(obs$x, X)
  
  return(list(mu = Y, Sigma = Var, loglik=loglik, Ef = Ef, Cf = Cf)
  }