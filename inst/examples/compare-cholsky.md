

## Compare Cholsky and Direct methods


```r
gp_fit_direct <- function(obs, X, pars=c(sigma_n=1, l=1)){
  
  sigma_n <- pars["sigma_n"]
  l <- pars["l"]
  
  ## Parameterization-specific
  SE <- function(Xi,Xj, l) exp(-0.5 * (Xi - Xj) ^ 2 / l ^ 2)
  cov <- function(X, Y) outer(X, Y, SE, l) 
  
  n <- length(obs$x)
  K <- cov(obs$x, obs$x)
  I <-  diag(1, n)
  
  ## Direct method 
  cov_xx_inv <- solve(K + sigma_n ^ 2 * I)
  Ef <- cov(X, obs$x) %*% cov_xx_inv %*% obs$y
  Cf <- cov(X, X) - cov(X, obs$x)  %*% cov_xx_inv %*% cov(obs$x, X)
  llik <- -.5 * t(obs$y) %*% cov_xx_inv %*% obs$y - 0.5 * log(det(cov_xx_inv)) - n * log(2 * pi) / 2
  
  out <- list(mu = Ef, var = Cf, loglik=llik)
  class(out) = "gp"
  out
}
```



```r
gp_fit_cholsky <- function(obs, X, pars=c(sigma_n=1, l=1)){
  
  sigma_n <- pars["sigma_n"]
  l <- pars["l"]
  
  ## Parameterization-specific
  SE <- function(Xi,Xj, l) exp(-0.5 * (Xi - Xj) ^ 2 / l ^ 2)
  cov <- function(X, Y) outer(X, Y, SE, l) 
  
  n <- length(obs$x)
  K <- cov(obs$x, obs$x)
  I <-  diag(1, n)
  
  ## Cholksy method -- Consider implementing in C
  L <- chol(K + sigma_n ^ 2 * I)
  alpha <- solve(t(L), solve(L, obs$y))
  loglik <- -.5 * t(obs$y) %*% alpha - sum(log(diag(L))) - n * log(2 * pi) / 2
  k_star <- cov(obs$x, X)
  mu <- t(k_star) %*% alpha
  v <- solve(L, k_star)
  var <- diag( cov(X,X) - t(v) %*% v )

  
  ## Pt by pt following Rassmussen & WIlliams, but there's no reason to do it this way
#  tmp <- lapply(X, function(x_i){  
#    k_star <- cov(obs$x, x_i)
#    Y <- t(k_star) %*% alpha
#    v <- solve(L, k_star)
#    Var <- cov(x_i,x_i) - t(v) %*% v
#    list(Y=Y, Var=Var)
#  })
#  mu <- sapply(tmp,`[[`, "Y")
#  var <- sapply(tmp,`[[`, "Var")
    
  out <- list(mu = mu, var = var, loglik=loglik)
  class(out) = "gp"
  out
}
```


Run the comparison


```r
X <- seq(-5,5,len=50)
obs <- data.frame(x = c(-4, -3, -1,  0,  2),
                  y = c(-2,  0,  1,  2, -1))
system.time(gp <- gp_fit_cholsky(obs, X, c(sigma_n=0.8, l=1)))
```

```
   user  system elapsed 
  0.004   0.000   0.005 
```

```r
system.time(gp_d <- gp_fit_direct(obs, X, c(sigma_n=0.8, l=1)))
```

```
   user  system elapsed 
  0.000   0.000   0.002 
```


Plot the results


```r
require(reshape2)
df <- data.frame(x=X, y=gp$mu, ymin=(gp$mu-2*sqrt(gp$var)), ymax=(gp$mu+2*sqrt(gp$var)))
df2 <- data.frame(x=X, y=gp_d$mu, ymin=(gp_d$mu-2*sqrt(diag(gp_d$var))), ymax=(gp_d$mu+2*sqrt(diag(gp_d$var))))
df2 <- melt(df2, id="x")
require(ggplot2)
ggplot(df)  + geom_ribbon(aes(x,y,ymin=ymin,ymax=ymax), fill="gray80") + 
  geom_line(aes(x,y)) + geom_point(data=obs, aes(x,y)) + 
  geom_line(data=df2, aes(x,value,color=variable))
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 


Looks fine for this data, but not for Reed example:




```r
require(pdgControl)
require(ggplot2)
f <- BevHolt
p <- c(1.2,.02)
K <- (p[1]-1)/p[2] 
z_g <- function(sigma_g) rlnorm(1, 0, sigma_g) #1+(2*runif(1, 0,  1)-1)*sigma_g #
x_grid = seq(0, 1.5 * K, length=100)
T <- 40
sigma_g <- 0.01
x <- numeric(T)
x[1] <- 1
for(t in 1:(T-1))
  x[t+1] = z_g(sigma_g) * f(x[t], h=0, p=p)
plot(x)
```

![plot of chunk libraries](figure/libraries.png) 




Predict the function over the target grid


```r
obs <- data.frame(x=x[1:(T-1)],y=x[2:T])
X <- x_grid
system.time(gp <- gp_fit_cholsky(obs, X, c(sigma_n=0.8, l=1)))
```

```
   user  system elapsed 
  0.004   0.000   0.007 
```

```r
system.time(gp_d <- gp_fit_direct(obs, X, c(sigma_n=0.8, l=1)))
```

```
   user  system elapsed 
  0.008   0.000   0.006 
```


Plot the results


```r
require(reshape2)
df <- data.frame(x=X, y=gp$mu, ymin=(gp$mu-2*sqrt(gp$var)), ymax=(gp$mu+2*sqrt(gp$var)))
```

```
Warning: NaNs produced
```

```
Warning: NaNs produced
```

```r
df2 <- data.frame(x=X, y=gp_d$mu, ymin=(gp_d$mu-2*sqrt(diag(gp_d$var))), ymax=(gp_d$mu+2*sqrt(diag(gp_d$var))))
df2 <- melt(df2, id="x")
require(ggplot2)
ggplot(df)  + geom_ribbon(aes(x,y,ymin=ymin,ymax=ymax), fill="gray80") + 
  geom_line(aes(x,y)) + geom_point(data=obs, aes(x,y)) + 
  geom_line(data=df2, aes(x,value,color=variable))
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 



