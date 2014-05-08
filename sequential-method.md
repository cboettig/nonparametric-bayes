---
layout: page
---

## Sequential method


```r
require(MASS)
require(ggplot2)
require(kernlab)
opts_knit$set(upload.fun = socialR::notebook.url)
```


Parameterization-specific



```r
X <- seq(-5,5, length=40)
obs <- data.frame(x = c(-3, -1, 3),
                  y = c(0,  1, -1))
#obs <- data.frame(x = -5:5, y = sin(x) + rnorm(length(x),sd=.1))
l <- 1
sigma_n <- 0.5
```



Radial basis function/Gaussian kernel:


```r
  SE <- function(Xi,Xj, l=1) exp(-0.5 * (Xi - Xj) ^ 2 / l ^ 2)
  cov <- function(X, Y) outer(X, Y, SE, l) 
```

  
### Direct method 


```r
  n <- length(obs$x)
  K <- cov(obs$x, obs$x)
  I <-  diag(1, n)
  
  cov_xx_inv <- solve(K + sigma_n^2*I)
  Ef <- cov(X, obs$x) %*% cov_xx_inv %*% obs$y
  Cf <- cov(X, X) - cov(X, obs$x)  %*% cov_xx_inv %*% cov(obs$x, X)
```


### Direct sequential method, avoids matrix inverse instability


```r

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
ef <- t(mu_seq(X, length(obs$x)))
cf <- C_seq(X, X, length(obs$x))

```




Compare these results: 


```r
require(reshape2)
require(ggplot2)
df <- data.frame(x = X, direct = Ef, sequential = ef)
df <- melt(df, id = "x")
ggplot(df)+ geom_jitter(aes(x, value, color = variable)) + geom_point(data = obs, aes(x,y))
```

![plot of chunk unnamed-chunk-6](http://carlboettiger.info/assets/figures/2012-11-19-58073767b5-unnamed-chunk-6.png) 






