








```r
X <- seq(-5,5,len=50)
obs <- data.frame(x = c(-4, -3, -1,  0,  2, 3),
                  y = c(-2,  0,  1,  2, -1, -1))
l <- 1
sigma_n <- 0.1
```






### Cholesky method
  

```r
  SE <- function(Xi,Xj, l=1) exp(-0.5 * (Xi - Xj) ^ 2 / l ^ 2)
  cov <- function(X, Y) outer(X, Y, SE, l) 
  n <- length(obs$x)
  K <- cov(obs$x, obs$x)
  I <-  diag(1, n)
  L <- chol(K + sigma_n^2 * I)
  alpha <- forwardsolve(t(L), backsolve(L, obs$y)) # see http://stat.ethz.ch/R-manual/R-patched/library/base/html/chol2inv.html
  k_star <- cov(obs$x, X)
  Y <- t(k_star) %*% alpha
  v <- backsolve(L, k_star)
  Var <- cov(X,X) - t(v) %*% v
  loglik <- -.5 * t(obs$y) %*% alpha - sum(log(diag(L))) - n * log(2 * pi) / 2
```

  
### Direct method 


```r
  d <- gp_fit(obs, X, c(sigma_n=sigma_n, l=l, tau=1), "direct")
  Ef <- d$Ef
  Cf <- d$Cf
```


### Direct sequential method, avoids matrix inverse instability


```r
s <- gp_fit(obs, X, c(sigma_n=sigma_n, l=l, tau=1), "sequential")
ef <- s$Ef
cf <- s$Cf
```




```r
ggplot(data.frame(x=X, Ef=Ef, ef=ef))+ geom_point(aes(x,Ef), col='red') + geom_line(aes(x,ef))
```

![plot of chunk unnamed-chunk-7](http://carlboettiger.info/assets/figures/2012-12-07-c431834819-unnamed-chunk-7.png) 



### kernlab method


```r
k <- gp_fit(obs, X, c(sigma_n=sigma_n, l=l, tau=1), "kernlab")
```

```
Using automatic sigma estimation (sigest) for RBF or laplace kernel 
```

```r
y_p <- k$Ef
```


Compare these results: 


```r
require(reshape2)
df <- data.frame(x = X, direct = Ef, Cholesky = Y, kernlab = y_p, sequential = ef)
df <- melt(df, id = "x")
ggplot(df)+ geom_jitter(aes(x, value, color = variable)) + geom_point(data = obs, aes(x,y))
```

![plot of chunk unnamed-chunk-9](http://carlboettiger.info/assets/figures/2012-12-07-c431834819-unnamed-chunk-9.png) 


