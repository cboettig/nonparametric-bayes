




```r
require(MASS)
require(ggplot2)
require(kernlab)
opts_knit$set(upload.fun = socialR::notebook.url)
```


Parameterization-specific



```r
X <- seq(-5,5,len=50)
obs <- data.frame(x = c(-4, -3, -1,  0,  2, 3),
                  y = c(-2,  0,  1,  2, -1, -1))
#x <- -5:5
#obs <- data.frame(x = -5:5, y = sin(x) + rnorm(length(x),sd=.1))
l <- 1
sigma_n <- 0.8
```



Radial basis function/Gaussian kernel:


```r
  SE <- function(Xi,Xj, l) exp(-0.5 * (Xi - Xj) ^ 2 / l ^ 2)
  cov <- function(X, Y) outer(X, Y, SE, l) 
```



Cholesky method
  

```r
  n <- length(obs$x)
  K <- cov(obs$x, obs$x)
  I <-  diag(1, n)

  L <- chol(K + sigma_n^2 * I)
  alpha <- solve(t(L), solve(L, obs$y))
  k_star <- cov(obs$x, X)
  Y <- t(k_star) %*% alpha
  v <- solve(L, k_star)
  Var <- cov(X,X) - t(v) %*% v
  loglik <- -.5 * t(obs$y) %*% alpha - sum(log(diag(L))) - n * log(2 * pi) / 2
```

  
Direct method 


```r
  cov_xx_inv <- solve(K + sigma_n^2*I)
  Ef <- cov(X, obs$x) %*% cov_xx_inv %*% obs$y
  Cf <- cov(X, X) - cov(X, obs$x)  %*% cov_xx_inv %*% cov(obs$x, X)
```


  Direct sequential method, avoids matrix inverse instability

```r
ef <- numeric(length(X))
cf <- matrix(0, nrow=length(X), ncol=length(X))
for(i in 1:length(obs$x)){
  ef <- ef + cov(X, obs$x[i]) * obs$y[i] / as.numeric(cov(obs$x[i],obs$x[i]) + sigma_n^2)
  cf <- cf + cov(X, X) - cov(X, obs$x[i])  %*% cov(obs$x[i], X) / as.numeric(cov(obs$x[i],obs$x[i]) + sigma_n^2)
}
```




```r
ggplot(data.frame(x=X, Ef=Ef, ef=ef))+ geom_point(aes(x,Ef), col='red') + geom_line(aes(x,ef))
```

![plot of chunk unnamed-chunk-7](http://carlboettiger.info/assets/figures/2012-11-09-d5cb8c43b1-unnamed-chunk-7.png) 



kernlab method


```r
gp <- gausspr(obs$x, obs$y, kernel="rbfdot", kpar=list(sigma=1/(2*l^2)), fit=FALSE, scaled=FALSE, var=0.8)
y_p <- predict(gp, X)
```



mlegp method


```r
require(mlegp)
out <- mlegp(obs$x, obs$y, nugget = 1)
```

```

========== FITTING GP # 1 ==============================
intial_scaled nugget is 0.461538
running simplex # 1...
...done
...simplex #1 complete, loglike = -9.812179 (convergence)
running simplex # 2...
...done
...simplex #2 complete, loglike = -9.812189 (convergence)
running simplex # 3...
...done
...simplex #3 complete, loglike = -9.812179 (convergence)
running simplex # 4...
...done
...simplex #4 complete, loglike = -9.812203 (convergence)
running simplex # 5...
...done
...simplex #5 complete, loglike = -9.812179 (convergence)

using L-BFGS method from simplex #1...
...L-BFGS method complete

Maximum likelihood estimates found, log like =  -9.812179
addNuggets...
creating gp object......done

```

```r
y_m <- predict(out, as.matrix(X))
```


Compare these results: 


```r
require(reshape2)
df <- data.frame(x = X, Ef = Ef, Y = Y, y_p = y_p, ef = ef, y_m = y_m)
df <- melt(df, id = "x")
ggplot(df)+ geom_jitter(aes(x, value, color = variable)) + geom_point(data = obs, aes(x,y))
```

![plot of chunk unnamed-chunk-10](http://carlboettiger.info/assets/figures/2012-11-09-d5cb8c43b1-unnamed-chunk-10.png) 




Compare alpha from Cholesky and alpha from `kernlab` to $K(x,x)^{-1} \vec y$


```r
cov_xx_inv %*% obs$y
```

```
        [,1]
[1,] -1.4059
[2,]  0.5008
[3,]  0.1319
[4,]  1.2165
[5,] -0.5589
[6,] -0.4113
```

```r
alpha
```

```
[1] -1.21284  0.46035  0.08053  1.42053 -0.49258 -0.52138
```

```r
alpha(gp)
```

```
        [,1]
[1,] -1.2475
[2,]  0.4010
[3,]  0.1681
[4,]  1.0927
[5,] -0.5068
[6,] -0.3916
```


