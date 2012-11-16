




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
sigma_n <- 0.
```



Radial basis function/Gaussian kernel:


```r
  SE <- function(Xi,Xj, l=1) exp(-0.5 * (Xi - Xj) ^ 2 / l ^ 2)
  cov <- function(X, Y) outer(X, Y, SE, l) 
```



### Cholesky method
  

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

  
### Direct method 


```r
  cov_xx_inv <- solve(K + sigma_n^2*I)
  Ef <- cov(X, obs$x) %*% cov_xx_inv %*% obs$y
  Cf <- cov(X, X) - cov(X, obs$x)  %*% cov_xx_inv %*% cov(obs$x, X)
```


### Direct sequential method, avoids matrix inverse instability


```r
ef <- numeric(length(X))
cf <- matrix(0, nrow=length(X), ncol=length(X))

A <- as.numeric( cov(obs$x[1], obs$x[1]) )
mu <- obs$y[1] 

for(i in 2:length(obs$x)){
  mu <- obs$y[i] + cov(obs$x[i], obs$x[1:(i-1)]) %*% (obs$y[1:(i-1)] - mu) / (A + sigma_n^2)
  A <- as.numeric( cov(obs$x[i], obs$x[i]) - cov(obs$x[i], obs$x[1:(i-1)]) %*% cov(obs$x[1:(i-1)], obs$x[i]) / (A + sigma_n^2) )
}
  
ef <- cov(X, obs$x) %*% (obs$y - mu) / A
cf <- cov(X,X) - cov(X, obs$x) %*% cov(obs$x, X) / A
```




```r
ggplot(data.frame(x=X, Ef=Ef, ef=ef))+ geom_point(aes(x,Ef), col='red') + geom_line(aes(x,ef))
```

![plot of chunk unnamed-chunk-7](http://carlboettiger.info/assets/figures/2012-11-15-e732f74985-unnamed-chunk-7.png) 



### kernlab method


```r
gp <- gausspr(obs$x, obs$y, kernel="rbfdot", kpar=list(sigma=1/(2*l^2)), fit=FALSE, scaled=FALSE, var=0.8)
y_p <- predict(gp, X)
```



These are the same as `y_p`, explicitly showing the $E_y = K(x, x_i) \alpha$ multiplication


```r
y_mult <-  kernelMult(kernelf(gp), X, xmatrix(gp), as.matrix(alpha(gp))) 
y_se <- sapply(X, function(x) cov(x, obs$x) %*% alpha(gp))
rbind(as.numeric(y_p), as.numeric(y_mult), as.numeric(y_se))[,1:10]
```

```
        [,1]    [,2]    [,3]   [,4]   [,5]    [,6]    [,7]    [,8]    [,9]
[1,] -0.7023 -0.8287 -0.9338 -1.003 -1.026 -0.9967 -0.9151 -0.7894 -0.6324
[2,] -0.7023 -0.8287 -0.9338 -1.003 -1.026 -0.9967 -0.9151 -0.7894 -0.6324
[3,] -0.7023 -0.8287 -0.9338 -1.003 -1.026 -0.9967 -0.9151 -0.7894 -0.6324
       [,10]
[1,] -0.4598
[2,] -0.4598
[3,] -0.4598
```


Which demonstrates the kernel trick `alpha = cov_xx_inv %*% obs$y`

### mlegp method


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
df <- data.frame(x = X, direct = Ef, Cholesky = Y, kernlab = y_p, sequential = ef)
df <- melt(df, id = "x")
ggplot(df)+ geom_jitter(aes(x, value, color = variable)) + geom_point(data = obs, aes(x,y))
```

![plot of chunk unnamed-chunk-11](http://carlboettiger.info/assets/figures/2012-11-15-e732f74985-unnamed-chunk-11.png) 




Compare alpha from Cholesky and alpha from `kernlab` to $K(x,x)^{-1} \vec y$


```r
cov_xx_inv %*% obs$y
```

```
        [,1]
[1,] -3.2303
[2,]  2.0421
[3,] -0.8287
[4,]  2.6398
[5,] -1.1452
[6,] -0.3344
```

```r
alpha
```

```
[1] -2.0407  1.6531 -0.8527  3.9213 -0.8673 -0.9967
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


