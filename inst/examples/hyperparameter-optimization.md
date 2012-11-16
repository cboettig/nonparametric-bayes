

## Optimization over hyperparameters:


```r
opts_knit$set(upload.fun = socialR::notebook.url)
```




Static definition of gp function


```r
gp_fit <- function(obs, X, pars = c(sigma_n = 1, l = 1)) {
    sigma_n <- pars["sigma_n"]
    l <- pars["l"]
    
    SE <- function(Xi, Xj, l = l) exp(-0.5 * (Xi - Xj)^2/l^2)
    cov <- function(X, Y) outer(X, Y, SE, l)
    
    n <- length(obs$x)
    K <- cov(obs$x, obs$x)
    I <- diag(1, n)
    
    cov_xx_inv <- solve(K + sigma_n^2 * I)
    Ef <- cov(X, obs$x) %*% cov_xx_inv %*% obs$y
    Cf <- cov(X, X) - cov(X, obs$x) %*% cov_xx_inv %*% cov(obs$x, X)
    llik <- -0.5 * t(obs$y) %*% cov_xx_inv %*% obs$y - 0.5 * log(det(cov_xx_inv)) - 
        n * log(2 * pi)/2
    
    out <- list(Ef = Ef, Cf = Cf, llik = llik, obs = obs, X = X, pars = pars, 
        K = K, inv = cov_xx_inv)
    class(out) = "gpfit"
    out
}


plot.gpfit <- function(gp) {
    dat <- data.frame(x = gp$X, y = (gp$Ef), ymin = (gp$Ef - 2 * sqrt(diag(gp$Cf))), 
        ymax = (gp$Ef + 2 * sqrt(diag(gp$Cf))))
    ggplot(dat) + geom_ribbon(aes(x = x, y = y, ymin = ymin, ymax = ymax), fill = "grey80") + 
        geom_line(aes(x = x, y = y), size = 1) + geom_point(data = gp$obs, aes(x = x, 
        y = y)) + labs(title = paste("llik =", prettyNum(gp$llik)))
}
```




Estimate a GP from data with given hyperparameters, then redraw data from an actual multivariate normal.  Presumably optimal hyperparameters might be near the given values then.  


```r
require(MASS)
X <- seq(-5, 5, len = 20)
obs <- data.frame(x = c(-4, -3, -1, 0, 2, 3), y = c(-2, 0, 1, 2, 
    -1, -1))

fit <- gp_fit(obs, X, pars = c(sigma_n = 0.3, l = 0.6))
gp_dat <- mvrnorm(1, fit$Ef, fit$Cf)
obs <- data.frame(x = X, y = gp_dat)
X <- seq(-5, 5, len = 50)
```



Example plot


```r
plot.gpfit(gp_fit(obs, X, pars = c(sigma_n = 0.3, l = 0.6)))
```

![plot of chunk unnamed-chunk-4](http://carlboettiger.info/assets/figures/2012-11-15-28b3256dfe-unnamed-chunk-4.png) 




estimates non-boundary solution of `l`:


```r
minusloglik <- function(par) {
    gp <- gp_fit(obs, X, c(sigma_n = 0.3, l = par))
    -gp$llik
}
minusloglik(1)
```

```
##       [,1]
## [1,] 37.76
```

```r
o <- optimize(minusloglik, c(0, 40))
o
```

```
## $minimum
## [1] 0.3614
## 
## $objective
##       [,1]
## [1,] 26.87
```



But makes `sigma_n` arbitrarily large (boundary solution).  


```r
minusloglik <- function(par) {
    gp <- gp_fit(obs, X, c(sigma_n = par, l = 0.6))
    -gp$llik
}
minusloglik(1)
```

```
##       [,1]
## [1,] 16.88
```

```r
o <- optimize(minusloglik, c(0, 40))
o
```

```
## $minimum
## [1] 40
## 
## $objective
##       [,1]
## [1,] -55.4
```





```r
gp_k <- gausspr(obs$x, obs$y, kernel = "rbfdot")
```

```
## Using automatic sigma estimation (sigest) for RBF or laplace kernel
```

```r
gp_k
```

```
## Gaussian Processes object of class "gausspr" 
## Problem type: regression 
## 
## Gaussian Radial Basis kernel function. 
##  Hyperparameter : sigma =  4.79808910472974 
## 
## Number of training instances learned : 20 
## Train error : 0.131915781
```


Convert to `l`, `sigma=1/(2l^2)`



```r
1/sqrt(2 * as.numeric(gp_k@kernelf@kpar))
```

```
## [1] 0.3228
```


Weirdly inconsistent?


```r
gp_k <- gausspr(obs$x, obs$y, kernel = "rbfdot")
```

```
## Using automatic sigma estimation (sigest) for RBF or laplace kernel
```

```r
1/sqrt(2 * as.numeric(gp_k@kernelf@kpar))
```

```
## [1] 0.1684
```

```r
gp_k <- gausspr(obs$x, obs$y, kernel = "rbfdot")
```

```
## Using automatic sigma estimation (sigest) for RBF or laplace kernel
```

```r
1/sqrt(2 * as.numeric(gp_k@kernelf@kpar))
```

```
## [1] 0.3206
```





### a second example, using the beverton holt data:


```r
knit("beverton_holt_data.Rmd")
```

```
## 
## 
## processing file: beverton_holt_data.Rmd
```

![plot of chunk unnamed-chunk-10](http://carlboettiger.info/assets/figures/2012-11-15-28b3256dfe-unnamed-chunk-10.png) 

```
## output file:
## /home/cboettig/Documents/code/nonparametric-bayes/inst/examples/beverton_holt_data.md
```






```r
minusloglik <- function(par) {
    if (any(par < 0)) 
        Inf else {
        gp <- gp_fit(obs, X, par)
        -gp$llik
    }
}
par <- c(sigma_n = 1, l = 1)
minusloglik(par)
```

```
##       [,1]
## [1,] 184.9
```

```r
# o <- optim(par, minusloglik)
o <- optim(par, minusloglik, method = "L", lower = c(0, 0), upper = c(20, 
    20))
o
```

```
## $par
## sigma_n       l 
##      20      20 
## 
## $value
## [1] -77.64
## 
## $counts
## function gradient 
##        2        2 
## 
## $convergence
## [1] 0
## 
## $message
## [1] "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL"
```

```r
hyperpars <- o$par
```


**Yikes.**  

let's try 1-D optimization:



```r
minusloglik <- function(par) {
    gp <- gp_fit(obs, X, c(sigma_n = par, l = 5.7))
    -gp$llik
}
minusloglik(1)
```

```
##      [,1]
## [1,]  105
```

```r
o <- optimize(minusloglik, c(0, 50))
o
```

```
## $minimum
## [1] 50
## 
## $objective
##        [,1]
## [1,] -116.2
```

```r
hyperpars <- c(sigma_n = o$minimum, l = 1)
```


**Yikes, that shouldn't happen either**


Plot the optimal hyperparameter solution:


```r
gp <- gp_fit(obs, X, hyperpars)
df <- data.frame(x = X, y = gp$Ef, ymin = (gp$Ef - 2 * sqrt(abs(diag(gp$Cf)))), 
    ymax = (gp$Ef + 2 * sqrt(abs(diag(gp$Cf)))))
true <- data.frame(x = X, y = sapply(X, f, 0, p))
ggplot(df) + geom_ribbon(aes(x, y, ymin = ymin, ymax = ymax), fill = "gray80") + 
    geom_line(aes(x, y)) + geom_point(data = obs, aes(x, y)) + geom_line(data = true, 
    aes(x, y), col = "red", lty = 2)
```

![plot of chunk unnamed-chunk-13](http://carlboettiger.info/assets/figures/2012-11-15-28b3256dfe-unnamed-chunk-13.png) 




