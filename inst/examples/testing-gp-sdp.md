Quick trial SDP approaches with GP function
==============================================



```r
knit("beverton_holt_data.Rmd")
```

```
## 
## 
## processing file: beverton_holt_data.Rmd
```

```
## output file:
## /home/cboettig/Documents/code/nonparametric-bayes/inst/examples/beverton_holt_data.md
```

```r
opts_knit$set(upload.fun = socialR::notebook.url)
```



Radial basis function/Gaussian kernel:


```r
SE <- function(Xi, Xj, l) exp(-0.5 * (Xi - Xj)^2/l^2)
cov <- function(X, Y) outer(X, Y, SE, l)
```



Cholesky method
  

```r
n <- length(obs$x)
K <- cov(obs$x, obs$x)
I <- diag(1, n)

L <- chol(K + sigma_n^2 * I)
alpha <- solve(t(L), solve(L, obs$y))
k_star <- cov(obs$x, X)
Y <- t(k_star) %*% alpha
v <- solve(L, k_star)
Var <- cov(X, X) - t(v) %*% v
loglik <- -0.5 * t(obs$y) %*% alpha - sum(log(diag(L))) - n * log(2 * 
    pi)/2
```

  
Direct method 


```r
cov_xx_inv <- solve(K + sigma_n^2 * I)
Ef <- cov(X, obs$x) %*% cov_xx_inv %*% obs$y
Cf <- cov(X, X) - cov(X, obs$x) %*% cov_xx_inv %*% cov(obs$x, X)
```


  Direct sequential method, avoids matrix inverse instability

```r
ef <- numeric(length(X))
cf <- matrix(0, nrow = length(X), ncol = length(X))
for (i in 1:length(obs$x)) {
    ef <- ef + cov(X, obs$x[i]) * obs$y[i]/as.numeric(cov(obs$x[i], obs$x[i]) + 
        sigma_n^2)
    cf <- cf + cov(X, X) - cov(X, obs$x[i]) %*% cov(obs$x[i], X)/as.numeric(cov(obs$x[i], 
        obs$x[i]) + sigma_n^2)
}
```




```r
ggplot(data.frame(x = X, Ef = Ef, ef = ef)) + geom_point(aes(x, Ef), 
    col = "red") + geom_line(aes(x, ef))
```

![plot of chunk unnamed-chunk-6](http://carlboettiger.info/assets/figures/2012-11-09-b2bfbd5a18-unnamed-chunk-6.png) 



kernlab method


```r
gp <- gausspr(obs$x, obs$y, kernel = "rbfdot", kpar = list(sigma = 1/(2 * 
    l^2)), fit = FALSE, scaled = FALSE, var = 0.8)
y_p <- predict(gp, X)
```



mlegp method


```r
require(mlegp)
out <- mlegp(obs$x, obs$y, nugget = 1)
```

```
## 
## ========== FITTING GP # 1 ==============================
## intial_scaled nugget is 0.111782
## running simplex # 1...
## ...done
## ...simplex #1 complete, loglike = -60.566054 (convergence)
## running simplex # 2...
## ...done
## ...simplex #2 complete, loglike = -60.566055 (convergence)
## running simplex # 3...
## ...done
## ...simplex #3 complete, loglike = -60.566053 (convergence)
## running simplex # 4...
## ...done
## ...simplex #4 complete, loglike = -60.566053 (convergence)
## running simplex # 5...
## ...done
## ...simplex #5 complete, loglike = -60.566053 (convergence)
## 
## using L-BFGS method from simplex #5...
## 	iteration: 1,loglike = -60.566053
## ...L-BFGS method complete
## 
## Maximum likelihood estimates found, log like =  -60.566053
## addNuggets...
## creating gp object......done
```

```r
y_m <- predict(out, as.matrix(X))
```


Compare these results: 


```r
require(reshape2)
df <- data.frame(x = X, direct = Ef, Cholesky = Y, kernlab = y_p, 
    sequential = ef, mlegp = y_m)
df <- melt(df, id = "x")
ggplot(df) + geom_jitter(aes(x, value, color = variable)) + geom_point(data = obs, 
    aes(x, y))
```

![plot of chunk unnamed-chunk-9](http://carlboettiger.info/assets/figures/2012-11-09-b2bfbd5a18-unnamed-chunk-9.png) 



