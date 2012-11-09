Quick trial SDP approaches with GP function
==============================================


```r
require(pdgControl)
```

```
## Loading required package: pdgControl
```

```r
require(nonparametricbayes)
```

```
## Loading required package: nonparametricbayes
```

```r
require(ggplot2)
```

```
## Loading required package: ggplot2
```

```r
opts_knit$set(upload.fun = socialR::notebook.url)
```



### Beverton-Holt function

Simulate some training data under a stochastic growth function with standard parameterization,



```r
f <- BevHolt
p <- c(1.5, 0.05)
K <- (p[1] - 1)/p[2]
```




Parameter definitions


```r
x_grid = seq(0, 1.5 * K, length = 101)
T <- 40
sigma_g <- 0.1
x <- numeric(T)
x[1] <- 1
```


Noise function, profit function

```r
z_g <- function() rlnorm(1, 0, sigma_g)  #1+(2*runif(1, 0,  1)-1)*sigma_g #
profit <- profit_harvest(1, 0, 0)
```



Simulation 


```r
for (t in 1:(T - 1)) x[t + 1] = z_g() * f(x[t], h = 0, p = p)
```




Predict the function over the target grid


```r
obs <- data.frame(x = x[1:(T - 1)], y = x[2:T])
X <- x_grid
library(nonparametricbayes)
gp <- gp_fit(obs, X, c(sigma_n = 0.5, l = 1.5))
```


Gaussian Process inference from this model.  True model shown in red.  


```r
df <- data.frame(x = X, y = gp$Ef, ymin = (gp$Ef - 2 * sqrt(abs(diag(gp$Cf)))), 
    ymax = (gp$Ef + 2 * sqrt(abs(diag(gp$Cf)))))
true <- data.frame(x = X, y = sapply(X, f, 0, p))
require(ggplot2)
ggplot(df) + geom_ribbon(aes(x, y, ymin = ymin, ymax = ymax), fill = "gray80") + 
    geom_line(aes(x, y)) + geom_point(data = obs, aes(x, y)) + geom_line(data = true, 
    aes(x, y), col = "red", lty = 2)
```

![plot of chunk unnamed-chunk-5](http://carlboettiger.info/assets/figures/2012-11-09-f312095efb-unnamed-chunk-5.png) 




## Optimization over hyperparameters:



```r
minusloglik <- function(par) {
    gp <- gp_fit(obs, X, par)
    -gp$llik
}
par <- c(sigma_n = 1, l = 1)
minusloglik(par)
```

```
## [1] 2108
```

```r
o <- optim(par, minusloglik)
o
```

```
## $par
## sigma_n       l 
##   8.394  -2.314 
## 
## $value
## [1] 1500
## 
## $counts
## function gradient 
##       49       NA 
## 
## $convergence
## [1] 0
## 
## $message
## NULL
```

```r
hyperpars <- o$par
```


**Yikes.**  

let's try 1-D optimization:



```r
minusloglik <- function(par) {
    gp <- gp_fit(obs, X, c(sigma_n = par, l = 1))
    -gp$llik
}
minusloglik(1)
```

```
## [1] 2108
```

```r
o <- optimize(minusloglik, c(0, 150))
o
```

```
## $minimum
## [1] 8.394
## 
## $objective
## [1] 1500
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

![plot of chunk unnamed-chunk-8](http://carlboettiger.info/assets/figures/2012-11-09-f312095efb-unnamed-chunk-8.png) 






### Misc ways to calc harvest matrices


```r
# escapement <- Vectorize(function(i,j) pmax(0, x_grid[i] -
# h_grid[D[j,(Tmax-t+1)]])) S <- outer(1:n_x, 1:n_h, escapement)
```


Old-school calculation method:



```r
profit <- profit_harvest(price = 1, c0 = 0, c1 = 0)
pdfn <- function(P, s) dlnorm(P, 0, s)
sdp_mat <- determine_SDP_matrix(f, p, x_grid, h_grid, sigma_g, pdfn)
```

```
## Error: object 'h_grid' not found
```

```r
opt <- find_dp_optim(sdp_mat, x_grid = x_grid, h_grid = h_grid, OptTime = 20, 
    xT = 0, profit = profit, delta = 0.05)
```

```
## Error: object 'h_grid' not found
```

```r
plot(opt$D[, 1])
```

```
## Error: object 'opt' not found
```





Somewhat silly way to adjust the matrix by harvest level:


```r
n <- length(h_grid)
```

```
## Error: object 'h_grid' not found
```

```r
silly <- function(F, n) lapply(0:(n - 1), function(i) {
    d <- diag(1, nrow = (n - i))
    top <- matrix(0, nrow = i, ncol = n)
    top[, 1] <- 1
    side <- matrix(0, nrow = (n - i), ncol = i)
    out <- F %*% rbind(top, cbind(d, side))
    t(apply(out, 1, function(x) x/sum(x)))  # normalize
})
```




binning function, not quite right: 


```r
bin <- function(x, grid_centers) {
    match <- which.min(abs(grid_centers - x))
    out <- integer(length(grid_centers))
    out[match] <- 1
    out
}
per_harvest <- function(F, n) {
    lapply(h_grid, function(h) {
        t(sapply(x_grid, function(x) bin(x - h, x_grid)))
    })
}
```

