Quick trial SDP approaches with GP function
==============================================


```r
require(pdgControl)
require(nonparametricbayes)
require(ggplot2)
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

![plot of chunk unnamed-chunk-5](http://carlboettiger.info/assets/figures/2012-11-09-31878488d5-unnamed-chunk-5.png) 




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
## [1] 2196
```

```r
o <- optim(par, minusloglik)
o
```

```
## $par
## sigma_n       l 
##   8.918  -2.544 
## 
## $value
## [1] 1503
## 
## $counts
## function gradient 
##       53       NA 
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
## [1] 2196
```

```r
o <- optimize(minusloglik, c(0, 150))
o
```

```
## $minimum
## [1] 8.915
## 
## $objective
## [1] 1503
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

![plot of chunk unnamed-chunk-8](http://carlboettiger.info/assets/figures/2012-11-09-31878488d5-unnamed-chunk-8.png) 




