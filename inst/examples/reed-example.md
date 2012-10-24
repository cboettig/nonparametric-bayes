Reed (1979) with Gaussian processes 
========================================================


```r
require(pdgControl)
```

```
## Loading required package: pdgControl
```

```r
require(ggplot2)
```

```
## Loading required package: ggplot2
```

```r
opts_knit$set(upload.fun = socialR::flickr.url)
```


Simulate some training data under a stochastic growth function with standard parameterization,



```r
f <- BevHolt
p <- c(1.5, 0.05)
K <- (p[1] - 1)/p[2]
```



Noise function 

```r
z_g <- function(sigma_g) rlnorm(1, 0, sigma_g)  #1+(2*runif(1, 0,  1)-1)*sigma_g #
```



Parameter definitions


```r
x_grid = seq(0, 1.5 * K, length = 100)
T <- 40
sigma_g <- 0.1
x <- numeric(T)
x[1] <- 1
```


Simulation 


```r
for (t in 1:(T - 1)) x[t + 1] = z_g(sigma_g) * f(x[t], h = 0, p = p)

plot(x)
```

![plot of chunk unnamed-chunk-3](http://farm9.staticflickr.com/8465/8117743934_a153fbabd3_o.png) 




Predict the function over the target grid


```r
obs <- data.frame(x = x[1:(T - 1)], y = x[2:T])
X <- x_grid
library(nonparametricbayes)
gp <- gp_fit(obs, X, c(sigma_n = 1, l = 1))
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

![plot of chunk unnamed-chunk-5](http://farm9.staticflickr.com/8333/8117733641_c30742893a_o.png) 









## Another example using the May model


```r
f <- May
p <- c(r = 0.75, k = 10, a = 1.2, H = 1, Q = 3)
K <- 8  # approx
```


Model dynamics look like this:


```r
birth <- function(x) p["r"] * (1 - x/p["k"])
death <- function(x) p["a"] * x^(p["Q"] - 1)/(x^p["Q"] + p["H"])
df <- data.frame(x = x_grid, b = sapply(x_grid, birth), d = sapply(x_grid, 
    death))
ggplot(df) + geom_line(aes(x, b), col = "blue") + geom_line(aes(x, 
    d), col = "red")
```

![plot of chunk unnamed-chunk-7](http://farm9.staticflickr.com/8052/8117733757_0a75a31cc7_o.png) 



Simulation 


```r
x[1] = 2.5
for (t in 1:(T - 1)) x[t + 1] = z_g(sigma_g) * f(x[t], h = 0, p = p)
plot(x)
```

![plot of chunk unnamed-chunk-8](http://farm9.staticflickr.com/8184/8117744294_07f74257a5_o.png) 


Predict the function over the target grid


```r
obs <- data.frame(x = x[1:(T - 1)], y = x[2:T])
X <- x_grid
gp <- gp_fit(obs, X, c(sigma_n = 1, l = 1))
```


Gaussian Process inference from this model


```r
df <- data.frame(x = X, y = gp$Ef, ymin = (gp$Ef - 2 * sqrt(abs(diag(gp$Cf)))), 
    ymax = (gp$Ef + 2 * sqrt(abs(diag(gp$Cf)))))
true <- data.frame(x = X, y = sapply(X, f, 0, p))
ggplot(df) + geom_ribbon(aes(x, y, ymin = ymin, ymax = ymax), fill = "gray80") + 
    geom_line(aes(x, y)) + geom_point(data = obs, aes(x, y)) + geom_line(data = true, 
    aes(x, y), col = "red", lty = 2)
```

![plot of chunk unnamed-chunk-10](http://farm9.staticflickr.com/8184/8117744504_dc8a6f906d_o.png) 



## Simple optimization over hyperparameters:



```r
minusloglik <- function(par) {
    gp <- gp_fit(obs, X, par)
    -gp$llik
}
par <- c(sigma_n = 1, l = 1)
minusloglik(par)
```

```
##       [,1]
## [1,] 120.6
```

```r
o <- optim(par, minusloglik)
o
```

```
## $par
## sigma_n       l 
##   14023    2613 
## 
## $value
## [1] -336.4
## 
## $counts
## function gradient 
##       67       NA 
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


Yikes.  let's try 1-D optimization:



```r
minusloglik <- function(par) {
    gp <- gp_fit(obs, X, c(sigma_n = par, l = 1))
    -gp$llik
}
minusloglik(1)
```

```
##       [,1]
## [1,] 120.6
```

```r
o <- optimize(minusloglik, c(0, 150))
o
```

```
## $minimum
## [1] 150
## 
## $objective
##        [,1]
## [1,] -159.5
```

```r
hyperpars <- c(sigma_n = o$minimum, l = 1)
```




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

![plot of chunk unnamed-chunk-13](http://farm9.staticflickr.com/8327/8117744938_f7f3cf82e2_o.png) 


