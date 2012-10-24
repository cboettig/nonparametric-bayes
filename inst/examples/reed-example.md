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

![plot of chunk unnamed-chunk-3](http://farm9.staticflickr.com/8190/8117649167_f5002b7024_o.png) 




Predict the function over the target grid


```r
obs <- data.frame(x = x[1:(T - 1)], y = x[2:T])
X <- x_grid
library(nonparametricbayes)
gp <- gp_fit(obs, X, c(sigma_n = 0.05, l = 10))
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

![plot of chunk unnamed-chunk-5](http://farm9.staticflickr.com/8336/8117660178_beeca82760_o.png) 









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

![plot of chunk unnamed-chunk-7](http://farm9.staticflickr.com/8046/8117660240_31c8d052c4_o.png) 



Simulation 


```r
x[1] = 2.5
for (t in 1:(T - 1)) x[t + 1] = z_g(sigma_g) * f(x[t], h = 0, p = p)
plot(x)
```

![plot of chunk unnamed-chunk-8](http://farm9.staticflickr.com/8330/8117649515_59c703d880_o.png) 


Predict the function over the target grid


```r
obs <- data.frame(x = x[1:(T - 1)], y = x[2:T])
X <- x_grid
library(nonparametricbayes)
gp <- gp_fit(obs, X, c(sigma_n = 0.05, l = 10))
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

![plot of chunk unnamed-chunk-10](http://farm9.staticflickr.com/8325/8117660460_a53eb7c10c_o.png) 


