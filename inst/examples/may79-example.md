GP Example using the May (1979) bistable model
==============================================



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
# opts_knit$set(upload.fun = socialR::flickr.url)
```







```r
f <- May
p <- c(r = 0.75, k = 10, a = 1.3, H = 1, Q = 3)
K <- 8  # approx
```


Model dynamics look like this:


```r
birth <- function(x) p["r"] * (1 - x/p["k"])
death <- function(x) p["a"] * x^(p["Q"] - 1)/(x^p["Q"] + p["H"])
df <- data.frame(x = x_grid, b = sapply(x_grid, birth), d = sapply(x_grid, 
    death))
```

```
## Error: object 'x_grid' not found
```

```r
ggplot(df) + geom_line(aes(x, b), col = "blue") + geom_line(aes(x, 
    d), col = "red")
```

```
## Error: ggplot2 doesn't know how to deal with data of class function
```




Noise function, profit function

```r
z_g <- function(sigma_g) rlnorm(1, 0, sigma_g)  #1+(2*runif(1, 0,  1)-1)*sigma_g #
profit <- profit_harvest(1, 0, 0)
```



Parameter definitions


```r
x_grid = seq(0, 1.5 * K, length = 101)
T <- 40
sigma_g <- 0.1
x <- numeric(T)
x[1] <- 1
```




Simulation 


```r
x[1] = 2.5
for (t in 1:(T - 1)) x[t + 1] = z_g(sigma_g) * f(x[t], h = 0, p = p)
plot(x)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 


Predict the function over the target grid


```r
obs <- data.frame(x = x[1:(T - 1)], y = x[2:T])
X <- x_grid
gp <- gp_fit(obs, X, c(sigma_n = 1, l = 1))
```

```
## Error: could not find function "gp_fit"
```


Gaussian Process inference from this model


```r
df <- data.frame(x = X, y = gp$Ef, ymin = (gp$Ef - 2 * sqrt(abs(diag(gp$Cf)))), 
    ymax = (gp$Ef + 2 * sqrt(abs(diag(gp$Cf)))))
```

```
## Error: object 'gp' not found
```

```r
true <- data.frame(x = X, y = sapply(X, f, 0, p))
ggplot(df) + geom_ribbon(aes(x, y, ymin = ymin, ymax = ymax), fill = "gray80") + 
    geom_line(aes(x, y)) + geom_point(data = obs, aes(x, y)) + geom_line(data = true, 
    aes(x, y), col = "red", lty = 2)
```

```
## Error: ggplot2 doesn't know how to deal with data of class function
```

