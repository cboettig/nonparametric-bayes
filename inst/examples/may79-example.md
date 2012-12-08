GP Example using the May (1979) bistable model
==============================================


Use my custom libraries and the `ggplot` plotting tools


```r
require(nonparametric - bayes)
```

```
## Warning: the condition has length > 1 and only the first element will be
## used
```

```
## c("Loading required package: -", "Loading required package:
## nonparametric", "Loading required package: bayes")
```

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







