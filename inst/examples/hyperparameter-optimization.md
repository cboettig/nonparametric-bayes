

## Optimization over hyperparameters:


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
minusloglik <- function(par) {
    gp <- gp_fit(obs, X, par)
    -gp$llik
}
par <- c(sigma_n = 1, l = 1)
minusloglik(par)
```

```
## [1] 2158
```

```r
o <- optim(par, minusloglik)
o
```

```
## $par
## sigma_n       l 
##   8.687  -2.442 
## 
## $value
## [1] 1502
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
## [1] 2158
```

```r
o <- optimize(minusloglik, c(0, 150))
o
```

```
## $minimum
## [1] 8.692
## 
## $objective
## [1] 1502
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

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 




