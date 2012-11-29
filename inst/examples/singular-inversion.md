testing singular matrix condition



```r
require(ggplot2)
```

```
## Loading required package: ggplot2
```

```r
require(kernlab)
```

```
## Loading required package: kernlab
```

```r
require(nonparametricbayes)
```

```
## Loading required package: nonparametricbayes
```

```r
opts_knit$set(upload.fun = socialR::notebook.url)
opts_chunk$set(dev.args = list(bg = "transparent"), comment = NA, 
    tidy = FALSE)
theme_set(theme_bw())
theme_update(panel.background = element_rect(fill = "transparent", 
    colour = NA))
```


Parameterization-specific



```r
X <- seq(0,1, length=1000)
obs <- data.frame(x = X,
                  y = sin(2*pi*10*X)+.3*rnorm(length(X)))
```



  
### Direct method 

Fit with reasonable choice of length scale `l`


```r
d <- gp_fit(obs, X, c(sigma_n=.5, l=.02, tau=.3), "direct")
det(d$K + .5^2 * diag(1, length(X)))
```

```
[1] 0
```

```r
plot.gpfit(d)
```

![plot of chunk unnamed-chunk-3](http://carlboettiger.info/assets/figures/2012-11-28-7faf6791cb-unnamed-chunk-3.png) 


Fit with long length scales.  


```r
d <- gp_fit(obs, X, c(sigma_n=.5, l=10, tau=.3), "direct")
det(d$K + .5^2 * diag(1, length(X)))
```

```
[1] 0
```

```r
plot.gpfit(d)
```

![plot of chunk unnamed-chunk-4](http://carlboettiger.info/assets/figures/2012-11-28-7faf6791cb-unnamed-chunk-4.png) 

