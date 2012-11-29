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


```r
d <- gp_fit(obs, X, c(sigma_n=1, l=.0001, tau=.3), "direct")
plot.gpfit(d)
```

![plot of chunk unnamed-chunk-3](http://carlboettiger.info/assets/figures/2012-11-28-113df7a98a-unnamed-chunk-3.png) 




```r
d <- gp_fit(obs, X, c(sigma_n=1, l=10, tau=.3), "direct")
plot.gpfit(d)
```

![plot of chunk unnamed-chunk-4](http://carlboettiger.info/assets/figures/2012-11-28-113df7a98a-unnamed-chunk-4.png) 

