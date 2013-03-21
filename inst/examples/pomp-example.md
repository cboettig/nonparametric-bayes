


```r
library(pomp)
```

```
## Loading required package: mvtnorm
```

```
## Loading required package: subplex
```

```
## Loading required package: deSolve
```

```r
opts_knit$set(upload.fun = socialR::flickr.url)
```

```r
process.sim <- function(x, t, params, delta.t, ...) {
    # r = 0.75, K = 10, Q = 3, H = 1, a = 1.3
    r <- params["r"]
    K <- params["K"]
    Q <- params["Q"]
    H <- params["H"]
    a <- params["a"]
    sigma <- params["sigma"]
    X <- x["X"]
    epsilon <- exp(rlnorm(n = 1, mean = 0, sd = sigma))
    
    c(X = unname(X * exp(r * delta.t * (1 - X/K) - a * X^(Q - 1)/(X^Q + H^Q)) * 
        epsilon))
}
```



```r
measure.sim <- function(x, t, params, ...) {
    tau <- params["tau"]
    X <- x["X"]
    y <- c(Y = unname(rlnorm(n = 1, meanlog = log(X), sd = tau)))
}
```


Create a container of class `pomp`


```r
may <- pomp(data = data.frame(time = 1:100, Y = NA), times = "time", 
    rprocess = discrete.time.sim(step.fun = process.sim, delta.t = 1), rmeasure = measure.sim, 
    t0 = 0)
```



Parameters and an initial condition






```r
theta <- c(r = 0.75, K = 10, Q = 3, H = 1, a = 1.3, sigma = 0.1, 
    tau = 0.1, X.0 = 5)
```

```

Simulate


```r
dat <- simulate(may, params = theta, obs = TRUE, states = TRUE)
library(reshape2)
library(ggplot2)
dat <- melt(dat)[3:5]
names(dat) = c("time", "value", "variable")
ggplot(dat, aes(time, value, lty = variable)) + geom_line()
```

![plot of chunk unnamed-chunk-6](http://farm9.staticflickr.com/8252/8576510722_aa95db3a10_o.png) 


Note that simulate returns an instance of our object, with the data column filled in.  This is not true if we give it arguments such as state or obs = TRUE.  

Measurement density


```r
measure.dens <- function(y, x, t, params, log, ...) {
    tau <- params["tau"]
    ## state at time t:
    X <- x["X"]
    ## observation at time t:
    Y <- y["Y"]
    ## compute the likelihood of Y|X,tau
    dlnorm(x = Y, meanlog = log(X), sdlog = tau, log = log)
}
```


Add to our container


```r
may <- pomp(may, dmeasure = measure.dens)
```





Optionally, work in the transformed variable space instead by providing parameter and inverse parameter transforms.  


```r
may <- pomp(may, parameter.transform = function(params, ...) {
    exp(params)
}, parameter.inv.transform = function(params, ...) {
    log(params)
})
```


Run the particle filter


```r
pf <- pfilter(may, params = theta, Np = 1000)
```

```
## Error: 'pfilter' error: 'dmeasure' returns non-finite value
```

```r
logLik(pf)
```

```
## Error: no applicable method for 'logLik' applied to an object of class
## "function"
```






