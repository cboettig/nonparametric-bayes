Reed (1979) with Gaussian processes 
========================================================






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


Predict the function over the target grid


```r
library(nonparametricbayes)
gp <- gp_fit(obs, X, c(sigma_n = 1, l = 1.5))
plot.gpfit(gp)
```

![plot of chunk unnamed-chunk-4](http://carlboettiger.info/assets/figures/2012-11-16-e670f9bebf-unnamed-chunk-4.png) 


Solution using `kernlab`


```r
library(kernlab)
# gp <- gausspr(obs$x, obs$y, kernel='rbfdot', kpar=list(sigma=0.5),
# fit=FALSE, scaled=FALSE, var=.09)
gp <- gausspr(obs$x, obs$y, kernel = "rbfdot", scaled = FALSE, var = 0.01)
```

```
## Using automatic sigma estimation (sigest) for RBF or laplace kernel
```

```r
Ef <- predict(gp, X)
```


Manually get the covariance, $K(x_*, x_*) - K(x_*, x)(K(x, x) - \sigma_n^2\mathbb{I})^{-1}K(x,x_*)$. 


```r
var <- gp@kcall$var
var <- 9
Inv <- solve(kernelMatrix(kernelf(gp), xmatrix(gp)) + diag(rep(var, 
    length = length(xmatrix(gp)))))
Cf <- kernelMatrix(kernelf(gp), as.matrix(X)) - kernelMult(kernelf(gp), 
    as.matrix(X), xmatrix(gp), Inv) %*% kernelMatrix(kernelf(gp), xmatrix(gp), 
    as.matrix(X))
```



```r
require(ggplot2)
dat <- data.frame(x = X, y = (Ef), ymin = (Ef - 2 * sqrt(diag(Cf))), 
    ymax = (Ef + 2 * sqrt(diag(Cf))))
ggplot(dat) + geom_ribbon(aes(x = x, y = y, ymin = ymin, ymax = ymax), 
    fill = "grey80") + geom_line(aes(x = x, y = y), size = 1) + geom_point(data = obs, 
    aes(x = x, y = y)) + scale_y_continuous(name = "output, f(x)") + xlab("input, x")
```

![plot of chunk unnamed-chunk-7](http://carlboettiger.info/assets/figures/2012-11-16-e670f9bebf-unnamed-chunk-7.png) 




## Stochastic Dynamic programming solution based on the posterior Gaussian process:


```r
rownorm <- function(M) t(apply(M, 1, function(x) x/sum(x)))
h_grid <- x_grid
```


Define a transition matrix $F$ from the Gaussian process, giving the probability of going from state $x_t$ to $x_{t+1}$.
We already have the Gaussian process mean and variance predicted for each point $x$ on our grid, so this is simply:



```r
V <- sqrt(diag(Cf))
matrices_gp <- lapply(h_grid, function(h) {
    F <- sapply(x_grid, function(x) dnorm(x, Ef - h, V))
    F <- rownorm(F)
})
```


True $f(x)$


```r
matrices_true <- lapply(h_grid, function(h) {
    mu <- sapply(x_grid, f, h, p)
    F_true <- sapply(x_grid, function(x) dnorm(x, mu, sigma_g))
    F_true <- rownorm(F_true)
})
```



Calculate the policy function using the true F:


```r
opt_true <- find_dp_optim(matrices_true, x_grid, h_grid, 20, 0, profit, 
    delta = 0.01)
```


Calculate using the inferred GP


```r
opt_gp <- find_dp_optim(matrices_gp, x_grid, h_grid, 20, 0, profit, 
    delta = 0.01)
```





```r
require(reshape2)
```

```
## Loading required package: reshape2
```

```r
policies <- melt(data.frame(stock = x_grid, GP = x_grid[opt_gp$D[, 
    1]], Exact = x_grid[opt_true$D[, 1]]), id = "stock")
q1 <- ggplot(policies, aes(stock, stock - value, color = variable)) + 
    geom_point() + xlab("stock size") + ylab("escapement")
q1
```

![plot of chunk unnamed-chunk-13](http://carlboettiger.info/assets/figures/2012-11-16-e670f9bebf-unnamed-chunk-13.png) 


We can see what happens when we attempt to manage a stock using this:


```r
z_g <- function() rlnorm(1, 0, sigma_g)
```



```r
set.seed(1)
sim_gp <- ForwardSimulate(f, p, x_grid, h_grid, K, opt_gp$D, z_g, 
    profit = profit)
set.seed(1)
sim_true <- ForwardSimulate(f, p, x_grid, h_grid, K, opt_true$D, 
    z_g, profit = profit)
```



```r
df <- data.frame(time = sim_gp$time, stock_gp = sim_gp$fishstock, 
    stock_true = sim_true$fishstock, harvest_gp = sim_gp$harvest, havest_true = sim_true$harvest)
df <- melt(df, id = "time")
ggplot(df) + geom_line(aes(time, value, color = variable))
```

![plot of chunk simplot](http://carlboettiger.info/assets/figures/2012-11-16-e670f9bebf-simplot.png) 


Total Profits


```r
sum(sim_gp$profit)
```

```
## [1] 27.56
```

```r
sum(sim_true$profit)
```

```
## [1] 29.95
```


