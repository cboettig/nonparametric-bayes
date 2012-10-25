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

![plot of chunk unnamed-chunk-5](http://farm9.staticflickr.com/8464/8123212898_904aa88749_o.png) 



## Stochastic Dynamic programming solution based on the posterior Gaussian process:


```r
rownorm <- function(M) t(apply(M, 1, function(x) x/sum(x)))
h_grid <- x_grid
```


Define a transition matrix $F$ from the Gaussian process, giving the probability of going from state $x_t$ to $x_{t+1}$.
We already have the Gaussian process mean and variance predicted for each point $x$ on our grid, so this is simply:



```r
V <- sqrt(diag(gp$Cf))
matrices_gp <- lapply(h_grid, function(h) {
    F <- sapply(x_grid, function(x) dnorm(x, gp$Ef - h, V))
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

![plot of chunk unnamed-chunk-11](http://farm9.staticflickr.com/8053/8123196485_a3dc43986f_o.png) 


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

![plot of chunk simplot](http://farm9.staticflickr.com/8191/8123196619_5b8ab8d81d_o.png) 


Total Profits


```r
sum(sim_gp$profit)
```

```
## [1] 25.43
```

```r
sum(sim_true$profit)
```

```
## [1] 29.95
```


