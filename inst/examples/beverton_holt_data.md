Generate data under Beverton Holt model
==============================================


```r
require(pdgControl)
require(nonparametricbayes)
require(ggplot2)
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
```







