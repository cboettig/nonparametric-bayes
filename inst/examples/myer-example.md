





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
require(reshape2)
```

```
## Loading required package: reshape2
```

```r
require(ggplot2)
```

```
## Loading required package: ggplot2
```

```r
require(data.table)
```

```
## Loading required package: data.table
```

```r
library(kernlab)
```








```r
f <- Myer_harvest
pars <- c(1, 2, 6) 
p <- pars # shorthand 
K <- p[1] * p[3] / 2 + sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2
xT <- p[1] * p[3] / 2 - sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2 # allee threshold
e_star <- (p[1] * sqrt(p[3]) - 2) / 2 ## Bifurcation point 
```



Noise function, profit function

```r
z_g <- function(sigma_g) rlnorm(1, 0, sigma_g) #1+(2*runif(1, 0,  1)-1)*sigma_g #
profit <- profit_harvest(1,0,0)
```



Parameter definitions


```r
x_grid = seq(0, 1.5 * K, length=101)
T <- 40
sigma_g <- 0.2
x <- numeric(T)
x[1] <- 1
seed <- 123
```




Simulation 


```r
x[1] = 2.5
for(t in 1:(T-1))
  x[t+1] = z_g(sigma_g) * f(x[t], h=0, p=p)
plot(x)
```

![plot of chunk unnamed-chunk-4](http://carlboettiger.info/assets/figures/2012-12-07-92988ba826-unnamed-chunk-4.png) 


Predict the function over the target grid


```r
obs <- data.frame(x=x[1:(T-1)],y=x[2:T])
X <- x_grid
gp <- gp_fit(obs, X, c(sigma_n=.9, tau = 1., l=.1), "kernlab")
```

```
Using automatic sigma estimation (sigest) for RBF or laplace kernel 
```


Gaussian Process inference from this model


```r
df <- data.frame(x=X, y=gp$Ef, ymin=(gp$Ef-2*sqrt(abs(diag(gp$Cf)))), ymax=(gp$Ef+2*sqrt(abs(diag(gp$Cf)))))
true <- data.frame(x=X, y=sapply(X,f, 0, p))
ggplot(df)  + geom_ribbon(aes(x,y,ymin=ymin,ymax=ymax), fill="gray80") +
  geom_line(aes(x,y)) + geom_point(data=obs, aes(x,y)) +
  geom_line(data=true, aes(x,y), col='red', lty=2)
```

![plot of chunk unnamed-chunk-6](http://carlboettiger.info/assets/figures/2012-12-07-92988ba826-unnamed-chunk-6.png) 




## Stochastic Dynamic programming solution based on the posterior Gaussian process:


```r
h_grid <- x_grid
```


GP Solution


```r
matrices_gp <- gp_transition_matrix(gp)
opt_gp <- find_dp_optim(matrices_gp, x_grid, h_grid, 20, 0, profit, delta=.01)
```


Optimal solution


```r
matrices_true <- f_transition_matrix(f, p, x_grid, h_grid, sigma_g)
opt_true <- find_dp_optim(matrices_true, x_grid, h_grid, 20, 0, profit, delta=.01)
```


Optimal parametric solution in which a parameter is estimated with error


```r
p_est = c(1.5,2,6)
matrices_estimated <- f_transition_matrix(f, p_est, x_grid, h_grid, sigma_g=sigma_g)
opt_estimated <- find_dp_optim(matrices_estimated, x_grid, h_grid, 20, 0, profit, delta=.01)
```




```r
policies <- melt(data.frame(stock=x_grid, 
                            GP = x_grid[opt_gp$D[,1]], 
                            Exact = x_grid[opt_true$D[,1]],
                            Approx = x_grid[opt_estimated$D[,1]]),
                  id="stock")

policy_plot <- ggplot(policies, aes(stock, stock - value, color=variable)) +
  geom_point() + xlab("stock size") + ylab("escapement") 
policy_plot
```

![plot of chunk policy_plot](http://carlboettiger.info/assets/figures/2012-12-07-92988ba826-policy_plot.png) 


We can see what happens when we attempt to manage a stock using this:


```r
z_g <- function() rlnorm(1,0, sigma_g)
```



```r
set.seed(seed)
sim_gp <- lapply(1:20, function(i) ForwardSimulate(f, p, x_grid, h_grid, K, opt_gp$D, z_g, profit=profit))
set.seed(seed)
sim_true <- lapply(1:20, function(i) ForwardSimulate(f, p, x_grid, h_grid, K, opt_true$D, z_g, profit=profit))
set.seed(seed)
sim_est <- lapply(1:20, function(i) ForwardSimulate(f, p, x_grid, h_grid, K, opt_estimated$D, z_g, profit=profit))
```




```r
dat <- list(est = sim_est, gp = sim_gp, true = sim_true)
dat <- melt(dat, id=names(dat[[1]][[1]]))
dt <- data.table(dat)
setnames(dt, c("L1", "L2"),  c("method", "reps")) 
```


And plot the resulting stock dynamcis across the replicates [^1]


```r
ggplot(dt) + geom_line(aes(time,fishstock, group=reps, color=method)) + facet_wrap(~method)
```

![plot of chunk simplot](http://carlboettiger.info/assets/figures/2012-12-07-92988ba826-simplot.png) 



[^1]: Not sure why that doesn't plot correctly as overlay when I drop the `facet_wrap` option from the command.


Total Profits


```r
c( gp = sum(sim_gp$profit), true = sum(sim_true$profit), est = sum(sim_est$profit))
```

```
  gp true  est 
   0    0    0 
```


