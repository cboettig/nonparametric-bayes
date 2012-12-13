Reed (1979) with Gaussian processes 
========================================================







```r
knit("~/Documents/code/nonparametric-bayes/inst/examples/beverton_holt_data.Rmd")
```

```


processing file:
~/Documents/code/nonparametric-bayes/inst/examples/beverton_holt_data.Rmd
```

```
output file:
/home/cboettig/Documents/code/nonparametric-bayes/inst/examples/beverton_holt_data.md
```

```r
true <- data.frame(x = X, y = sapply(X, f, 0, p))
```


Predict the function over the target grid


```r
gp <- gp_fit(obs, X, c(sigma_n=2, tau=10, l=8), "conditional")
plot.gpfit(gp, true)
```

![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-12-63cd87e9ac-unnamed-chunk-2.png) 



## Stochastic Dynamic programming solution based on the posterior Gaussian process:


```r
h_grid <- x_grid
```


GP Solution


```r
matrices_gp <- gp_transition_matrix(gp$Ef, diag(gp$Cf), x_grid)
opt_gp <- find_dp_optim(matrices_gp, x_grid, h_grid, 20, 0, profit, delta=.01)
```


Optimal solution


```r
matrices_true <- f_transition_matrix(f, p, x_grid, h_grid, sigma_g)
opt_true <- find_dp_optim(matrices_true, x_grid, h_grid, 20, 0, profit, delta=.01)
```


Optimal parametric solution in which a parameter is estimated with error


```r
k <- 15
A <- 1.5
B <- (A - 1)/k
matrices_estimated <- f_transition_matrix(f, c(A, B), x_grid, h_grid, sigma_g)
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

![plot of chunk policy_plot](http://carlboettiger.info/assets/figures/2012-12-12-63cd87e9ac-policy_plot.png) 


We can see what happens when we attempt to manage a stock using this:


```r
z_g <- function() rlnorm(1,0, sigma_g)
```



```r
set.seed(1)
sim_gp <- ForwardSimulate(f, p, x_grid, h_grid, K, opt_gp$D, z_g, profit=profit)
set.seed(1)
sim_true <- ForwardSimulate(f, p, x_grid, h_grid, K, opt_true$D, z_g, profit=profit)
set.seed(1)
sim_est <- ForwardSimulate(f, p, x_grid, h_grid, K, opt_estimated$D, z_g, profit=profit)
```




```r
dat <- list(est = sim_est, gp = sim_gp, true = sim_true)
dat <- melt(dat, id=names(dat[[1]]))
dt <- data.table(dat)
setnames(dt, "L1", "method") 
```



```r
ggplot(dt) + geom_line(aes(time,fishstock, color=method))
```

![plot of chunk simplot](http://carlboettiger.info/assets/figures/2012-12-12-63cd87e9ac-simplot1.png) 

```r
ggplot(dt) + geom_line(aes(time,harvest, color=method))
```

![plot of chunk simplot](http://carlboettiger.info/assets/figures/2012-12-12-63cd87e9ac-simplot2.png) 


Total Profits


```r
c( gp = sum(sim_gp$profit), true = sum(sim_true$profit), est = sum(sim_est$profit))
```

```
   gp  true   est 
27.60 29.95 26.78 
```


