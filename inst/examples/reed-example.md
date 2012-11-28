Reed (1979) with Gaussian processes 
========================================================







```r
knit("beverton_holt_data.Rmd")
```

```


processing file: beverton_holt_data.Rmd
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
gp <- gp_fit(obs, X, c(sigma_n=.1, l=8), method="kernlab", fit=FALSE)
plot.gpfit(gp, true)
```

![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-11-28-158fae0254-unnamed-chunk-2.png) 



## Stochastic Dynamic programming solution based on the posterior Gaussian process:


```r
h_grid <- x_grid
matrices_gp <- gp_transition_matrix(gp)
matrices_true <- f_transition_matrix(f, p, x_grid, h_grid)
opt_true <- find_dp_optim(matrices_true, x_grid, h_grid, 20, 0, profit, delta=.01)
opt_gp <- find_dp_optim(matrices_gp, x_grid, h_grid, 20, 0, profit, delta=.01)
```



```r
policies <- melt(data.frame(stock=x_grid, GP = x_grid[opt_gp$D[,1]], Exact = x_grid[opt_true$D[,1]]), id="stock")
policy_plot <- ggplot(policies, aes(stock, stock - value, color=variable)) +
  geom_point() + xlab("stock size") + ylab("escapement") 
policy_plot
```

![plot of chunk policy_plot](http://carlboettiger.info/assets/figures/2012-11-28-158fae0254-policy_plot.png) 


We can see what happens when we attempt to manage a stock using this:


```r
z_g <- function() rlnorm(1,0, sigma_g)
```



```r
set.seed(1)
sim_gp <- ForwardSimulate(f, p, x_grid, h_grid, K, opt_gp$D, z_g, profit=profit)
set.seed(1)
sim_true <- ForwardSimulate(f, p, x_grid, h_grid, K, opt_true$D, z_g, profit=profit)
```



```r
df <- data.frame(time = sim_gp$time, stock_gp = sim_gp$fishstock, stock_true = sim_true$fishstock, harvest_gp = sim_gp$harvest, havest_true = sim_true$harvest)
df <- melt(df, id="time")
simplot <- ggplot(df) + geom_line(aes(time,value, color=variable))
simplot
```

![plot of chunk simplot](http://carlboettiger.info/assets/figures/2012-11-28-158fae0254-simplot.png) 


Total Profits


```r
sum(sim_gp$profit)
```

```
[1] 29.33
```

```r
sum(sim_true$profit)
```

```
[1] 29.95
```


