









```r
knit("beverton_holt_data.Rmd")
```


```r
true <- data.frame(x = X, y = sapply(X, f, 0, p))
```





```r
gp <- bgp(X=obs$x, XX=X, Z=obs$y)
```

```

burn in:
r=1000 d=[0.157057]; n=39

Sampling @ nn=101 pred locs:
r=1000 d=[0.0910807]; mh=1 n=39
r=2000 d=[1.03849]; mh=1 n=39
r=3000 d=[1.30774]; mh=1 n=39

```


Extract the posterior Gaussian process mean and the $\pm 2$ standard deviations over the predicted grid from the fit:


```r
center <- tgp:::tgp.choose.center(gp, "mean")
Z.mean <- center$Z
X_vals <- center$X[, 1]
o <- order(X_vals)
Zb.q1 <- Z.mean + 1.96 * sqrt(c(gp$Zp.ks2, gp$ZZ.ks2))
Zb.q2 <- Z.mean - 1.96 * sqrt(c(gp$Zp.ks2, gp$ZZ.ks2))
V <- c(gp$Zp.ks2, gp$ZZ.ks2)[o]
dat <- data.frame(x=X_vals[o], y=Z.mean[o], ymin=Zb.q1[o], ymax=Zb.q2[o])
```


Plot the posterior Gaussian Process:


```r
  p1 <- ggplot(dat) +
    geom_ribbon(aes(x=x,y=y, ymin=ymin, ymax=ymax), fill="grey80") + # Var
    geom_line(aes(x=x,y=y), size=1) + #MEAN
    geom_point(data=gp$obs,aes(x=x,y=y)) + 
    labs(title=paste("llik =", prettyNum(gp$llik)))
  if(!is.null(true))
    p1 <- p1 + geom_line(data = true, aes(x, y), col = "red", lty = 2) 
p1 + theme_notebook
```

![plot of chunk gp-plot](http://carlboettiger.info/assets/figures/2012-11-29-71f2443754-gp-plot.png) 



Find the actual optimum policy from the true parametric model


```r
x_grid <- dat$x
h_grid <- x_grid
matrices_true <- f_transition_matrix(f, p, x_grid, h_grid)
opt_true <- find_dp_optim(matrices_true, x_grid, h_grid, 20, 0, profit, delta=.01)
```



Estimate a policy from the Gaussian Process:


```r
rownorm <- function(M) t(apply(M, 1, function(x) x/sum(x)))
matrices_gp <- lapply(h_grid, function(h){
    S <- dat$y - h
    F_ <- t(sapply(1:length(S), function(i){
      if(S[i]>0) {
        out <- dlnorm(x_grid/S[i], 0, V[i])
      } else {
        out <- numeric(length(x_grid))
        out[1] <- 1
        out
      }
    }))
    F_ <- rownorm(F_)
  })
opt_gp <- find_dp_optim(matrices_gp, x_grid, h_grid, 20, 0, profit, delta=.01)
```




```r
matrices_true[[1]][1:10,1:10]
```

```
      [,1]       [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
 [1,]    1  0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00
 [2,]    0  5.361e-02 9.464e-01 1.078e-09 2.276e-20 3.993e-31 1.614e-41
 [3,]    0  9.446e-26 7.322e-04 9.915e-01 7.731e-03 5.996e-07 6.519e-12
 [4,]    0  8.376e-48 6.260e-14 8.679e-04 6.371e-01 3.564e-01 5.501e-03
 [5,]    0  1.195e-67 2.460e-25 2.951e-10 6.897e-04 2.011e-01 5.152e-01
 [6,]    0  3.004e-85 1.951e-36 1.482e-17 1.725e-08 6.220e-04 8.162e-02
 [7,]    0 6.940e-101 8.422e-47 7.759e-25 1.392e-13 2.499e-07 7.986e-04
 [8,]    0 1.421e-110 1.840e-53 1.001e-29 3.243e-17 5.493e-10 1.099e-05
 [9,]    0 3.645e-115 1.177e-56 4.205e-32 5.179e-19 2.471e-11 1.152e-06
[10,]    0 4.393e-128 9.076e-66 5.459e-39 2.554e-24 2.047e-15 9.566e-10
           [,8]      [,9]     [,10]
 [1,] 0.000e+00 0.000e+00 0.000e+00
 [2,] 3.496e-48 1.969e-51 7.453e-61
 [3,] 1.939e-15 3.099e-17 1.112e-22
 [4,] 1.085e-04 1.210e-05 8.833e-09
 [5,] 1.950e-01 8.535e-02 2.634e-03
 [6,] 3.003e-01 3.769e-01 2.078e-01
 [7,] 1.860e-02 5.485e-02 3.134e-01
 [8,] 7.383e-04 3.557e-03 7.788e-02
 [9,] 1.262e-04 7.628e-04 3.103e-02
[10,] 3.971e-07 4.447e-06 9.789e-04
```

```r
matrices_gp[[1]][1:10,1:10]
```

```
      [,1]    [,2]    [,3]    [,4]    [,5]    [,6]    [,7]    [,8]    [,9]
 [1,]    0 0.12773 0.09762 0.07549 0.06015 0.04916 0.04101 0.03668 0.03477
 [2,]    0 0.11133 0.09047 0.07220 0.05872 0.04870 0.04109 0.03698 0.03516
 [3,]    0 0.09043 0.07983 0.06658 0.05576 0.04726 0.04055 0.03684 0.03518
 [4,]    0 0.08944 0.08042 0.06749 0.05665 0.04805 0.04122 0.03744 0.03574
 [5,]    0 0.06990 0.06888 0.06072 0.05268 0.04579 0.04004 0.03676 0.03527
 [6,]    0 0.05668 0.06033 0.05539 0.04937 0.04377 0.03887 0.03600 0.03467
 [7,]    0 0.05130 0.05700 0.05346 0.04827 0.04319 0.03862 0.03590 0.03463
 [8,]    0 0.02259 0.03933 0.04504 0.04556 0.04380 0.04110 0.03912 0.03811
 [9,]    0 0.04476 0.05239 0.05045 0.04635 0.04199 0.03790 0.03541 0.03424
[10,]    0 0.03832 0.04739 0.04698 0.04399 0.04040 0.03687 0.03465 0.03359
        [,10]
 [1,] 0.02988
 [2,] 0.03044
 [3,] 0.03081
 [4,] 0.03127
 [5,] 0.03127
 [6,] 0.03107
 [7,] 0.03116
 [8,] 0.03511
 [9,] 0.03101
[10,] 0.03064
```




```r
policies <- melt(data.frame(stock=x_grid, GP = x_grid[opt_gp$D[,1]], Exact = x_grid[opt_true$D[,1]]), id="stock")
policy_plot <- ggplot(policies, aes(stock, stock - value, color=variable)) +
  geom_point() + xlab("stock size") + ylab("escapement") 
policy_plot + theme_notebook
```

![plot of chunk policy_plot](http://carlboettiger.info/assets/figures/2012-11-29-71f2443754-policy_plot.png) 


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
simplot + theme_notebook
```

![plot of chunk simplot](http://carlboettiger.info/assets/figures/2012-11-29-71f2443754-simplot.png) 


Total Profits


```r
sum(sim_gp$profit)
```

```
[1] 17.66
```

```r
sum(sim_true$profit)
```

```
[1] 29.64
```





