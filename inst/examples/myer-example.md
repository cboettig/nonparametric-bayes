
```
## Loading required package: bibtex
```

```
## Warning: replacing previous import 'write.bib' when loading 'pkgmaker'
```










```r
f <- Myer_harvest
pars <- c(1, 2, 6)
p <- pars  # shorthand
K <- p[1] * p[3]/2 + sqrt((p[1] * p[3])^2 - 4 * p[3])/2
allee <- p[1] * p[3]/2 - sqrt((p[1] * p[3])^2 - 4 * p[3])/2  # allee threshold
e_star <- (p[1] * sqrt(p[3]) - 2)/2  ## Bifurcation point
```


We use the model of Myers _et. al._ (1995).



```r
sigma_g <- 0.02
z_g <- function(sigma_g) rlnorm(1, 0, sigma_g)  #1+(2*runif(1, 0,  1)-1)*sigma_g #
x_grid <- seq(0, 1.5 * K, length = 101)
h_grid <- x_grid
profit = function(x, h) pmin(x, h)
delta <- 0.01
OptTime = 20
xT = allee
reward = profit(x_grid[length(x_grid)], x_grid[length(x_grid)]) + 
    1/(1 - delta)^OptTime
## x_0_observed is starting condition for simulation of the observed data.
## It should be in preferred state for bistable model, above Allee
## threshold for Allee model, and near zero for BH or Ricker models
x_0_observed <- allee + x_grid[1]
```


With parameters `1, 2, 6`. 


```r
Tobs <- 100
x <- numeric(Tobs)
x[1] <- x_0_observed
for (t in 1:(Tobs - 1)) x[t + 1] = z_g(sigma_g) * f(x[t], h = 0, 
    p = p)
plot(x)
```

![plot of chunk sim-obs](figure/sim-obs.png) 



```r
obs <- data.frame(x = x[1:(Tobs - 1)], y = x[2:Tobs])
```



```r
K = mean(x)
r = 0.5 * sigma_g^2/var(x)
f_alt <- Ricker
p_alt <- p
```





```r
gp <- bgp(X = obs$x, XX = x_grid, Z = obs$y, verb = 0, meanfn = "linear", 
    bprior = "b0", BTE = c(2000, 6000, 2), m0r1 = FALSE, corr = "exp", trace = TRUE, 
    beta = c(0, 0), s2.p = c(50, 50), d.p = c(10, 1/0.01, 10, 1/0.01), nug.p = c(10, 
        1/0.01, 10, 1/0.01), s2.lam = "fixed", d.lam = "fixed", nug.lam = "fixed", 
    tau2.lam = "fixed", tau2.p = c(50, 1))
```


We fit a Gaussian process with 


```r
V <- gp$ZZ.ks2
Ef = gp$ZZ.km
tgp_dat <- data.frame(x = gp$XX[[1]], y = gp$ZZ.km, ymin = gp$ZZ.km - 
    1.96 * sqrt(gp$ZZ.ks2), ymax = gp$ZZ.km + 1.96 * sqrt(gp$ZZ.ks2))
```



```r
true <- data.frame(x = x_grid, y = sapply(x_grid, f, 0, p))
ggplot(tgp_dat) + geom_ribbon(aes(x, y, ymin = ymin, ymax = ymax), 
    fill = "gray80") + geom_line(aes(x, y)) + geom_point(data = obs, aes(x, 
    y)) + geom_line(data = true, aes(x, y), col = "red", lty = 2)
```

![plot of chunk gp-plot](figure/gp-plot.png) 



## The transition matrix of the inferred process


```r
X <- numeric(length(x_grid))
X[38] = 1
h <- 0
```



```r
F_ <- gp_F(h, Ef, V, x_grid)
xt1 <- X %*% F_
xt10 <- xt1
for (s in 1:OptTime) xt10 <- xt10 %*% F_
qplot(x_grid, xt10[1, ]) + geom_point(aes(y = xt1[1, ]), col = "grey")
```

![plot of chunk gp-F-sim](figure/gp-F-sim.png) 



```r
F_true <- par_F(h, f, p, x_grid, sigma_g)
yt1 <- X %*% F_true
yt10 <- yt1
for (s in 1:OptTime) yt10 <- yt10 %*% F_true
qplot(x_grid, yt10[1, ]) + geom_point(aes(y = yt1[1, ]), col = "grey")
```

![plot of chunk par-F-sim](figure/par-F-sim.png) 



```r
transition <- melt(data.frame(x = x_grid, gp = xt1[1, ], parametric = yt1[1, 
    ]), id = "x")
ggplot(transition) + geom_point(aes(x, value, col = variable))
```

![plot of chunk F-sim-plot](figure/F-sim-plot.png) 





```r
matrices_gp <- gp_transition_matrix(Ef, V, x_grid, h_grid)
opt_gp <- find_dp_optim(matrices_gp, x_grid, h_grid, OptTime, xT, 
    profit, delta, reward = reward)
```



```r
matrices_true <- f_transition_matrix(f, p, x_grid, h_grid, sigma_g)
opt_true <- find_dp_optim(matrices_true, x_grid, h_grid, OptTime, 
    xT, profit, delta = delta, reward = reward)
```



```r
matrices_estimated <- f_transition_matrix(f_alt, p_alt, x_grid, h_grid, 
    sigma_g)
opt_estimated <- find_dp_optim(matrices_estimated, x_grid, h_grid, 
    OptTime, xT, profit, delta = delta, reward = reward)
```




```r
policies <- melt(data.frame(stock = x_grid, GP = x_grid[opt_gp$D[, 
    1]], Exact = x_grid[opt_true$D[, 1]], Approx = x_grid[opt_estimated$D[, 
    1]]), id = "stock")

policy_plot <- ggplot(policies, aes(stock, stock - value, color = variable)) + 
    geom_point() + xlab("stock size") + ylab("escapement")
policy_plot
```

![plot of chunk policy_plot](figure/policy_plot.png) 




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
set.seed(1)
sim_est <- ForwardSimulate(f, p, x_grid, h_grid, K, opt_estimated$D, 
    z_g, profit = profit)
```




```r
dat <- list(est = sim_est, gp = sim_gp, true = sim_true)
dat <- melt(dat, id = names(dat[[1]]))
dt <- data.table(dat)
setnames(dt, "L1", "method")
```




```r
ggplot(dt) + geom_line(aes(time, fishstock, color = method))
```

![plot of chunk sim-fish](figure/sim-fish.png) 




```r
ggplot(dt) + geom_line(aes(time, harvest, color = method))
```

![plot of chunk sim-harvest](figure/sim-harvest.png) 



```r
c(gp = sum(sim_gp$profit), true = sum(sim_true$profit), est = sum(sim_est$profit))
```

```
##     gp   true    est 
## 13.983 14.196  3.691
```


<p>Myers RA, Barrowman NJ, Hutchings JA and Rosenberg AA (1995).
&ldquo;Population Dynamics of Exploited Fish Stocks at Low Population Levels.&rdquo;
<EM>Science</EM>, <B>269</B>.
ISSN 0036-8075, <a href="http://dx.doi.org/10.1126/science.269.5227.1106">http://dx.doi.org/10.1126/science.269.5227.1106</a>.

