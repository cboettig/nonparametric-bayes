










```r
f <- Myer_harvest
pars <- c(1, 2, 4.5) 
p <- pars # shorthand 
K <- p[1] * p[3] / 2 + sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2
allee <- p[1] * p[3] / 2 - sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2 # allee threshold
e_star <- (p[1] * sqrt(p[3]) - 2) / 2 ## Bifurcation point 
```


We use the model of Myers _et. al._ (1995).



```r
sigma_g <- 0.05
sigma_m <- 0.2
z_g <- function(sigma_g) rlnorm(1, 0, sigma_g) #1+(2*runif(1, 0,  1)-1)*sigma_g #
x_grid <- seq(0, 1.5 * K, length=101)
h_grid <- x_grid
profit = function(x,h) pmin(x, h)
delta <- 0.01
OptTime = 20
reward = profit(x_grid[length(x_grid)], x_grid[length(x_grid)]) + 1 / (1 - delta) ^ OptTime 
xT <- 0
```


With parameters `1, 2, 4.5`. 


```r
x_0_observed <- allee + x_grid[15]
xT <- 0
set.seed(1)
```



```r
Tobs <- 30
x <- numeric(Tobs)
x[1] <- x_0_observed
for(t in 1:(Tobs-1))
  x[t+1] = z_g(sigma_g) * f(x[t], h=0, p=p)
plot(x)
```

![plot of chunk sim-obs](http://carlboettiger.info/assets/figures/2012-12-19-36a9985baf-sim-obs.png) 


We simulate data under this model, starting from a size of `2.13`.  



```r
obs <- data.frame(x=c(0,x[1:(Tobs-1)]),y=c(0,x[2:Tobs]))
```


We consider the observations as ordered pairs of observations of current stock size $x_t$ and observed stock in the following year, $x_{t+1}$.  We add the pseudo-observation of $0,0$.  Alternatively we could condition strictly on solutions passing through the origin, though in practice the weaker assumption is often sufficient. 


```r
estf <- function(p){
  mu <- log(obs$x) + p["r"]*(1-obs$x/p["K"])
  -sum(dlnorm(obs$y, mu, p["s"]), log=TRUE)
}
o <- optim(par = c(r=1,K=mean(x),s=1), estf, method="L", lower=c(1e-3,1e-3,1e-3))
f_alt <- Ricker
p_alt <- c(o$par['r'], o$par['K'])
sigma_g_alt <- o$par['s']
```





```r
gp <- bgp(X=obs$x, XX=x_grid, Z=obs$y, verb=0,
          meanfn="linear", bprior="b0", BTE=c(2000,6000,2), m0r1=FALSE, 
          corr="exp", trace=TRUE, beta = c(0,0),
          s2.p = c(50,50), d.p = c(10, 1/0.01, 10, 1/0.01), nug.p = c(10, 1/0.01, 10, 1/0.01),
          s2.lam = "fixed", d.lam = "fixed", nug.lam = "fixed", 
          tau2.lam = "fixed", tau2.p = c(50,1))
```


We fit a Gaussian process with 


```r
V <- gp$ZZ.ks2
Ef = gp$ZZ.km
tgp_dat <- data.frame(x   = gp$XX[[1]], 
                  y   = gp$ZZ.km, 
                 ymin = gp$ZZ.km - 1.96 * sqrt(gp$ZZ.ks2), 
                 ymax = gp$ZZ.km + 1.96 * sqrt(gp$ZZ.ks2))
```



```r
true <- sapply(x_grid, f, 0, p)
est <- sapply(x_grid, f_alt, 0, p_alt)
models <- data.frame(x=x_grid, GP=tgp_dat$y, Parametric=est, True=true)
models <- melt(models, id="x")
names(models) <- c("x", "method", "value")
# plot
ggplot(tgp_dat)  + geom_ribbon(aes(x,y,ymin=ymin,ymax=ymax), fill="gray80") +
  geom_line(data=models, aes(x, value, col=method), lwd=2, alpha=0.8) + 
  geom_point(data=obs, aes(x,y), alpha=0.8) + 
  xlab(expression(X[t])) + ylab(expression(X[t+1])) +
  scale_colour_manual(values=cbPalette)
```

![plot of chunk gp-plot](http://carlboettiger.info/assets/figures/2012-12-19-36a9985baf-gp-plot.png) 



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
for(s in 1:OptTime)
  xt10 <- xt10 %*% F_
qplot(x_grid, xt10[1,]) + geom_point(aes(y=xt1[1,]), col="grey")
```

![plot of chunk gp-F-sim](http://carlboettiger.info/assets/figures/2012-12-19-36a9985baf-gp-F-sim.png) 



```r
F_true <- par_F(h, f, p, x_grid, sigma_g)
yt1 <- X %*% F_true
yt10 <- yt1
for(s in 1:OptTime)
  yt10 <- yt10 %*% F_true
qplot(x_grid, yt10[1,]) + geom_point(aes(y=yt1[1,]), col="grey")
```

![plot of chunk par-F-sim](http://carlboettiger.info/assets/figures/2012-12-19-36a9985baf-par-F-sim.png) 



```r
transition <- melt(data.frame(x = x_grid, gp = xt1[1,], parametric = yt1[1,]), id="x")
ggplot(transition) + geom_point(aes(x,value, col=variable))
```

![plot of chunk F-sim-plot](http://carlboettiger.info/assets/figures/2012-12-19-36a9985baf-F-sim-plot.png) 





```r
matrices_gp <- gp_transition_matrix(Ef, V, x_grid, h_grid)
opt_gp <- find_dp_optim(matrices_gp, x_grid, h_grid, OptTime, xT, profit, delta, reward=reward)
```



```r
matrices_true <- f_transition_matrix(f, p, x_grid, h_grid, sigma_g)
opt_true <- find_dp_optim(matrices_true, x_grid, h_grid, OptTime, xT, profit, delta=delta, reward = reward)
```



```r
matrices_estimated <- f_transition_matrix(f_alt, p_alt, x_grid, h_grid, sigma_g_alt)
opt_estimated <- find_dp_optim(matrices_estimated, x_grid, h_grid, OptTime, xT, profit, delta=delta, reward = reward)
```




```r
policies <- melt(data.frame(stock=x_grid, 
                            GP = x_grid[opt_gp$D[,1]], 
                            Parametric = x_grid[opt_estimated$D[,1]],
                            True = x_grid[opt_true$D[,1]]),
                  id="stock")
names(policies) <- c("stock", "method", "value")
policy_plot <- ggplot(policies, aes(stock, stock - value, color=method)) +
  geom_line(lwd=2, alpha=0.8) + 
  xlab("stock size") + ylab("escapement")  +
  scale_colour_manual(values=cbPalette)
policy_plot
```

![plot of chunk policy_plot](http://carlboettiger.info/assets/figures/2012-12-19-36a9985baf-policy_plot.png) 




```r
z_g = function() rlnorm(1, 0, sigma_g)
z_m = function() 1+(2*runif(1, 0,  1)-1) * sigma_m
```



```r
m <- sapply(1:OptTime, function(i) opt_gp$D[,1])
opt_gp$D <- m
mm <- sapply(1:OptTime, function(i) opt_true$D[,1])
opt_true$D <- mm
mmm <- sapply(1:OptTime, function(i) opt_estimated$D[,1])
opt_estimated$D <- mmm
```





```r
set.seed(1)
sim_gp <- lapply(1:100, function(i) ForwardSimulate(f, p, x_grid, h_grid, K, opt_gp$D, z_g, profit=profit))
set.seed(1)
sim_true <- lapply(1:100, function(i) ForwardSimulate(f, p, x_grid, h_grid, K, opt_true$D, z_g, profit=profit))
set.seed(1)
sim_est <- lapply(1:100, function(i) ForwardSimulate(f, p, x_grid, h_grid, K, opt_estimated$D, z_g, profit=profit))
```




```r
dat <- list(GP = sim_gp, Parametric = sim_est, True = sim_true)
dat <- melt(dat, id=names(dat[[1]][[1]]))
dt <- data.table(dat)
setnames(dt, c("L1", "L2"), c("method", "reps")) 
```




```r
ggplot(dt) + 
  geom_line(aes(time, fishstock, group=interaction(reps,method), color=method), alpha=.1) +
  scale_colour_manual(values=cbPalette, guide = guide_legend(override.aes = list(alpha = 1)))
```

![plot of chunk sim-fish](http://carlboettiger.info/assets/figures/2012-12-19-36a9985baf-sim-fish.png) 




```r
ggplot(dt) +
  geom_line(aes(time, harvest, group=interaction(reps,method), color=method), alpha=.1) +
  scale_colour_manual(values=cbPalette, guide = guide_legend(override.aes = list(alpha = 1)))
```

![plot of chunk sim-harvest](http://carlboettiger.info/assets/figures/2012-12-19-36a9985baf-sim-harvest.png) 



```r
profits <- dt[, sum(profit), by = c("reps", "method")]
means <- profits[, mean(V1), by = method]
sds <- profits[, sd(V1), by = method]
cbind(means, sd = sds$V1)
```

```
       method    V1     sd
1:         GP 3.000 0.0000
2: Parametric 0.000 0.0000
3:       True 3.114 0.4534
```


<p>Myers RA, Barrowman NJ, Hutchings JA and Rosenberg AA (1995).
&ldquo;Population Dynamics of Exploited Fish Stocks at Low Population Levels.&rdquo;
<EM>Science</EM>, <B>269</B>.
ISSN 0036-8075, <a href="http://dx.doi.org/10.1126/science.269.5227.1106">http://dx.doi.org/10.1126/science.269.5227.1106</a>.

