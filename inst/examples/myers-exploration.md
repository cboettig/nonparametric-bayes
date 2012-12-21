









Fixed priors on hyperparameters, fixed model type.


```r
s2.p <- c(50,50)
tau2.p <- c(20,1)
d.p = c(10, 1/0.01, 10, 1/0.01)
nug.p = c(10, 1/0.01, 10, 1/0.01)
s2_prior <- function(x) dinvgamma(x, s2.p[1], s2.p[2])
tau2_prior <- function(x) dinvgamma(x, tau2.p[1], tau2.p[2])
d_prior <- function(x) dgamma(x, d.p[1], scale = d.p[2]) + dgamma(x, d.p[3], scale = d.p[4])
nug_prior <- function(x) dgamma(x, nug.p[1], scale = nug.p[2]) + dgamma(x, nug.p[3], scale = nug.p[4])
beta0_prior <- function(x, tau) dnorm(x, 0, tau)
beta = c(0)
```



```r
profit = function(x,h) pmin(x, h)
delta <- 0.01
OptTime = 20
reward = 0
xT <- 0
z_g = function() rlnorm(1, 0, sigma_g)
z_m = function() 1+(2*runif(1, 0,  1)-1) * sigma_m
```




```r
f <- Myer_harvest
pars <- c(1, 2, 4.5) 
p <- pars # shorthand 
K <- p[1] * p[3] / 2 + sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2
allee <- p[1] * p[3] / 2 - sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2 # allee threshold
e_star <- (p[1] * sqrt(p[3]) - 2) / 2 ## Bifurcation point 
```




```r
sigma_g <- 0.05
sigma_m <- 0.2
x_grid <- seq(0, 1.5 * K, length=101)
h_grid <- x_grid
```


With parameters 1, 2, 4.5. 


```r

sim_obs <- function(Xo, z_g, f, p, Tobs = 40, seed = 1){
  x <- numeric(Tobs)
  x[1] <- Xo
  set.seed(seed)
  for(t in 1:(Tobs-1))
    x[t+1] = z_g() * f(x[t], h=0, p=p)
  data_plot <- qplot(1:Tobs, x) + ggtitle(paste("seed", seed))
  print(data_plot)
  obs <- data.frame(x=c(0,x[1:(Tobs-1)]),y=c(0,x[2:Tobs]))
}

par_est <- function(obs){
  estf <- function(p){
    mu <- log(obs$x) + p["r"]*(1-obs$x/p["K"])
    llik <- -sum(dlnorm(obs$y, mu, p["s"]), log=TRUE)
    if(!is.numeric(llik) | is.nan(llik) | !(llik < Inf)){
      warning("possible error in llik")
      llik <- 1e30
    }
    llik
  }
  o <- optim(par = c(r=1,K=mean(obs$x),s=1), estf, method="L", lower=c(1e-3,1e-3,1e-3))
  f_alt <- Ricker
  p_alt <- c(o$par['r'], o$par['K'])
  sigma_g_alt <- o$par['s']
  list(f_alt = f_alt, p_alt = p_alt, sigma_g_alt = sigma_g_alt)
}
optimal_policy <- function(gp, f, f_alt, p, p_alt, x_grid, h_grid, sigma_g, sigma_g_alt, delta, xT, profit, reward, OptTime){
  matrices_gp <- gp_transition_matrix(gp$ZZ.ks2, gp$ZZ.km, x_grid, h_grid)
  opt_gp <- find_dp_optim(matrices_gp, x_grid, h_grid, OptTime, xT, profit, delta, reward=reward)
  matrices_true <- f_transition_matrix(f, p, x_grid, h_grid, sigma_g)
  opt_true <- find_dp_optim(matrices_true, x_grid, h_grid, OptTime, xT, profit, delta=delta, reward = reward)
  matrices_estimated <- f_transition_matrix(f_alt, p_alt, x_grid, h_grid, sigma_g_alt)
  opt_estimated <- find_dp_optim(matrices_estimated, x_grid, h_grid, OptTime, xT, profit, delta=delta, reward = reward)
  list(gp_D = opt_gp$D[,1], true_D = opt_true$D[,1], est_D = opt_estimated$D[,1])
}

simulate_opt <- function(OPT, f, p, x_grid, h_grid, K, z_g, profit){
  gp_D <- sapply(1:OptTime, function(i) OPT$gp_D)
  true_D <- sapply(1:OptTime, function(i) OPT$true_D)
  est_D <- sapply(1:OptTime, function(i) OPT$est_D)
  set.seed(1)
  sim_gp <- lapply(1:100, function(i) ForwardSimulate(f, p, x_grid, h_grid, K, gp_D, z_g, profit=profit))
  set.seed(1)
  sim_true <- lapply(1:100, function(i) ForwardSimulate(f, p, x_grid, h_grid, K, true_D, z_g, profit=profit))
  set.seed(1)
  sim_est <- lapply(1:100, function(i) ForwardSimulate(f, p, x_grid, h_grid, K, est_D, z_g, profit=profit))
  dat <- list(GP = sim_gp, Parametric = sim_est, True = sim_true)
  dat <- melt(dat, id=names(dat[[1]][[1]]))
  dt <- data.table(dat)
  setnames(dt, c("L1", "L2"), c("method", "reps")) 
  dt
}

sim_plots <- function(dt, seed=1){
  fish_plot <- ggplot(dt) + 
    geom_line(aes(time, fishstock, group=interaction(reps,method), color=method), alpha=.1) +
    scale_colour_manual(values=cbPalette, guide = guide_legend(override.aes = list(alpha = 1)))+
    ggtitle(paste("seed", seed))
  print(fish_plot)
  harvest_plot <- ggplot(dt) +
    geom_line(aes(time, harvest, group=interaction(reps,method), color=method), alpha=.1) +
    scale_colour_manual(values=cbPalette, guide = guide_legend(override.aes = list(alpha = 1))) +
    ggtitle(paste("seed", seed))
  print(harvest_plot)
}

profits_stats <- function(dt){
  profits <- dt[, sum(profit), by = c("reps", "method")]
  means <- profits[, mean(V1), by = method]
  sds <- profits[, sd(V1), by = method]
  yield <- cbind(means, sd = sds$V1)
  yield
}
```







```r
require(snowfall)
sfInit(par=TRUE, cpu=8)
sfExportAll()
sfLibrary(pdgControl)
sfLibrary(nonparametricbayes)
sfLibrary(reshape2)
sfLibrary(ggplot2)
sfLibrary(data.table)
sfLibrary(tgp)
sfLibrary(kernlab)
sfLibrary(MCMCpack)
sfLibrary(plyr)
```




```r
seed <- round(runif(16) * 1e6)
seed
```

```
 [1] 284002 748255 612767 512801 874437  13282 555496 924097 535430 227849
[11]   4267 674742 111760 114786 773117 725971
```

```r
yields <- lapply(seed, 
function(seed_i){
  Xo <- allee + x_grid[5]
  obs <- sim_obs(Xo, z_g, f, p, Tobs=50, seed = seed_i)
  est <- par_est(obs)
  gp <- bgp(X=obs$x, XX=x_grid, Z=obs$y, verb=0,
          meanfn="constant", bprior="b0", BTE=c(2000,6000,2),
          m0r1=FALSE, corr="exp", trace=FALSE, 
          beta = beta, s2.p = s2.p, d.p = d.p, nug.p = nug.p, tau2.p = tau2.p,
          s2.lam = "fixed", d.lam = "fixed", nug.lam = "fixed", tau2.lam = "fixed")      
  OPT <- optimal_policy(gp, f, est$f_alt, p, est$p_alt, x_grid, h_grid, sigma_g, est$sigma_g_alt, delta, xT, profit, reward, OptTime)
  dt <- simulate_opt(OPT, f, p, x_grid, h_grid, K, z_g, profit)
  sim_plots(dt, seed=seed_i)
  profits_stats(dt)
})
```

![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-25-278d66f2af-unnamed-chunk-21.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-26-278d66f2af-unnamed-chunk-22.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-26-278d66f2af-unnamed-chunk-23.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-28-278d66f2af-unnamed-chunk-24.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-28-278d66f2af-unnamed-chunk-25.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-29-278d66f2af-unnamed-chunk-26.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-30-278d66f2af-unnamed-chunk-27.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-30-278d66f2af-unnamed-chunk-28.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-31-278d66f2af-unnamed-chunk-29.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-37-278d66f2af-unnamed-chunk-210.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-38-278d66f2af-unnamed-chunk-211.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-38-278d66f2af-unnamed-chunk-212.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-39-278d66f2af-unnamed-chunk-213.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-40-278d66f2af-unnamed-chunk-214.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-41-278d66f2af-unnamed-chunk-215.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-41-278d66f2af-unnamed-chunk-216.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-42-278d66f2af-unnamed-chunk-217.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-43-278d66f2af-unnamed-chunk-218.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-44-278d66f2af-unnamed-chunk-219.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-44-278d66f2af-unnamed-chunk-220.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-45-278d66f2af-unnamed-chunk-221.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-46-278d66f2af-unnamed-chunk-222.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-47-278d66f2af-unnamed-chunk-223.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-52-278d66f2af-unnamed-chunk-224.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-53-278d66f2af-unnamed-chunk-225.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-54-278d66f2af-unnamed-chunk-226.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-55-278d66f2af-unnamed-chunk-227.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-55-278d66f2af-unnamed-chunk-228.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-56-278d66f2af-unnamed-chunk-229.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-57-278d66f2af-unnamed-chunk-230.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-58-278d66f2af-unnamed-chunk-231.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-58-278d66f2af-unnamed-chunk-232.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-33-59-278d66f2af-unnamed-chunk-233.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-34-00-278d66f2af-unnamed-chunk-234.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-34-01-278d66f2af-unnamed-chunk-235.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-34-01-278d66f2af-unnamed-chunk-236.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-34-02-278d66f2af-unnamed-chunk-237.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-34-03-278d66f2af-unnamed-chunk-238.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-34-09-278d66f2af-unnamed-chunk-239.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-34-09-278d66f2af-unnamed-chunk-240.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-34-10-278d66f2af-unnamed-chunk-241.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-34-11-278d66f2af-unnamed-chunk-242.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-34-12-278d66f2af-unnamed-chunk-243.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-34-12-278d66f2af-unnamed-chunk-244.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-34-13-278d66f2af-unnamed-chunk-245.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-34-14-278d66f2af-unnamed-chunk-246.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-34-14-278d66f2af-unnamed-chunk-247.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-16-34-15-278d66f2af-unnamed-chunk-248.png) 



```r
yields <- melt(yields, id=c("method", "V1", "sd"))
yields         
```

```
       method     V1      sd L1
1          GP 3.0000 0.00000  1
2  Parametric 0.0000 0.00000  1
3        True 3.1140 0.45340  1
4          GP 3.0000 0.00000  2
5  Parametric 1.8684 0.11215  2
6        True 3.1140 0.45340  2
7          GP 3.0000 0.00000  3
8  Parametric 1.7221 0.26300  3
9        True 3.1140 0.45340  3
10         GP 3.0000 0.00000  4
11 Parametric 1.8847 0.13847  4
12       True 3.1140 0.45340  4
13         GP 3.0000 0.00000  5
14 Parametric 1.8503 0.05496  5
15       True 3.1140 0.45340  5
16         GP 3.0000 0.00000  6
17 Parametric 1.9166 0.19041  6
18       True 3.1140 0.45340  6
19         GP 3.0000 0.00000  7
20 Parametric 1.7435 0.06622  7
21       True 3.1140 0.45340  7
22         GP 3.0000 0.00000  8
23 Parametric 1.9053 0.12432  8
24       True 3.1140 0.45340  8
25         GP 3.0000 0.00000  9
26 Parametric 1.6916 0.13359  9
27       True 3.1140 0.45340  9
28         GP 3.0000 0.00000 10
29 Parametric 1.9237 0.04701 10
30       True 3.1140 0.45340 10
31         GP 3.0000 0.00000 11
32 Parametric 0.0315 0.06538 11
33       True 3.1140 0.45340 11
34         GP 3.0000 0.00000 12
35 Parametric 1.7746 0.09434 12
36       True 3.1140 0.45340 12
37         GP 3.0000 0.00000 13
38 Parametric 2.0467 0.34409 13
39       True 3.1140 0.45340 13
40         GP 3.0000 0.00000 14
41 Parametric 0.0000 0.00000 14
42       True 3.1140 0.45340 14
43         GP 3.0000 0.00000 15
44 Parametric 1.7403 0.06704 15
45       True 3.1140 0.45340 15
46         GP 3.0000 0.00000 16
47 Parametric 1.6750 0.11505 16
48       True 3.1140 0.45340 16
```




NULL

