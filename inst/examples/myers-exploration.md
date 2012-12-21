









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
OptTime = 20  # stationarity with unstable models is tricky thing
reward = 0
xT <- 0
z_g = function() rlnorm(1, 0, sigma_g)
z_m = function() 1+(2*runif(1, 0,  1)-1) * sigma_m
```




```r
f <- RickerAllee
p <- c(1.5, 10, 5) 
K <- 10
allee <- 5
# valid only for p[2] = 2
#K <- p[1] * p[3] / 2 + sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2 #
#allee <- p[1] * p[3] / 2 - sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2 # allee threshold
#e_star <- (p[1] * sqrt(p[3]) - 2) / 2 ## Bifurcation point 
```




```r
sigma_g <- 0.05
sigma_m <- 0.2
x_grid <- seq(0, 1.5 * K, length=101)
h_grid <- x_grid
```


With parameters 1.5, 10, 5. 


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
  obs
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
  matrices_gp <- gp_transition_matrix(gp$ZZ.km, gp$ZZ.ks2, x_grid, h_grid)
  opt_gp <- find_dp_optim(matrices_gp, x_grid, h_grid, OptTime, xT, profit, delta, reward=reward)
  matrices_true <- f_transition_matrix(f, p, x_grid, h_grid, sigma_g)
  opt_true <- find_dp_optim(matrices_true, x_grid, h_grid, OptTime, xT, profit, delta=delta, reward = reward)
  matrices_estimated <- f_transition_matrix(f_alt, p_alt, x_grid, h_grid, sigma_g_alt)
  opt_estimated <- find_dp_optim(matrices_estimated, x_grid, h_grid, OptTime, xT, profit, delta=delta, reward = reward)
  list(gp_D = opt_gp$D[,1], true_D = opt_true$D[,1], est_D = opt_estimated$D[,1])
}

simulate_opt <- function(OPT, f, p, x_grid, h_grid, x0, z_g, profit){
  gp_D <- sapply(1:OptTime, function(i) OPT$gp_D)
  true_D <- sapply(1:OptTime, function(i) OPT$true_D)
  est_D <- sapply(1:OptTime, function(i) OPT$est_D)
  set.seed(1)
  sim_gp <- lapply(1:100, function(i) ForwardSimulate(f, p, x_grid, h_grid, x0, gp_D, z_g, profit=profit))
  set.seed(1)
  sim_true <- lapply(1:100, function(i) ForwardSimulate(f, p, x_grid, h_grid, x0, true_D, z_g, profit=profit))
  set.seed(1)
  sim_est <- lapply(1:100, function(i) ForwardSimulate(f, p, x_grid, h_grid, x0, est_D, z_g, profit=profit))
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


gp_plot <- function(gp, f, p, f_alt, p_alt, x_grid, obs, seed){
  tgp_dat <- 
  data.frame(  x = gp$XX[[1]], 
               y = gp$ZZ.km, 
             ymin = gp$ZZ.km - 1.96 * sqrt(gp$ZZ.ks2), 
             ymax = gp$ZZ.km + 1.96 * sqrt(gp$ZZ.ks2))
  true <- sapply(x_grid, f, 0, p)
  est <- sapply(x_grid, f_alt, 0, p_alt)
  models <- data.frame(x=x_grid, GP=tgp_dat$y, Parametric=est, True=true)
  models <- melt(models, id="x")
  names(models) <- c("x", "method", "value")
  plot_gp <- ggplot(tgp_dat) + geom_ribbon(aes(x,y,ymin=ymin,ymax=ymax), fill="gray80") +
    geom_line(data=models, aes(x, value, col=method), lwd=2, alpha=0.8) + 
    geom_point(data=obs, aes(x,y), alpha=0.8) + 
    xlab(expression(X[t])) + ylab(expression(X[t+1])) +
    scale_colour_manual(values=cbPalette) +
    ggtitle(paste("seed", seed))
  print(plot_gp)
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
seed <- round(runif(32) * 1e6)
seed[1] <- 1
seed
```

```
 [1]      1 497051 697457 301496 394253  66714 909713 807430 589571 153604
[11] 122498 193816 971561 798732 658520 813171 538883 824588 486796 991560
[21] 602497 309354 389353  42820 976787 607364 969921 197732 876143 169407
[31] 421958 290102
```

```r
yields <- lapply(seed, 
function(seed_i){
  Xo <- allee + x_grid[5] # observations start from
  x0 <- 2 # simulation under policy starts from
  obs <- sim_obs(Xo, z_g, f, p, Tobs=40, seed = seed_i)
  est <- par_est(obs)
  gp <- bgp(X=obs$x, XX=x_grid, Z=obs$y, verb=0,
          meanfn="constant", bprior="b0", BTE=c(2000,16000,2),
          m0r1=FALSE, corr="exp", trace=FALSE, 
          beta = beta, s2.p = s2.p, d.p = d.p, nug.p = nug.p, tau2.p = tau2.p,
          s2.lam = "fixed", d.lam = "fixed", nug.lam = "fixed", tau2.lam = "fixed")      
  gp_plot(gp, f, p, est$f_alt, est$p_alt, x_grid, obs, seed_i)
  OPT <- optimal_policy(gp, f, est$f_alt, p, est$p_alt, x_grid, h_grid, sigma_g, est$sigma_g_alt, delta, xT, profit, reward, OptTime)
  dt <- simulate_opt(OPT, f, p, x_grid, h_grid, x0, z_g, profit)
  sim_plots(dt, seed=seed_i)
  profits_stats(dt)
})
```

![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-10-ee54f02a6d-unnamed-chunk-21.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-11-ee54f02a6d-unnamed-chunk-22.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-12-ee54f02a6d-unnamed-chunk-23.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-12-ee54f02a6d-unnamed-chunk-24.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-13-ee54f02a6d-unnamed-chunk-25.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-14-ee54f02a6d-unnamed-chunk-26.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-14-ee54f02a6d-unnamed-chunk-27.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-15-ee54f02a6d-unnamed-chunk-28.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-16-ee54f02a6d-unnamed-chunk-29.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-16-ee54f02a6d-unnamed-chunk-210.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-17-ee54f02a6d-unnamed-chunk-211.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-18-ee54f02a6d-unnamed-chunk-212.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-18-ee54f02a6d-unnamed-chunk-213.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-19-ee54f02a6d-unnamed-chunk-214.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-20-ee54f02a6d-unnamed-chunk-215.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-21-ee54f02a6d-unnamed-chunk-216.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-21-ee54f02a6d-unnamed-chunk-217.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-22-ee54f02a6d-unnamed-chunk-218.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-23-ee54f02a6d-unnamed-chunk-219.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-28-ee54f02a6d-unnamed-chunk-220.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-29-ee54f02a6d-unnamed-chunk-221.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-30-ee54f02a6d-unnamed-chunk-222.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-30-ee54f02a6d-unnamed-chunk-223.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-31-ee54f02a6d-unnamed-chunk-224.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-32-ee54f02a6d-unnamed-chunk-225.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-32-ee54f02a6d-unnamed-chunk-226.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-33-ee54f02a6d-unnamed-chunk-227.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-34-ee54f02a6d-unnamed-chunk-228.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-34-ee54f02a6d-unnamed-chunk-229.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-35-ee54f02a6d-unnamed-chunk-230.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-36-ee54f02a6d-unnamed-chunk-231.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-36-ee54f02a6d-unnamed-chunk-232.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-37-ee54f02a6d-unnamed-chunk-233.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-38-ee54f02a6d-unnamed-chunk-234.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-43-ee54f02a6d-unnamed-chunk-235.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-44-ee54f02a6d-unnamed-chunk-236.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-45-ee54f02a6d-unnamed-chunk-237.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-45-ee54f02a6d-unnamed-chunk-238.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-46-ee54f02a6d-unnamed-chunk-239.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-47-ee54f02a6d-unnamed-chunk-240.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-47-ee54f02a6d-unnamed-chunk-241.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-48-ee54f02a6d-unnamed-chunk-242.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-49-ee54f02a6d-unnamed-chunk-243.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-50-ee54f02a6d-unnamed-chunk-244.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-50-ee54f02a6d-unnamed-chunk-245.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-51-ee54f02a6d-unnamed-chunk-246.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-52-ee54f02a6d-unnamed-chunk-247.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-52-ee54f02a6d-unnamed-chunk-248.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-58-ee54f02a6d-unnamed-chunk-249.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-59-ee54f02a6d-unnamed-chunk-250.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-29-59-ee54f02a6d-unnamed-chunk-251.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-00-ee54f02a6d-unnamed-chunk-252.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-01-ee54f02a6d-unnamed-chunk-253.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-01-ee54f02a6d-unnamed-chunk-254.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-02-ee54f02a6d-unnamed-chunk-255.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-03-ee54f02a6d-unnamed-chunk-256.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-04-ee54f02a6d-unnamed-chunk-257.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-04-ee54f02a6d-unnamed-chunk-258.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-05-ee54f02a6d-unnamed-chunk-259.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-06-ee54f02a6d-unnamed-chunk-260.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-06-ee54f02a6d-unnamed-chunk-261.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-07-ee54f02a6d-unnamed-chunk-262.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-08-ee54f02a6d-unnamed-chunk-263.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-13-ee54f02a6d-unnamed-chunk-264.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-14-ee54f02a6d-unnamed-chunk-265.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-15-ee54f02a6d-unnamed-chunk-266.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-15-ee54f02a6d-unnamed-chunk-267.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-16-ee54f02a6d-unnamed-chunk-268.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-17-ee54f02a6d-unnamed-chunk-269.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-17-ee54f02a6d-unnamed-chunk-270.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-18-ee54f02a6d-unnamed-chunk-271.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-19-ee54f02a6d-unnamed-chunk-272.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-19-ee54f02a6d-unnamed-chunk-273.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-20-ee54f02a6d-unnamed-chunk-274.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-21-ee54f02a6d-unnamed-chunk-275.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-22-ee54f02a6d-unnamed-chunk-276.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-22-ee54f02a6d-unnamed-chunk-277.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-23-ee54f02a6d-unnamed-chunk-278.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-24-ee54f02a6d-unnamed-chunk-279.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-24-ee54f02a6d-unnamed-chunk-280.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-25-ee54f02a6d-unnamed-chunk-281.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-26-ee54f02a6d-unnamed-chunk-282.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-26-ee54f02a6d-unnamed-chunk-283.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-27-ee54f02a6d-unnamed-chunk-284.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-28-ee54f02a6d-unnamed-chunk-285.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-28-ee54f02a6d-unnamed-chunk-286.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-29-ee54f02a6d-unnamed-chunk-287.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-30-ee54f02a6d-unnamed-chunk-288.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-30-ee54f02a6d-unnamed-chunk-289.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-31-ee54f02a6d-unnamed-chunk-290.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-32-ee54f02a6d-unnamed-chunk-291.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-32-ee54f02a6d-unnamed-chunk-292.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-33-ee54f02a6d-unnamed-chunk-293.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-34-ee54f02a6d-unnamed-chunk-294.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-34-ee54f02a6d-unnamed-chunk-295.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-35-ee54f02a6d-unnamed-chunk-296.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-36-ee54f02a6d-unnamed-chunk-297.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-36-ee54f02a6d-unnamed-chunk-298.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-38-ee54f02a6d-unnamed-chunk-299.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-43-ee54f02a6d-unnamed-chunk-2100.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-44-ee54f02a6d-unnamed-chunk-2101.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-45-ee54f02a6d-unnamed-chunk-2102.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-45-ee54f02a6d-unnamed-chunk-2103.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-46-ee54f02a6d-unnamed-chunk-2104.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-47-ee54f02a6d-unnamed-chunk-2105.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-47-ee54f02a6d-unnamed-chunk-2106.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-48-ee54f02a6d-unnamed-chunk-2107.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-49-ee54f02a6d-unnamed-chunk-2108.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-49-ee54f02a6d-unnamed-chunk-2109.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-50-ee54f02a6d-unnamed-chunk-2110.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-51-ee54f02a6d-unnamed-chunk-2111.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-51-ee54f02a6d-unnamed-chunk-2112.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-52-ee54f02a6d-unnamed-chunk-2113.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-53-ee54f02a6d-unnamed-chunk-2114.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-53-ee54f02a6d-unnamed-chunk-2115.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-54-ee54f02a6d-unnamed-chunk-2116.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-55-ee54f02a6d-unnamed-chunk-2117.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-55-ee54f02a6d-unnamed-chunk-2118.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-56-ee54f02a6d-unnamed-chunk-2119.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-57-ee54f02a6d-unnamed-chunk-2120.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-57-ee54f02a6d-unnamed-chunk-2121.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-58-ee54f02a6d-unnamed-chunk-2122.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-30-59-ee54f02a6d-unnamed-chunk-2123.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-31-00-ee54f02a6d-unnamed-chunk-2124.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-31-00-ee54f02a6d-unnamed-chunk-2125.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-31-01-ee54f02a6d-unnamed-chunk-2126.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-31-02-ee54f02a6d-unnamed-chunk-2127.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-31-02-ee54f02a6d-unnamed-chunk-2128.png) 



```r
yields <- melt(yields, id=c("method", "V1", "sd"))
ggplot(yields) + geom_histogram(aes(V1)) + facet_wrap(~method, scale='free')
```

```
stat_bin: binwidth defaulted to range/30. Use 'binwidth = x' to adjust
this.
```

```
stat_bin: binwidth defaulted to range/30. Use 'binwidth = x' to adjust
this.
```

```
stat_bin: binwidth defaulted to range/30. Use 'binwidth = x' to adjust
this.
```

![plot of chunk unnamed-chunk-3](http://carlboettiger.info/assets/figures/2012-12-20-20-31-03-ee54f02a6d-unnamed-chunk-3.png) 

```r
yields
```

```
       method     V1       sd L1
1          GP 1.9500 0.000000  1
2  Parametric 0.8117 0.091989  1
3        True 1.9500 0.000000  1
4          GP 1.9500 0.000000  2
5  Parametric 1.1788 0.204966  2
6        True 1.9500 0.000000  2
7          GP 1.9500 0.000000  3
8  Parametric 0.3725 0.143559  3
9        True 1.9500 0.000000  3
10         GP 1.8000 0.000000  4
11 Parametric 0.4041 0.068383  4
12       True 1.9500 0.000000  4
13         GP 1.9500 0.000000  5
14 Parametric 0.1395 0.128373  5
15       True 1.9500 0.000000  5
16         GP 1.7631 0.014326  6
17 Parametric 0.2896 0.085062  6
18       True 1.9500 0.000000  6
19         GP 1.9500 0.000000  7
20 Parametric 0.4041 0.068383  7
21       True 1.9500 0.000000  7
22         GP 1.9500 0.000000  8
23 Parametric 1.3806 0.080586  8
24       True 1.9500 0.000000  8
25         GP 1.8000 0.000000  9
26 Parametric 0.0000 0.000000  9
27       True 1.9500 0.000000  9
28         GP 1.9500 0.000000 10
29 Parametric 0.0000 0.000000 10
30       True 1.9500 0.000000 10
31         GP 1.8988 0.005595 11
32 Parametric 0.2572 0.061829 11
33       True 1.9500 0.000000 11
34         GP 1.9500 0.000000 12
35 Parametric 0.0000 0.000000 12
36       True 1.9500 0.000000 12
37         GP 1.9500 0.000000 13
38 Parametric 0.3138 0.085535 13
39       True 1.9500 0.000000 13
40         GP 1.9500 0.000000 14
41 Parametric 1.3746 0.082639 14
42       True 1.9500 0.000000 14
43         GP 1.7631 0.014326 15
44 Parametric 0.0000 0.000000 15
45       True 1.9500 0.000000 15
46         GP 1.9500 0.000000 16
47 Parametric 0.2558 0.067430 16
48       True 1.9500 0.000000 16
49         GP 1.9500 0.000000 17
50 Parametric 0.1362 0.012719 17
51       True 1.9500 0.000000 17
52         GP 1.9500 0.000000 18
53 Parametric 0.5536 0.086189 18
54       True 1.9500 0.000000 18
55         GP 1.7631 0.014326 19
56 Parametric 0.2558 0.067430 19
57       True 1.9500 0.000000 19
58         GP 1.9500 0.000000 20
59 Parametric 0.1753 0.071906 20
60       True 1.9500 0.000000 20
61         GP 1.9500 0.000000 21
62 Parametric 0.2480 0.047577 21
63       True 1.9500 0.000000 21
64         GP 1.8000 0.000000 22
65 Parametric 0.5179 0.068171 22
66       True 1.9500 0.000000 22
67         GP 1.9500 0.000000 23
68 Parametric 0.2558 0.067430 23
69       True 1.9500 0.000000 23
70         GP 1.9500 0.000000 24
71 Parametric 0.2480 0.047577 24
72       True 1.9500 0.000000 24
73         GP 1.9500 0.000000 25
74 Parametric 0.0000 0.000000 25
75       True 1.9500 0.000000 25
76         GP 1.9500 0.000000 26
77 Parametric 0.0000 0.000000 26
78       True 1.9500 0.000000 26
79         GP 1.9500 0.000000 27
80 Parametric 0.2558 0.067430 27
81       True 1.9500 0.000000 27
82         GP 1.9500 0.000000 28
83 Parametric 0.4041 0.068383 28
84       True 1.9500 0.000000 28
85         GP 1.9500 0.000000 29
86 Parametric 0.1362 0.012719 29
87       True 1.9500 0.000000 29
88         GP 1.9500 0.000000 30
89 Parametric 0.2572 0.061829 30
90       True 1.9500 0.000000 30
91         GP 1.9500 0.000000 31
92 Parametric 0.3516 0.089726 31
93       True 1.9500 0.000000 31
94         GP 1.8988 0.005595 32
95 Parametric 0.3405 0.066929 32
96       True 1.9500 0.000000 32
```




NULL

