









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
f <- Myers
p <- c(1, 2, 6) 
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


With parameters 1, 2, 6. 


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
 [1]      1 564401 397742 749419 447478 452644 693258 541466 639686 383114
[11] 630719 120662 652032 502257 418747 396757 276956 222107 408188 404791
[21] 836402 983857  31590 738651 567066 660329 271326 712728 568539 279225
[31]  40863 149409
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

![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-01-8e5a4cd1ab-unnamed-chunk-21.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-02-8e5a4cd1ab-unnamed-chunk-22.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-03-8e5a4cd1ab-unnamed-chunk-23.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-04-8e5a4cd1ab-unnamed-chunk-24.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-04-8e5a4cd1ab-unnamed-chunk-25.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-05-8e5a4cd1ab-unnamed-chunk-26.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-06-8e5a4cd1ab-unnamed-chunk-27.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-07-8e5a4cd1ab-unnamed-chunk-28.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-07-8e5a4cd1ab-unnamed-chunk-29.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-08-8e5a4cd1ab-unnamed-chunk-210.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-09-8e5a4cd1ab-unnamed-chunk-211.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-10-8e5a4cd1ab-unnamed-chunk-212.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-10-8e5a4cd1ab-unnamed-chunk-213.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-11-8e5a4cd1ab-unnamed-chunk-214.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-12-8e5a4cd1ab-unnamed-chunk-215.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-13-8e5a4cd1ab-unnamed-chunk-216.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-13-8e5a4cd1ab-unnamed-chunk-217.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-19-8e5a4cd1ab-unnamed-chunk-218.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-20-8e5a4cd1ab-unnamed-chunk-219.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-21-8e5a4cd1ab-unnamed-chunk-220.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-21-8e5a4cd1ab-unnamed-chunk-221.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-22-8e5a4cd1ab-unnamed-chunk-222.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-23-8e5a4cd1ab-unnamed-chunk-223.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-24-8e5a4cd1ab-unnamed-chunk-224.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-24-8e5a4cd1ab-unnamed-chunk-225.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-25-8e5a4cd1ab-unnamed-chunk-226.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-26-8e5a4cd1ab-unnamed-chunk-227.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-27-8e5a4cd1ab-unnamed-chunk-228.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-27-8e5a4cd1ab-unnamed-chunk-229.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-28-8e5a4cd1ab-unnamed-chunk-230.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-29-8e5a4cd1ab-unnamed-chunk-231.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-30-8e5a4cd1ab-unnamed-chunk-232.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-30-8e5a4cd1ab-unnamed-chunk-233.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-31-8e5a4cd1ab-unnamed-chunk-234.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-32-8e5a4cd1ab-unnamed-chunk-235.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-33-8e5a4cd1ab-unnamed-chunk-236.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-33-8e5a4cd1ab-unnamed-chunk-237.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-34-8e5a4cd1ab-unnamed-chunk-238.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-35-8e5a4cd1ab-unnamed-chunk-239.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-36-8e5a4cd1ab-unnamed-chunk-240.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-36-8e5a4cd1ab-unnamed-chunk-241.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-37-8e5a4cd1ab-unnamed-chunk-242.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-38-8e5a4cd1ab-unnamed-chunk-243.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-39-8e5a4cd1ab-unnamed-chunk-244.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-39-8e5a4cd1ab-unnamed-chunk-245.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-40-8e5a4cd1ab-unnamed-chunk-246.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-41-8e5a4cd1ab-unnamed-chunk-247.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-42-8e5a4cd1ab-unnamed-chunk-248.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-47-8e5a4cd1ab-unnamed-chunk-249.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-48-8e5a4cd1ab-unnamed-chunk-250.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-49-8e5a4cd1ab-unnamed-chunk-251.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-50-8e5a4cd1ab-unnamed-chunk-252.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-50-8e5a4cd1ab-unnamed-chunk-253.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-51-8e5a4cd1ab-unnamed-chunk-254.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-52-8e5a4cd1ab-unnamed-chunk-255.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-53-8e5a4cd1ab-unnamed-chunk-256.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-53-8e5a4cd1ab-unnamed-chunk-257.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-54-8e5a4cd1ab-unnamed-chunk-258.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-55-8e5a4cd1ab-unnamed-chunk-259.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-56-8e5a4cd1ab-unnamed-chunk-260.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-56-8e5a4cd1ab-unnamed-chunk-261.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-57-8e5a4cd1ab-unnamed-chunk-262.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-51-58-8e5a4cd1ab-unnamed-chunk-263.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-00-8e5a4cd1ab-unnamed-chunk-264.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-00-8e5a4cd1ab-unnamed-chunk-265.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-01-8e5a4cd1ab-unnamed-chunk-266.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-02-8e5a4cd1ab-unnamed-chunk-267.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-02-8e5a4cd1ab-unnamed-chunk-268.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-03-8e5a4cd1ab-unnamed-chunk-269.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-04-8e5a4cd1ab-unnamed-chunk-270.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-05-8e5a4cd1ab-unnamed-chunk-271.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-05-8e5a4cd1ab-unnamed-chunk-272.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-06-8e5a4cd1ab-unnamed-chunk-273.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-07-8e5a4cd1ab-unnamed-chunk-274.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-08-8e5a4cd1ab-unnamed-chunk-275.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-08-8e5a4cd1ab-unnamed-chunk-276.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-09-8e5a4cd1ab-unnamed-chunk-277.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-10-8e5a4cd1ab-unnamed-chunk-278.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-11-8e5a4cd1ab-unnamed-chunk-279.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-11-8e5a4cd1ab-unnamed-chunk-280.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-12-8e5a4cd1ab-unnamed-chunk-281.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-13-8e5a4cd1ab-unnamed-chunk-282.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-14-8e5a4cd1ab-unnamed-chunk-283.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-19-8e5a4cd1ab-unnamed-chunk-284.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-20-8e5a4cd1ab-unnamed-chunk-285.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-21-8e5a4cd1ab-unnamed-chunk-286.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-22-8e5a4cd1ab-unnamed-chunk-287.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-22-8e5a4cd1ab-unnamed-chunk-288.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-23-8e5a4cd1ab-unnamed-chunk-289.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-24-8e5a4cd1ab-unnamed-chunk-290.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-25-8e5a4cd1ab-unnamed-chunk-291.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-25-8e5a4cd1ab-unnamed-chunk-292.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-26-8e5a4cd1ab-unnamed-chunk-293.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-27-8e5a4cd1ab-unnamed-chunk-294.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-28-8e5a4cd1ab-unnamed-chunk-295.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-28-8e5a4cd1ab-unnamed-chunk-296.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-34-8e5a4cd1ab-unnamed-chunk-297.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-35-8e5a4cd1ab-unnamed-chunk-298.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-36-8e5a4cd1ab-unnamed-chunk-299.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-36-8e5a4cd1ab-unnamed-chunk-2100.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-37-8e5a4cd1ab-unnamed-chunk-2101.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-38-8e5a4cd1ab-unnamed-chunk-2102.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-39-8e5a4cd1ab-unnamed-chunk-2103.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-39-8e5a4cd1ab-unnamed-chunk-2104.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-40-8e5a4cd1ab-unnamed-chunk-2105.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-41-8e5a4cd1ab-unnamed-chunk-2106.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-42-8e5a4cd1ab-unnamed-chunk-2107.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-42-8e5a4cd1ab-unnamed-chunk-2108.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-43-8e5a4cd1ab-unnamed-chunk-2109.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-44-8e5a4cd1ab-unnamed-chunk-2110.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-45-8e5a4cd1ab-unnamed-chunk-2111.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-45-8e5a4cd1ab-unnamed-chunk-2112.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-46-8e5a4cd1ab-unnamed-chunk-2113.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-47-8e5a4cd1ab-unnamed-chunk-2114.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-48-8e5a4cd1ab-unnamed-chunk-2115.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-48-8e5a4cd1ab-unnamed-chunk-2116.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-49-8e5a4cd1ab-unnamed-chunk-2117.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-50-8e5a4cd1ab-unnamed-chunk-2118.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-51-8e5a4cd1ab-unnamed-chunk-2119.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-52-8e5a4cd1ab-unnamed-chunk-2120.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-52-8e5a4cd1ab-unnamed-chunk-2121.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-53-8e5a4cd1ab-unnamed-chunk-2122.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-54-8e5a4cd1ab-unnamed-chunk-2123.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-55-8e5a4cd1ab-unnamed-chunk-2124.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-55-8e5a4cd1ab-unnamed-chunk-2125.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-56-8e5a4cd1ab-unnamed-chunk-2126.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-57-8e5a4cd1ab-unnamed-chunk-2127.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-19-52-58-8e5a4cd1ab-unnamed-chunk-2128.png) 



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

![plot of chunk unnamed-chunk-3](http://carlboettiger.info/assets/figures/2012-12-20-19-53-04-8e5a4cd1ab-unnamed-chunk-3.png) 

```r
yields
```

```
       method      V1      sd L1
1          GP 9.01740 0.57858  1
2  Parametric 0.00000 0.00000  1
3        True 9.69171 0.70437  1
4          GP 8.74767 0.55809  2
5  Parametric 8.62061 0.55830  2
6        True 9.69171 0.70437  2
7          GP 9.40566 0.61843  3
8  Parametric 8.27849 0.54308  3
9        True 9.69171 0.70437  3
10         GP 9.55046 0.62870  4
11 Parametric 9.64487 0.70904  4
12       True 9.69171 0.70437  4
13         GP 7.73122 0.49750  5
14 Parametric 0.00000 0.00000  5
15       True 9.69171 0.70437  5
16         GP 9.42412 0.62059  6
17 Parametric 8.70011 0.56084  6
18       True 9.69171 0.70437  6
19         GP 9.42483 0.61406  7
20 Parametric 9.02733 0.61105  7
21       True 9.69171 0.70437  7
22         GP 9.55969 0.62979  8
23 Parametric 8.69301 0.56420  8
24       True 9.69171 0.70437  8
25         GP 9.40566 0.61843  9
26 Parametric 0.00000 0.00000  9
27       True 9.69171 0.70437  9
28         GP 8.74696 0.55697 10
29 Parametric 8.71715 0.60966 10
30       True 9.69171 0.70437 10
31         GP 9.15652 0.58253 11
32 Parametric 8.62061 0.55830 11
33       True 9.69171 0.70437 11
34         GP 9.73217 0.66202 12
35 Parametric 8.44742 0.55511 12
36       True 9.69171 0.70437 12
37         GP 9.41773 0.61796 13
38 Parametric 9.20904 0.60462 13
39       True 9.69171 0.70437 13
40         GP 8.55673 0.54230 14
41 Parametric 0.00000 0.00000 14
42       True 9.69171 0.70437 14
43         GP 8.98829 0.58084 15
44 Parametric 9.37656 0.72369 15
45       True 9.69171 0.70437 15
46         GP 9.66119 0.64890 16
47 Parametric 0.00000 0.00000 16
48       True 9.69171 0.70437 16
49         GP 8.77251 0.55681 17
50 Parametric 0.03123 0.04658 17
51       True 9.69171 0.70437 17
52         GP 9.00675 0.57990 18
53 Parametric 9.37656 0.72369 18
54       True 9.69171 0.70437 18
55         GP 9.67539 0.65053 19
56 Parametric 7.98817 0.59202 19
57       True 9.69171 0.70437 19
58         GP 8.55602 0.54193 20
59 Parametric 9.70875 0.66987 20
60       True 9.69171 0.70437 20
61         GP 9.57318 0.63564 21
62 Parametric 0.02484 0.04772 21
63       True 9.69171 0.70437 21
64         GP 8.98404 0.57817 22
65 Parametric 7.59849 0.50579 22
66       True 9.69171 0.70437 22
67         GP 8.77038 0.55509 23
68 Parametric 0.00000 0.00000 23
69       True 9.69171 0.70437 23
70         GP 9.23176 0.60236 24
71 Parametric 6.07808 0.44111 24
72       True 9.69171 0.70437 24
73         GP 9.55046 0.62870 25
74 Parametric 8.37502 0.53996 25
75       True 9.69171 0.70437 25
76         GP 9.42199 0.62588 26
77 Parametric 8.27849 0.54308 26
78       True 9.69171 0.70437 26
79         GP 7.80221 0.50127 27
80 Parametric 9.20904 0.60462 27
81       True 9.69171 0.70437 27
82         GP 9.22537 0.59792 28
83 Parametric 0.00000 0.00000 28
84       True 9.69171 0.70437 28
85         GP 9.55046 0.62870 29
86 Parametric 9.20904 0.60462 29
87       True 9.69171 0.70437 29
88         GP 7.19248 0.47888 30
89 Parametric 9.43121 0.62025 30
90       True 9.69171 0.70437 30
91         GP 8.89886 0.58381 31
92 Parametric 0.00000 0.00000 31
93       True 9.69171 0.70437 31
94         GP 9.11677 0.60358 32
95 Parametric 0.00000 0.00000 32
96       True 9.69171 0.70437 32
```




NULL

