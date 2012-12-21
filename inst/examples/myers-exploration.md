









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
 [1]      1 765342 702046 974554 120214 550994 611193 284661 691415 359609
[11] 469373 726434 220686 466199 803742 166768  76519 911809 193433 808768
[21] 641369 177842 154051 866121 743123 800438 912086 139454 173031 285617
[31] 197764 361303
```

```r
yields <- lapply(seed, 
function(seed_i){
  Xo <- allee + x_grid[10] # observations start from
  x0 <- Xo # simulation under policy starts from
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

![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-48-47-3f7bd6aa34-unnamed-chunk-21.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-48-47-3f7bd6aa34-unnamed-chunk-22.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-48-48-3f7bd6aa34-unnamed-chunk-23.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-48-49-3f7bd6aa34-unnamed-chunk-24.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-48-50-3f7bd6aa34-unnamed-chunk-25.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-48-51-3f7bd6aa34-unnamed-chunk-26.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-48-51-3f7bd6aa34-unnamed-chunk-27.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-48-52-3f7bd6aa34-unnamed-chunk-28.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-48-53-3f7bd6aa34-unnamed-chunk-29.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-48-54-3f7bd6aa34-unnamed-chunk-210.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-48-55-3f7bd6aa34-unnamed-chunk-211.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-48-56-3f7bd6aa34-unnamed-chunk-212.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-48-56-3f7bd6aa34-unnamed-chunk-213.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-48-57-3f7bd6aa34-unnamed-chunk-214.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-48-58-3f7bd6aa34-unnamed-chunk-215.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-48-59-3f7bd6aa34-unnamed-chunk-216.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-48-59-3f7bd6aa34-unnamed-chunk-217.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-00-3f7bd6aa34-unnamed-chunk-218.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-01-3f7bd6aa34-unnamed-chunk-219.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-02-3f7bd6aa34-unnamed-chunk-220.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-02-3f7bd6aa34-unnamed-chunk-221.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-03-3f7bd6aa34-unnamed-chunk-222.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-04-3f7bd6aa34-unnamed-chunk-223.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-04-3f7bd6aa34-unnamed-chunk-224.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-05-3f7bd6aa34-unnamed-chunk-225.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-06-3f7bd6aa34-unnamed-chunk-226.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-07-3f7bd6aa34-unnamed-chunk-227.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-07-3f7bd6aa34-unnamed-chunk-228.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-08-3f7bd6aa34-unnamed-chunk-229.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-09-3f7bd6aa34-unnamed-chunk-230.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-10-3f7bd6aa34-unnamed-chunk-231.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-10-3f7bd6aa34-unnamed-chunk-232.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-11-3f7bd6aa34-unnamed-chunk-233.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-12-3f7bd6aa34-unnamed-chunk-234.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-12-3f7bd6aa34-unnamed-chunk-235.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-13-3f7bd6aa34-unnamed-chunk-236.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-14-3f7bd6aa34-unnamed-chunk-237.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-15-3f7bd6aa34-unnamed-chunk-238.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-15-3f7bd6aa34-unnamed-chunk-239.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-21-3f7bd6aa34-unnamed-chunk-240.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-22-3f7bd6aa34-unnamed-chunk-241.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-23-3f7bd6aa34-unnamed-chunk-242.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-23-3f7bd6aa34-unnamed-chunk-243.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-24-3f7bd6aa34-unnamed-chunk-244.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-25-3f7bd6aa34-unnamed-chunk-245.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-25-3f7bd6aa34-unnamed-chunk-246.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-26-3f7bd6aa34-unnamed-chunk-247.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-27-3f7bd6aa34-unnamed-chunk-248.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-28-3f7bd6aa34-unnamed-chunk-249.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-28-3f7bd6aa34-unnamed-chunk-250.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-29-3f7bd6aa34-unnamed-chunk-251.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-30-3f7bd6aa34-unnamed-chunk-252.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-31-3f7bd6aa34-unnamed-chunk-253.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-36-3f7bd6aa34-unnamed-chunk-254.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-37-3f7bd6aa34-unnamed-chunk-255.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-38-3f7bd6aa34-unnamed-chunk-256.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-39-3f7bd6aa34-unnamed-chunk-257.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-39-3f7bd6aa34-unnamed-chunk-258.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-40-3f7bd6aa34-unnamed-chunk-259.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-41-3f7bd6aa34-unnamed-chunk-260.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-42-3f7bd6aa34-unnamed-chunk-261.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-42-3f7bd6aa34-unnamed-chunk-262.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-43-3f7bd6aa34-unnamed-chunk-263.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-44-3f7bd6aa34-unnamed-chunk-264.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-45-3f7bd6aa34-unnamed-chunk-265.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-45-3f7bd6aa34-unnamed-chunk-266.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-46-3f7bd6aa34-unnamed-chunk-267.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-47-3f7bd6aa34-unnamed-chunk-268.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-47-3f7bd6aa34-unnamed-chunk-269.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-48-3f7bd6aa34-unnamed-chunk-270.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-49-3f7bd6aa34-unnamed-chunk-271.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-50-3f7bd6aa34-unnamed-chunk-272.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-50-3f7bd6aa34-unnamed-chunk-273.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-51-3f7bd6aa34-unnamed-chunk-274.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-52-3f7bd6aa34-unnamed-chunk-275.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-53-3f7bd6aa34-unnamed-chunk-276.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-53-3f7bd6aa34-unnamed-chunk-277.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-54-3f7bd6aa34-unnamed-chunk-278.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-55-3f7bd6aa34-unnamed-chunk-279.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-56-3f7bd6aa34-unnamed-chunk-280.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-56-3f7bd6aa34-unnamed-chunk-281.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-57-3f7bd6aa34-unnamed-chunk-282.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-58-3f7bd6aa34-unnamed-chunk-283.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-59-3f7bd6aa34-unnamed-chunk-284.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-49-59-3f7bd6aa34-unnamed-chunk-285.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-00-3f7bd6aa34-unnamed-chunk-286.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-01-3f7bd6aa34-unnamed-chunk-287.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-07-3f7bd6aa34-unnamed-chunk-288.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-07-3f7bd6aa34-unnamed-chunk-289.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-08-3f7bd6aa34-unnamed-chunk-290.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-09-3f7bd6aa34-unnamed-chunk-291.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-09-3f7bd6aa34-unnamed-chunk-292.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-10-3f7bd6aa34-unnamed-chunk-293.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-11-3f7bd6aa34-unnamed-chunk-294.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-12-3f7bd6aa34-unnamed-chunk-295.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-12-3f7bd6aa34-unnamed-chunk-296.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-13-3f7bd6aa34-unnamed-chunk-297.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-14-3f7bd6aa34-unnamed-chunk-298.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-14-3f7bd6aa34-unnamed-chunk-299.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-15-3f7bd6aa34-unnamed-chunk-2100.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-16-3f7bd6aa34-unnamed-chunk-2101.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-22-3f7bd6aa34-unnamed-chunk-2102.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-22-3f7bd6aa34-unnamed-chunk-2103.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-23-3f7bd6aa34-unnamed-chunk-2104.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-24-3f7bd6aa34-unnamed-chunk-2105.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-25-3f7bd6aa34-unnamed-chunk-2106.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-25-3f7bd6aa34-unnamed-chunk-2107.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-26-3f7bd6aa34-unnamed-chunk-2108.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-27-3f7bd6aa34-unnamed-chunk-2109.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-28-3f7bd6aa34-unnamed-chunk-2110.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-28-3f7bd6aa34-unnamed-chunk-2111.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-29-3f7bd6aa34-unnamed-chunk-2112.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-30-3f7bd6aa34-unnamed-chunk-2113.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-31-3f7bd6aa34-unnamed-chunk-2114.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-31-3f7bd6aa34-unnamed-chunk-2115.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-32-3f7bd6aa34-unnamed-chunk-2116.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-33-3f7bd6aa34-unnamed-chunk-2117.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-33-3f7bd6aa34-unnamed-chunk-2118.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-34-3f7bd6aa34-unnamed-chunk-2119.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-35-3f7bd6aa34-unnamed-chunk-2120.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-36-3f7bd6aa34-unnamed-chunk-2121.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-36-3f7bd6aa34-unnamed-chunk-2122.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-37-3f7bd6aa34-unnamed-chunk-2123.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-38-3f7bd6aa34-unnamed-chunk-2124.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-39-3f7bd6aa34-unnamed-chunk-2125.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-39-3f7bd6aa34-unnamed-chunk-2126.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-40-3f7bd6aa34-unnamed-chunk-2127.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-20-50-41-3f7bd6aa34-unnamed-chunk-2128.png) 



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

![plot of chunk unnamed-chunk-3](http://carlboettiger.info/assets/figures/2012-12-20-20-50-42-3f7bd6aa34-unnamed-chunk-3.png) 

```r
yields
```

```
       method      V1      sd L1
1          GP 11.4345 1.63844  1
2  Parametric  3.9598 1.67446  1
3        True 11.6025 1.78496  1
4          GP 11.7585 1.70572  2
5  Parametric  2.6915 0.66033  2
6        True 11.6025 1.78496  2
7          GP 11.4330 1.63576  3
8  Parametric  2.4028 0.07108  3
9        True 11.6025 1.78496  3
10         GP 11.1480 1.60990  4
11 Parametric  2.0145 0.14442  4
12       True 11.6025 1.78496  4
13         GP 11.7630 1.70976  5
14 Parametric  2.7404 0.55265  5
15       True 11.6025 1.78496  5
16         GP 11.6505 1.68524  6
17 Parametric  2.5231 0.24360  6
18       True 11.6025 1.78496  6
19         GP 11.0190 1.56632  7
20 Parametric  0.5235 0.33541  7
21       True 11.6025 1.78496  7
22         GP  9.7980 1.48337  8
23 Parametric  0.7995 0.40286  8
24       True 11.6025 1.78496  8
25         GP 11.1975 1.59379  9
26 Parametric  2.2595 0.32122  9
27       True 11.6025 1.78496  9
28         GP 11.4345 1.63844 10
29 Parametric  2.4175 0.73484 10
30       True 11.6025 1.78496 10
31         GP 10.3275 1.52589 11
32 Parametric  2.2045 0.38366 11
33       True 11.6025 1.78496 11
34         GP 11.6505 1.68834 12
35 Parametric  2.2815 0.31068 12
36       True 11.6025 1.78496 12
37         GP 11.2095 1.60154 13
38 Parametric  0.5235 0.33541 13
39       True 11.6025 1.78496 13
40         GP 11.7615 1.71014 14
41 Parametric  2.6511 0.25403 14
42       True 11.6025 1.78496 14
43         GP 11.1480 1.60990 15
44 Parametric  2.3789 0.25189 15
45       True 11.6025 1.78496 15
46         GP  9.9105 1.46741 16
47 Parametric  0.5235 0.33541 16
48       True 11.6025 1.78496 16
49         GP 11.6625 1.68326 17
50 Parametric  2.4028 0.07108 17
51       True 11.6025 1.78496 17
52         GP 11.6760 1.69702 18
53 Parametric  2.7380 0.55549 18
54       True 11.6025 1.78496 18
55         GP 11.7000 1.70080 19
56 Parametric  2.5238 0.20871 19
57       True 11.6025 1.78496 19
58         GP 11.1480 1.60990 20
59 Parametric  2.4768 0.21650 20
60       True 11.6025 1.78496 20
61         GP 11.7555 1.70383 21
62 Parametric  2.4598 0.43125 21
63       True 11.6025 1.78496 21
64         GP  6.2085 1.41162 22
65 Parametric  0.5235 0.33541 22
66       True 11.6025 1.78496 22
67         GP 11.6760 1.69702 23
68 Parametric  2.4975 0.06672 23
69       True 11.6025 1.78496 23
70         GP 11.6460 1.67760 24
71 Parametric  2.2356 0.30879 24
72       True 11.6025 1.78496 24
73         GP 11.4705 1.65220 25
74 Parametric  2.5807 0.24771 25
75       True 11.6025 1.78496 25
76         GP 11.1480 1.60990 26
77 Parametric  2.6348 0.10825 26
78       True 11.6025 1.78496 26
79         GP  9.8310 1.46960 27
80 Parametric  0.5235 0.33541 27
81       True 11.6025 1.78496 27
82         GP 11.4405 1.65074 28
83 Parametric 11.7030 1.77981 28
84       True 11.6025 1.78496 28
85         GP 11.4345 1.63844 29
86 Parametric  2.2231 0.12769 29
87       True 11.6025 1.78496 29
88         GP 11.7615 1.71014 30
89 Parametric  2.5517 0.54640 30
90       True 11.6025 1.78496 30
91         GP 10.8135 1.57308 31
92 Parametric  0.5235 0.33541 31
93       True 11.6025 1.78496 31
94         GP 11.6505 1.68834 32
95 Parametric  0.5235 0.33541 32
96       True 11.6025 1.78496 32
```




NULL

