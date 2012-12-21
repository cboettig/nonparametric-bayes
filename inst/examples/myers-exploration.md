










Fixed priors on hyperparameters, fixed model type.


```r
s2.p <- c(50,50)  # inverse gamma mean is 
tau2.p <- c(20,1)
d.p = c(10, 1/0.01, 10, 1/0.01)
nug.p = c(10, 1/0.01, 10, 1/0.01) # gamma mean
s2_prior <- function(x) dinvgamma(x, s2.p[1], s2.p[2])
tau2_prior <- function(x) dinvgamma(x, tau2.p[1], tau2.p[2])
d_prior <- function(x) dgamma(x, d.p[1], scale = d.p[2]) + dgamma(x, d.p[3], scale = d.p[4])
nug_prior <- function(x) dgamma(x, nug.p[1], scale = nug.p[2]) + dgamma(x, nug.p[3], scale = nug.p[4])
beta0_prior <- function(x, tau) dnorm(x, 0, tau)
beta = c(0)
priors <- list(s2 = s2_prior, tau2 = tau2_prior, beta0 = dnorm, nug = nug_prior, d = d_prior, ldetK = function(x) 0)
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
sigma_g <- 0.01
sigma_m <- 0.02
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



posteriors_plot <- function(gp, priors){
  hyperparameters <- c("index", "s2", "tau2", "beta0", "nug", "d", "ldetK")
  posteriors <- melt(gp$trace$XX[[1]][,hyperparameters], id="index")
  prior_curves <- ddply(posteriors, "variable", function(dd){
    grid <- seq(min(dd$value), max(dd$value), length = 100)
    data.frame(value = grid, density = priors[[dd$variable[1]]](grid))
  })
  plot_posteriors <- ggplot(posteriors) + 
    #geom_density(aes(value), lwd=2) +
    geom_histogram(aes(x=value, y=..density..), alpha=0.7) +
    geom_line(data=prior_curves, aes(x=value, y=density), col="red", lwd=2) +
    facet_wrap(~ variable, scale="free")
  print(plot_posteriors)
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
 [1]      1 749185 790609 845332 785431 174721 732117 755546 856947 314670
[11] 591214  24567 805565 886432 104630 309463 259650 329741 793992  97328
[21] 146720 860803 411101 855389 306821 516549 418233 536502 573028 797686
[31] 879221 195760
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
#  posteriors_plot(gp, priors) # needs trace=TRUE!
  OPT <- optimal_policy(gp, f, est$f_alt, p, est$p_alt, x_grid, h_grid, sigma_g, est$sigma_g_alt, delta, xT, profit, reward, OptTime)
  dt <- simulate_opt(OPT, f, p, x_grid, h_grid, x0, z_g, profit)
  sim_plots(dt, seed=seed_i)
  profits_stats(dt)
})
```

![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-12-96ad28dd6e-unnamed-chunk-21.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-12-96ad28dd6e-unnamed-chunk-22.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-13-96ad28dd6e-unnamed-chunk-23.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-14-96ad28dd6e-unnamed-chunk-24.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-15-96ad28dd6e-unnamed-chunk-25.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-16-96ad28dd6e-unnamed-chunk-26.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-17-96ad28dd6e-unnamed-chunk-27.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-17-96ad28dd6e-unnamed-chunk-28.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-18-96ad28dd6e-unnamed-chunk-29.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-19-96ad28dd6e-unnamed-chunk-210.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-20-96ad28dd6e-unnamed-chunk-211.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-20-96ad28dd6e-unnamed-chunk-212.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-21-96ad28dd6e-unnamed-chunk-213.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-22-96ad28dd6e-unnamed-chunk-214.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-22-96ad28dd6e-unnamed-chunk-215.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-23-96ad28dd6e-unnamed-chunk-216.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-24-96ad28dd6e-unnamed-chunk-217.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-25-96ad28dd6e-unnamed-chunk-218.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-25-96ad28dd6e-unnamed-chunk-219.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-26-96ad28dd6e-unnamed-chunk-220.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-27-96ad28dd6e-unnamed-chunk-221.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-27-96ad28dd6e-unnamed-chunk-222.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-28-96ad28dd6e-unnamed-chunk-223.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-29-96ad28dd6e-unnamed-chunk-224.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-30-96ad28dd6e-unnamed-chunk-225.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-30-96ad28dd6e-unnamed-chunk-226.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-31-96ad28dd6e-unnamed-chunk-227.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-32-96ad28dd6e-unnamed-chunk-228.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-33-96ad28dd6e-unnamed-chunk-229.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-33-96ad28dd6e-unnamed-chunk-230.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-34-96ad28dd6e-unnamed-chunk-231.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-35-96ad28dd6e-unnamed-chunk-232.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-35-96ad28dd6e-unnamed-chunk-233.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-36-96ad28dd6e-unnamed-chunk-234.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-37-96ad28dd6e-unnamed-chunk-235.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-38-96ad28dd6e-unnamed-chunk-236.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-38-96ad28dd6e-unnamed-chunk-237.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-39-96ad28dd6e-unnamed-chunk-238.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-40-96ad28dd6e-unnamed-chunk-239.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-40-96ad28dd6e-unnamed-chunk-240.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-41-96ad28dd6e-unnamed-chunk-241.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-42-96ad28dd6e-unnamed-chunk-242.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-43-96ad28dd6e-unnamed-chunk-243.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-43-96ad28dd6e-unnamed-chunk-244.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-44-96ad28dd6e-unnamed-chunk-245.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-45-96ad28dd6e-unnamed-chunk-246.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-45-96ad28dd6e-unnamed-chunk-247.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-46-96ad28dd6e-unnamed-chunk-248.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-47-96ad28dd6e-unnamed-chunk-249.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-47-96ad28dd6e-unnamed-chunk-250.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-48-96ad28dd6e-unnamed-chunk-251.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-49-96ad28dd6e-unnamed-chunk-252.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-50-96ad28dd6e-unnamed-chunk-253.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-50-96ad28dd6e-unnamed-chunk-254.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-51-96ad28dd6e-unnamed-chunk-255.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-52-96ad28dd6e-unnamed-chunk-256.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-53-96ad28dd6e-unnamed-chunk-257.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-53-96ad28dd6e-unnamed-chunk-258.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-54-96ad28dd6e-unnamed-chunk-259.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-55-96ad28dd6e-unnamed-chunk-260.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-55-96ad28dd6e-unnamed-chunk-261.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-56-96ad28dd6e-unnamed-chunk-262.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-57-96ad28dd6e-unnamed-chunk-263.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-58-96ad28dd6e-unnamed-chunk-264.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-58-96ad28dd6e-unnamed-chunk-265.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-56-59-96ad28dd6e-unnamed-chunk-266.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-00-96ad28dd6e-unnamed-chunk-267.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-01-96ad28dd6e-unnamed-chunk-268.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-01-96ad28dd6e-unnamed-chunk-269.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-02-96ad28dd6e-unnamed-chunk-270.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-03-96ad28dd6e-unnamed-chunk-271.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-03-96ad28dd6e-unnamed-chunk-272.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-04-96ad28dd6e-unnamed-chunk-273.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-05-96ad28dd6e-unnamed-chunk-274.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-11-96ad28dd6e-unnamed-chunk-275.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-11-96ad28dd6e-unnamed-chunk-276.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-12-96ad28dd6e-unnamed-chunk-277.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-13-96ad28dd6e-unnamed-chunk-278.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-13-96ad28dd6e-unnamed-chunk-279.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-14-96ad28dd6e-unnamed-chunk-280.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-15-96ad28dd6e-unnamed-chunk-281.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-16-96ad28dd6e-unnamed-chunk-282.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-16-96ad28dd6e-unnamed-chunk-283.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-17-96ad28dd6e-unnamed-chunk-284.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-18-96ad28dd6e-unnamed-chunk-285.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-18-96ad28dd6e-unnamed-chunk-286.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-19-96ad28dd6e-unnamed-chunk-287.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-20-96ad28dd6e-unnamed-chunk-288.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-21-96ad28dd6e-unnamed-chunk-289.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-21-96ad28dd6e-unnamed-chunk-290.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-22-96ad28dd6e-unnamed-chunk-291.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-23-96ad28dd6e-unnamed-chunk-292.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-23-96ad28dd6e-unnamed-chunk-293.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-24-96ad28dd6e-unnamed-chunk-294.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-25-96ad28dd6e-unnamed-chunk-295.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-26-96ad28dd6e-unnamed-chunk-296.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-26-96ad28dd6e-unnamed-chunk-297.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-27-96ad28dd6e-unnamed-chunk-298.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-28-96ad28dd6e-unnamed-chunk-299.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-29-96ad28dd6e-unnamed-chunk-2100.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-29-96ad28dd6e-unnamed-chunk-2101.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-30-96ad28dd6e-unnamed-chunk-2102.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-31-96ad28dd6e-unnamed-chunk-2103.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-32-96ad28dd6e-unnamed-chunk-2104.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-32-96ad28dd6e-unnamed-chunk-2105.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-33-96ad28dd6e-unnamed-chunk-2106.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-34-96ad28dd6e-unnamed-chunk-2107.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-34-96ad28dd6e-unnamed-chunk-2108.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-35-96ad28dd6e-unnamed-chunk-2109.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-41-96ad28dd6e-unnamed-chunk-2110.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-42-96ad28dd6e-unnamed-chunk-2111.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-42-96ad28dd6e-unnamed-chunk-2112.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-43-96ad28dd6e-unnamed-chunk-2113.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-44-96ad28dd6e-unnamed-chunk-2114.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-44-96ad28dd6e-unnamed-chunk-2115.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-45-96ad28dd6e-unnamed-chunk-2116.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-46-96ad28dd6e-unnamed-chunk-2117.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-47-96ad28dd6e-unnamed-chunk-2118.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-47-96ad28dd6e-unnamed-chunk-2119.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-48-96ad28dd6e-unnamed-chunk-2120.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-49-96ad28dd6e-unnamed-chunk-2121.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-50-96ad28dd6e-unnamed-chunk-2122.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-50-96ad28dd6e-unnamed-chunk-2123.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-51-96ad28dd6e-unnamed-chunk-2124.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-52-96ad28dd6e-unnamed-chunk-2125.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-52-96ad28dd6e-unnamed-chunk-2126.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-53-96ad28dd6e-unnamed-chunk-2127.png) ![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-20-22-57-54-96ad28dd6e-unnamed-chunk-2128.png) 



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

![plot of chunk unnamed-chunk-3](http://carlboettiger.info/assets/figures/2012-12-20-22-57-55-96ad28dd6e-unnamed-chunk-3.png) 

```r
yields
```

```
       method     V1      sd L1
1          GP 11.409 0.33698  1
2  Parametric  0.219 0.13386  1
3        True 11.695 0.34867  1
4          GP 11.106 0.32747  2
5  Parametric  2.507 0.20706  2
6        True 11.695 0.34867  2
7          GP 11.106 0.32747  3
8  Parametric  0.000 0.00000  3
9        True 11.695 0.34867  3
10         GP 11.409 0.33698  4
11 Parametric  0.219 0.13386  4
12       True 11.695 0.34867  4
13         GP 11.106 0.32747  5
14 Parametric  0.000 0.00000  5
15       True 11.695 0.34867  5
16         GP 11.418 0.33661  6
17 Parametric  0.219 0.13386  6
18       True 11.695 0.34867  6
19         GP 11.106 0.32747  7
20 Parametric  2.820 0.33209  7
21       True 11.695 0.34867  7
22         GP 11.409 0.33698  8
23 Parametric  0.219 0.13386  8
24       True 11.695 0.34867  8
25         GP 11.416 0.33432  9
26 Parametric  0.219 0.13386  9
27       True 11.695 0.34867  9
28         GP 11.409 0.33698 10
29 Parametric  0.000 0.00000 10
30       True 11.695 0.34867 10
31         GP 11.117 0.31615 11
32 Parametric  2.945 0.12709 11
33       True 11.695 0.34867 11
34         GP 11.409 0.33698 12
35 Parametric  0.219 0.13386 12
36       True 11.695 0.34867 12
37         GP 11.106 0.32747 13
38 Parametric  0.219 0.13386 13
39       True 11.695 0.34867 13
40         GP 11.416 0.33432 14
41 Parametric  0.000 0.00000 14
42       True 11.695 0.34867 14
43         GP 11.418 0.33661 15
44 Parametric  0.219 0.13386 15
45       True 11.695 0.34867 15
46         GP 11.106 0.32747 16
47 Parametric  2.882 0.25726 16
48       True 11.695 0.34867 16
49         GP 11.236 0.34384 17
50 Parametric  0.219 0.13386 17
51       True 11.695 0.34867 17
52         GP 11.418 0.33661 18
53 Parametric  2.341 0.06012 18
54       True 11.695 0.34867 18
55         GP 11.409 0.33698 19
56 Parametric  2.455 0.20387 19
57       True 11.695 0.34867 19
58         GP 11.614 0.33252 20
59 Parametric  0.219 0.13386 20
60       True 11.695 0.34867 20
61         GP 11.107 0.32709 21
62 Parametric  0.219 0.13386 21
63       True 11.695 0.34867 21
64         GP 11.409 0.33698 22
65 Parametric  2.284 0.19695 22
66       True 11.695 0.34867 22
67         GP 11.107 0.32709 23
68 Parametric  2.172 0.06204 23
69       True 11.695 0.34867 23
70         GP 11.416 0.33432 24
71 Parametric  0.219 0.13386 24
72       True 11.695 0.34867 24
73         GP 11.106 0.32747 25
74 Parametric  2.381 0.32609 25
75       True 11.695 0.34867 25
76         GP 11.418 0.33661 26
77 Parametric  2.392 0.11605 26
78       True 11.695 0.34867 26
79         GP 11.236 0.34384 27
80 Parametric  1.806 0.19317 27
81       True 11.695 0.34867 27
82         GP 11.409 0.33698 28
83 Parametric  0.219 0.13386 28
84       True 11.695 0.34867 28
85         GP 11.106 0.32747 29
86 Parametric  5.413 0.25504 29
87       True 11.695 0.34867 29
88         GP 11.614 0.33252 30
89 Parametric  2.505 0.12871 30
90       True 11.695 0.34867 30
91         GP 11.106 0.32747 31
92 Parametric  0.000 0.00000 31
93       True 11.695 0.34867 31
94         GP 11.418 0.33661 32
95 Parametric  0.000 0.00000 32
96       True 11.695 0.34867 32
```




NULL

