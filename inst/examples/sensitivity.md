










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
p <- c(1.1, 10, 5) 
K <- 10
allee <- 5
```




```r
sigma_g <- 0.1
sigma_m <- 0.0
x_grid <- seq(0, 1.5 * K, length=101)
h_grid <- x_grid
```


With parameters 1.1, 10, 5. 



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
 [1]      1 727529 331616 420763 650304 714859 715363 980890 392119 502935
[11] 163612 764559 550056 968107 108920 294786 310274  53317 463928 965482
[21] 496211 464446 462772 906454 421686 423394  59475 395987 551138 549884
[31]  47510 140824
```

```r
yields <- lapply(seed, 
function(seed_i){
  Xo <- K # observations start from
  x0 <- Xo # simulation under policy starts from
  obs <- sim_obs(Xo, z_g, f, p, Tobs=35, nz= 1, 
                 harvest = sort(rep(seq(0, .8, length=7), 5)), seed = seed_i)
  alt <- par_est(obs)
  est <- par_est_allee(obs, f, p)
  gp <- bgp(X=obs$x, XX=x_grid, Z=obs$y, verb=0,
          meanfn="constant", bprior="b0", BTE=c(2000,16000,2),
          m0r1=FALSE, corr="exp", trace=FALSE, 
          beta = beta, s2.p = s2.p, d.p = d.p, nug.p = nug.p, tau2.p = tau2.p,
          s2.lam = "fixed", d.lam = "fixed", nug.lam = "fixed", tau2.lam = "fixed")      
  gp_plot(gp, f, p, est$f, est$p, alt$f, alt$p, x_grid, obs, seed_i)
#  posteriors_plot(gp, priors) # needs trace=TRUE!
  OPT <- optimal_policy(gp, f, est$f, alt$f,
                        p, est$p, alt$p,
                        x_grid, h_grid, sigma_g, 
                        est$sigma_g, alt$sigma_g, 
                        delta, xT, profit, reward, OptTime)
  plot_policies(x_grid, OPT$gp_D, OPT$est_D, OPT$true_D, OPT$alt_D)
  dt <- simulate_opt(OPT, f, p, x_grid, h_grid, x0, z_g, profit)
  sim_plots(dt, seed=seed_i)
  profits_stats(dt)
})
```

![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-11-422b55d7bf-unnamed-chunk-11.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-15-422b55d7bf-unnamed-chunk-12.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-15-422b55d7bf-unnamed-chunk-13.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-16-422b55d7bf-unnamed-chunk-14.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-17-422b55d7bf-unnamed-chunk-15.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-18-422b55d7bf-unnamed-chunk-16.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-18-422b55d7bf-unnamed-chunk-17.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-19-422b55d7bf-unnamed-chunk-18.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-20-422b55d7bf-unnamed-chunk-19.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-20-422b55d7bf-unnamed-chunk-110.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-21-422b55d7bf-unnamed-chunk-111.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-22-422b55d7bf-unnamed-chunk-112.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-23-422b55d7bf-unnamed-chunk-113.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-23-422b55d7bf-unnamed-chunk-114.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-24-422b55d7bf-unnamed-chunk-115.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-25-422b55d7bf-unnamed-chunk-116.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-26-422b55d7bf-unnamed-chunk-117.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-26-422b55d7bf-unnamed-chunk-118.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-27-422b55d7bf-unnamed-chunk-119.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-28-422b55d7bf-unnamed-chunk-120.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-28-422b55d7bf-unnamed-chunk-121.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-29-422b55d7bf-unnamed-chunk-122.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-30-422b55d7bf-unnamed-chunk-123.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-31-422b55d7bf-unnamed-chunk-124.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-31-422b55d7bf-unnamed-chunk-125.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-32-422b55d7bf-unnamed-chunk-126.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-33-422b55d7bf-unnamed-chunk-127.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-34-422b55d7bf-unnamed-chunk-128.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-34-422b55d7bf-unnamed-chunk-129.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-35-422b55d7bf-unnamed-chunk-130.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-36-422b55d7bf-unnamed-chunk-131.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-37-422b55d7bf-unnamed-chunk-132.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-37-422b55d7bf-unnamed-chunk-133.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-38-422b55d7bf-unnamed-chunk-134.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-39-422b55d7bf-unnamed-chunk-135.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-39-422b55d7bf-unnamed-chunk-136.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-40-422b55d7bf-unnamed-chunk-137.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-41-422b55d7bf-unnamed-chunk-138.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-41-422b55d7bf-unnamed-chunk-139.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-42-422b55d7bf-unnamed-chunk-140.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-43-422b55d7bf-unnamed-chunk-141.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-44-422b55d7bf-unnamed-chunk-142.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-44-422b55d7bf-unnamed-chunk-143.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-45-422b55d7bf-unnamed-chunk-144.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-46-422b55d7bf-unnamed-chunk-145.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-46-422b55d7bf-unnamed-chunk-146.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-47-422b55d7bf-unnamed-chunk-147.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-48-422b55d7bf-unnamed-chunk-148.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-49-422b55d7bf-unnamed-chunk-149.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-49-422b55d7bf-unnamed-chunk-150.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-50-422b55d7bf-unnamed-chunk-151.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-51-422b55d7bf-unnamed-chunk-152.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-51-422b55d7bf-unnamed-chunk-153.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-52-422b55d7bf-unnamed-chunk-154.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-53-422b55d7bf-unnamed-chunk-155.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-54-422b55d7bf-unnamed-chunk-156.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-54-422b55d7bf-unnamed-chunk-157.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-55-422b55d7bf-unnamed-chunk-158.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-56-422b55d7bf-unnamed-chunk-159.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-57-422b55d7bf-unnamed-chunk-160.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-57-422b55d7bf-unnamed-chunk-161.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-58-422b55d7bf-unnamed-chunk-162.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-59-422b55d7bf-unnamed-chunk-163.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-11-59-422b55d7bf-unnamed-chunk-164.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-00-422b55d7bf-unnamed-chunk-165.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-01-422b55d7bf-unnamed-chunk-166.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-01-422b55d7bf-unnamed-chunk-167.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-02-422b55d7bf-unnamed-chunk-168.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-03-422b55d7bf-unnamed-chunk-169.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-04-422b55d7bf-unnamed-chunk-170.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-04-422b55d7bf-unnamed-chunk-171.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-05-422b55d7bf-unnamed-chunk-172.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-06-422b55d7bf-unnamed-chunk-173.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-06-422b55d7bf-unnamed-chunk-174.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-07-422b55d7bf-unnamed-chunk-175.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-08-422b55d7bf-unnamed-chunk-176.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-08-422b55d7bf-unnamed-chunk-177.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-09-422b55d7bf-unnamed-chunk-178.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-10-422b55d7bf-unnamed-chunk-179.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-11-422b55d7bf-unnamed-chunk-180.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-11-422b55d7bf-unnamed-chunk-181.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-12-422b55d7bf-unnamed-chunk-182.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-13-422b55d7bf-unnamed-chunk-183.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-14-422b55d7bf-unnamed-chunk-184.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-14-422b55d7bf-unnamed-chunk-185.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-15-422b55d7bf-unnamed-chunk-186.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-16-422b55d7bf-unnamed-chunk-187.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-16-422b55d7bf-unnamed-chunk-188.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-17-422b55d7bf-unnamed-chunk-189.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-18-422b55d7bf-unnamed-chunk-190.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-18-422b55d7bf-unnamed-chunk-191.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-19-422b55d7bf-unnamed-chunk-192.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-20-422b55d7bf-unnamed-chunk-193.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-21-422b55d7bf-unnamed-chunk-194.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-21-422b55d7bf-unnamed-chunk-195.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-22-422b55d7bf-unnamed-chunk-196.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-23-422b55d7bf-unnamed-chunk-197.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-24-422b55d7bf-unnamed-chunk-198.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-24-422b55d7bf-unnamed-chunk-199.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-25-422b55d7bf-unnamed-chunk-1100.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-26-422b55d7bf-unnamed-chunk-1101.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-27-422b55d7bf-unnamed-chunk-1102.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-27-422b55d7bf-unnamed-chunk-1103.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-28-422b55d7bf-unnamed-chunk-1104.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-29-422b55d7bf-unnamed-chunk-1105.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-29-422b55d7bf-unnamed-chunk-1106.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-30-422b55d7bf-unnamed-chunk-1107.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-31-422b55d7bf-unnamed-chunk-1108.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-32-422b55d7bf-unnamed-chunk-1109.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-32-422b55d7bf-unnamed-chunk-1110.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-33-422b55d7bf-unnamed-chunk-1111.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-34-422b55d7bf-unnamed-chunk-1112.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-34-422b55d7bf-unnamed-chunk-1113.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-35-422b55d7bf-unnamed-chunk-1114.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-36-422b55d7bf-unnamed-chunk-1115.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-37-422b55d7bf-unnamed-chunk-1116.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-37-422b55d7bf-unnamed-chunk-1117.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-38-422b55d7bf-unnamed-chunk-1118.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-39-422b55d7bf-unnamed-chunk-1119.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-39-422b55d7bf-unnamed-chunk-1120.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-40-422b55d7bf-unnamed-chunk-1121.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-41-422b55d7bf-unnamed-chunk-1122.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-42-422b55d7bf-unnamed-chunk-1123.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-43-422b55d7bf-unnamed-chunk-1124.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-43-422b55d7bf-unnamed-chunk-1125.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-44-422b55d7bf-unnamed-chunk-1126.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-45-422b55d7bf-unnamed-chunk-1127.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-13-12-46-422b55d7bf-unnamed-chunk-1128.png) 



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

```
stat_bin: binwidth defaulted to range/30. Use 'binwidth = x' to adjust
this.
```

![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-27-13-12-47-422b55d7bf-unnamed-chunk-2.png) 

```r
yields
```

```
        method      V1      sd L1
1           GP  8.1743 2.75359  1
2   Parametric 11.2512 3.52318  1
3         True 12.3135 3.53284  1
4   Structural  6.7277 0.11598  1
5           GP 12.4673 3.38532  2
6   Parametric 10.5957 2.96255  2
7         True 12.3135 3.53284  2
8   Structural  6.6315 0.12844  2
9           GP 10.5975 3.15331  3
10  Parametric 11.7120 3.39898  3
11        True 12.3135 3.53284  3
12  Structural  1.0350 0.70243  3
13          GP 11.2012 3.37961  4
14  Parametric 12.5258 3.42431  4
15        True 12.3135 3.53284  4
16  Structural  6.6045 0.48371  4
17          GP  9.8916 0.92801  5
18  Parametric  0.7275 0.41636  5
19        True 12.3135 3.53284  5
20  Structural  6.0168 0.41472  5
21          GP 12.1230 3.52055  6
22  Parametric 10.2726 3.30608  6
23        True 12.3135 3.53284  6
24  Structural  6.5506 0.71009  6
25          GP  8.6290 2.24891  7
26  Parametric 10.6570 3.42975  7
27        True 12.3135 3.53284  7
28  Structural 12.3165 3.62911  7
29          GP  6.7170 1.41535  8
30  Parametric  4.7092 2.51650  8
31        True 12.3135 3.53284  8
32  Structural  7.3424 2.13432  8
33          GP 11.6815 3.49706  9
34  Parametric 12.5720 3.43977  9
35        True 12.3135 3.53284  9
36  Structural  0.7275 0.41636  9
37          GP  6.7152 0.47269 10
38  Parametric  7.0349 1.62783 10
39        True 12.3135 3.53284 10
40  Structural  6.4578 0.95719 10
41          GP 12.3330 3.53975 11
42  Parametric 11.5639 3.44359 11
43        True 12.3135 3.53284 11
44  Structural  6.6037 0.86138 11
45          GP  9.4197 2.65413 12
46  Parametric 11.0692 3.72106 12
47        True 12.3135 3.53284 12
48  Structural  7.4790 2.99940 12
49          GP  7.1216 0.54003 13
50  Parametric 12.5585 3.68182 13
51        True 12.3135 3.53284 13
52  Structural  7.9442 2.37701 13
53          GP 11.9925 3.51266 14
54  Parametric 12.2634 3.44924 14
55        True 12.3135 3.53284 14
56  Structural  6.0615 1.45197 14
57          GP 11.8455 3.44512 15
58  Parametric 12.7011 3.58071 15
59        True 12.3135 3.53284 15
60  Structural 10.9445 3.74118 15
61          GP 11.3143 3.20557 16
62  Parametric 12.7008 3.64679 16
63        True 12.3135 3.53284 16
64  Structural  0.7275 0.41636 16
65          GP  6.7353 0.74380 17
66  Parametric 12.5585 3.68182 17
67        True 12.3135 3.53284 17
68  Structural  0.7275 0.41636 17
69          GP  8.4668 2.52583 18
70  Parametric 13.0837 3.54547 18
71        True 12.3135 3.53284 18
72  Structural  2.8530 1.69652 18
73          GP 11.5949 3.40676 19
74  Parametric 12.2068 3.44594 19
75        True 12.3135 3.53284 19
76  Structural  6.8082 0.28163 19
77          GP  7.4435 0.93481 20
78  Parametric 10.1315 3.36812 20
79        True 12.3135 3.53284 20
80  Structural  6.4497 0.84411 20
81          GP  7.8248 2.16258 21
82  Parametric 12.2430 3.58411 21
83        True 12.3135 3.53284 21
84  Structural  6.7633 0.77099 21
85          GP  7.2139 1.10800 22
86  Parametric 12.5585 3.68182 22
87        True 12.3135 3.53284 22
88  Structural  6.6075 0.07500 22
89          GP  6.9532 0.81113 23
90  Parametric 13.0572 3.72219 23
91        True 12.3135 3.53284 23
92  Structural  8.3906 2.74473 23
93          GP  9.2019 3.08996 24
94  Parametric 11.8464 3.67216 24
95        True 12.3135 3.53284 24
96  Structural  6.4065 0.36691 24
97          GP  8.8296 1.93650 25
98  Parametric 11.8236 3.45687 25
99        True 12.3135 3.53284 25
100 Structural  6.3948 0.40833 25
101         GP  6.1680 0.43331 26
102 Parametric 10.4232 3.46401 26
103       True 12.3135 3.53284 26
104 Structural  0.7395 0.42919 26
105         GP 12.0930 3.43933 27
106 Parametric 12.4813 3.51280 27
107       True 12.3135 3.53284 27
108 Structural  9.1838 3.02881 27
109         GP  7.7912 0.79381 28
110 Parametric  0.7275 0.41636 28
111       True 12.3135 3.53284 28
112 Structural  1.0350 0.70243 28
113         GP 11.9623 3.33112 29
114 Parametric 13.1633 3.56114 29
115       True 12.3135 3.53284 29
116 Structural  0.7275 0.41636 29
117         GP 10.8869 3.37941 30
118 Parametric  5.6995 2.76630 30
119       True 12.3135 3.53284 30
120 Structural  7.0854 1.48048 30
121         GP  6.7492 1.82674 31
122 Parametric 10.5915 3.45608 31
123       True 12.3135 3.53284 31
124 Structural  6.7953 0.39640 31
125         GP 10.2450 3.36631 32
126 Parametric 12.3673 3.43741 32
127       True 12.3135 3.53284 32
128 Structural  6.8852 0.08156 32
```




NULL

