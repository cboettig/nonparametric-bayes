










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
sigma_g <- 0.05
sigma_m <- 0.02
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
 [1]      1 980427   1000 990863 241078 901116 656932 814497 720135 318268
[11]  48137 574809 722140   9799 153557 201288 215812 324127 304851 791412
[21] 336145 183942 691199 613993 784611 136776 322899 434786 998680 750979
[31] 950859 816359
```

```r
yields <- lapply(seed, 
function(seed_i){
  Xo <- K # observations start from
  x0 <- Xo # simulation under policy starts from
  obs <- sim_obs(Xo, z_g, f, p, Tobs=40, seed = seed_i)
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

![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-21-40-4836f1e85b-unnamed-chunk-11.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-21-40-4836f1e85b-unnamed-chunk-12.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-21-41-4836f1e85b-unnamed-chunk-13.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-21-42-4836f1e85b-unnamed-chunk-14.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-21-42-4836f1e85b-unnamed-chunk-15.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-21-42-4836f1e85b-unnamed-chunk-16.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-21-43-4836f1e85b-unnamed-chunk-17.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-21-44-4836f1e85b-unnamed-chunk-18.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-21-44-4836f1e85b-unnamed-chunk-19.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-21-45-4836f1e85b-unnamed-chunk-110.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-21-45-4836f1e85b-unnamed-chunk-111.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-21-46-4836f1e85b-unnamed-chunk-112.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-21-46-4836f1e85b-unnamed-chunk-113.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-21-47-4836f1e85b-unnamed-chunk-114.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-21-47-4836f1e85b-unnamed-chunk-115.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-21-48-4836f1e85b-unnamed-chunk-116.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-21-49-4836f1e85b-unnamed-chunk-117.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-21-54-4836f1e85b-unnamed-chunk-118.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-21-55-4836f1e85b-unnamed-chunk-119.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-21-55-4836f1e85b-unnamed-chunk-120.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-21-56-4836f1e85b-unnamed-chunk-121.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-21-56-4836f1e85b-unnamed-chunk-122.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-21-57-4836f1e85b-unnamed-chunk-123.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-21-57-4836f1e85b-unnamed-chunk-124.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-21-58-4836f1e85b-unnamed-chunk-125.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-21-58-4836f1e85b-unnamed-chunk-126.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-21-59-4836f1e85b-unnamed-chunk-127.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-21-59-4836f1e85b-unnamed-chunk-128.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-00-4836f1e85b-unnamed-chunk-129.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-00-4836f1e85b-unnamed-chunk-130.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-01-4836f1e85b-unnamed-chunk-131.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-01-4836f1e85b-unnamed-chunk-132.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-02-4836f1e85b-unnamed-chunk-133.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-02-4836f1e85b-unnamed-chunk-134.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-03-4836f1e85b-unnamed-chunk-135.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-03-4836f1e85b-unnamed-chunk-136.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-04-4836f1e85b-unnamed-chunk-137.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-04-4836f1e85b-unnamed-chunk-138.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-05-4836f1e85b-unnamed-chunk-139.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-05-4836f1e85b-unnamed-chunk-140.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-06-4836f1e85b-unnamed-chunk-141.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-06-4836f1e85b-unnamed-chunk-142.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-07-4836f1e85b-unnamed-chunk-143.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-07-4836f1e85b-unnamed-chunk-144.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-08-4836f1e85b-unnamed-chunk-145.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-08-4836f1e85b-unnamed-chunk-146.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-09-4836f1e85b-unnamed-chunk-147.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-10-4836f1e85b-unnamed-chunk-148.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-10-4836f1e85b-unnamed-chunk-149.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-10-4836f1e85b-unnamed-chunk-150.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-11-4836f1e85b-unnamed-chunk-151.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-12-4836f1e85b-unnamed-chunk-152.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-12-4836f1e85b-unnamed-chunk-153.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-12-4836f1e85b-unnamed-chunk-154.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-13-4836f1e85b-unnamed-chunk-155.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-14-4836f1e85b-unnamed-chunk-156.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-14-4836f1e85b-unnamed-chunk-157.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-15-4836f1e85b-unnamed-chunk-158.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-15-4836f1e85b-unnamed-chunk-159.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-16-4836f1e85b-unnamed-chunk-160.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-16-4836f1e85b-unnamed-chunk-161.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-17-4836f1e85b-unnamed-chunk-162.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-17-4836f1e85b-unnamed-chunk-163.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-18-4836f1e85b-unnamed-chunk-164.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-18-4836f1e85b-unnamed-chunk-165.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-24-4836f1e85b-unnamed-chunk-166.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-24-4836f1e85b-unnamed-chunk-167.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-25-4836f1e85b-unnamed-chunk-168.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-25-4836f1e85b-unnamed-chunk-169.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-26-4836f1e85b-unnamed-chunk-170.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-26-4836f1e85b-unnamed-chunk-171.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-27-4836f1e85b-unnamed-chunk-172.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-27-4836f1e85b-unnamed-chunk-173.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-28-4836f1e85b-unnamed-chunk-174.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-28-4836f1e85b-unnamed-chunk-175.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-29-4836f1e85b-unnamed-chunk-176.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-29-4836f1e85b-unnamed-chunk-177.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-30-4836f1e85b-unnamed-chunk-178.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-30-4836f1e85b-unnamed-chunk-179.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-31-4836f1e85b-unnamed-chunk-180.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-31-4836f1e85b-unnamed-chunk-181.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-32-4836f1e85b-unnamed-chunk-182.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-32-4836f1e85b-unnamed-chunk-183.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-33-4836f1e85b-unnamed-chunk-184.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-33-4836f1e85b-unnamed-chunk-185.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-34-4836f1e85b-unnamed-chunk-186.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-39-4836f1e85b-unnamed-chunk-187.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-40-4836f1e85b-unnamed-chunk-188.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-40-4836f1e85b-unnamed-chunk-189.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-41-4836f1e85b-unnamed-chunk-190.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-41-4836f1e85b-unnamed-chunk-191.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-42-4836f1e85b-unnamed-chunk-192.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-42-4836f1e85b-unnamed-chunk-193.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-43-4836f1e85b-unnamed-chunk-194.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-44-4836f1e85b-unnamed-chunk-195.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-44-4836f1e85b-unnamed-chunk-196.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-45-4836f1e85b-unnamed-chunk-197.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-45-4836f1e85b-unnamed-chunk-198.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-46-4836f1e85b-unnamed-chunk-199.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-46-4836f1e85b-unnamed-chunk-1100.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-47-4836f1e85b-unnamed-chunk-1101.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-47-4836f1e85b-unnamed-chunk-1102.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-48-4836f1e85b-unnamed-chunk-1103.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-48-4836f1e85b-unnamed-chunk-1104.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-49-4836f1e85b-unnamed-chunk-1105.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-49-4836f1e85b-unnamed-chunk-1106.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-55-4836f1e85b-unnamed-chunk-1107.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-55-4836f1e85b-unnamed-chunk-1108.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-56-4836f1e85b-unnamed-chunk-1109.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-56-4836f1e85b-unnamed-chunk-1110.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-57-4836f1e85b-unnamed-chunk-1111.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-57-4836f1e85b-unnamed-chunk-1112.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-58-4836f1e85b-unnamed-chunk-1113.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-58-4836f1e85b-unnamed-chunk-1114.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-59-4836f1e85b-unnamed-chunk-1115.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-22-59-4836f1e85b-unnamed-chunk-1116.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-23-00-4836f1e85b-unnamed-chunk-1117.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-23-00-4836f1e85b-unnamed-chunk-1118.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-23-01-4836f1e85b-unnamed-chunk-1119.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-23-02-4836f1e85b-unnamed-chunk-1120.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-23-02-4836f1e85b-unnamed-chunk-1121.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-23-02-4836f1e85b-unnamed-chunk-1122.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-23-03-4836f1e85b-unnamed-chunk-1123.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-23-04-4836f1e85b-unnamed-chunk-1124.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-23-04-4836f1e85b-unnamed-chunk-1125.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-23-04-4836f1e85b-unnamed-chunk-1126.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-23-05-4836f1e85b-unnamed-chunk-1127.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-24-13-23-06-4836f1e85b-unnamed-chunk-1128.png) 



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

![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-24-13-23-07-4836f1e85b-unnamed-chunk-2.png) 

```r
yields
```

```
        method      V1      sd L1
1           GP  9.5130 1.45073  1
2   Parametric  0.4995 0.43540  1
3         True 12.0060 1.72349  1
4   Structural  9.7860 1.40779  1
5           GP 10.9770 1.48666  2
6   Parametric  0.2850 0.39629  2
7         True 12.0060 1.72349  2
8   Structural  6.1947 0.75683  2
9           GP  8.4538 1.49893  3
10  Parametric  5.1855 1.60429  3
11        True 12.0060 1.72349  3
12  Structural  7.1189 1.13594  3
13          GP  7.3660 1.47742  4
14  Parametric  0.0000 0.00000  4
15        True 12.0060 1.72349  4
16  Structural  8.2845 1.50406  4
17          GP  6.7679 0.96193  5
18  Parametric  0.7785 0.68014  5
19        True 12.0060 1.72349  5
20  Structural  8.7105 1.79814  5
21          GP  8.6535 1.51713  6
22  Parametric  0.0810 0.20385  6
23        True 12.0060 1.72349  6
24  Structural  6.7904 0.10826  6
25          GP 11.9400 1.58329  7
26  Parametric  5.5703 4.29215  7
27        True 12.0060 1.72349  7
28  Structural  6.1560 0.04729  7
29          GP  8.1030 1.53341  8
30  Parametric  0.2850 0.39629  8
31        True 12.0060 1.72349  8
32  Structural 11.5035 1.75961  8
33          GP  7.0002 1.15280  9
34  Parametric  9.4395 1.81962  9
35        True 12.0060 1.72349  9
36  Structural  5.3535 0.39063  9
37          GP 11.5650 1.54295 10
38  Parametric  7.6020 1.68966 10
39        True 12.0060 1.72349 10
40  Structural  9.5475 1.81313 10
41          GP  6.5178 0.76158 11
42  Parametric  7.9110 1.78193 11
43        True 12.0060 1.72349 11
44  Structural  6.9417 0.09921 11
45          GP  9.0630 1.54907 12
46  Parametric  1.7850 1.01523 12
47        True 12.0060 1.72349 12
48  Structural  8.7105 1.79814 12
49          GP 10.6320 1.46943 13
50  Parametric  6.8265 1.60206 13
51        True 12.0060 1.72349 13
52  Structural  0.9465 0.64410 13
53          GP 11.0055 1.50811 14
54  Parametric 10.0815 1.79953 14
55        True 12.0060 1.72349 14
56  Structural  0.0000 0.00000 14
57          GP  8.6580 1.46289 15
58  Parametric  7.0455 1.74947 15
59        True 12.0060 1.72349 15
60  Structural  9.5475 1.81313 15
61          GP 11.6400 1.54361 16
62  Parametric  0.0000 0.00000 16
63        True 12.0060 1.72349 16
64  Structural  0.0000 0.00000 16
65          GP  8.1135 1.49151 17
66  Parametric  0.3405 0.29454 17
67        True 12.0060 1.72349 17
68  Structural  0.3405 0.29454 17
69          GP  9.2925 1.44086 18
70  Parametric  9.0075 1.80270 18
71        True 12.0060 1.72349 18
72  Structural  6.5148 0.18002 18
73          GP  8.5740 1.50736 19
74  Parametric  7.9110 1.78193 19
75        True 12.0060 1.72349 19
76  Structural  8.7855 1.67986 19
77          GP  9.8310 1.43201 20
78  Parametric  9.5475 1.81313 20
79        True 12.0060 1.72349 20
80  Structural  6.5148 0.18002 20
81          GP  6.2029 0.41585 21
82  Parametric 10.0000 0.00000 21
83        True 12.0060 1.72349 21
84  Structural  0.0000 0.00000 21
85          GP 11.6040 1.54059 22
86  Parametric  0.0000 0.00000 22
87        True 12.0060 1.72349 22
88  Structural  0.3405 0.29454 22
89          GP 11.9940 1.58644 23
90  Parametric 10.6365 1.80543 23
91        True 12.0060 1.72349 23
92  Structural  9.1905 1.51760 23
93          GP  8.1090 1.54520 24
94  Parametric  5.1855 1.60429 24
95        True 12.0060 1.72349 24
96  Structural  6.5009 0.92136 24
97          GP  7.1721 1.23218 25
98  Parametric 10.6365 1.80543 25
99        True 12.0060 1.72349 25
100 Structural  9.1725 1.79334 25
101         GP 10.7295 1.46791 26
102 Parametric  2.6460 1.11814 26
103       True 12.0060 1.72349 26
104 Structural  6.5330 0.91451 26
105         GP  9.6450 1.45594 27
106 Parametric  0.0000 0.00000 27
107       True 12.0060 1.72349 27
108 Structural 11.1060 1.49528 27
109         GP  8.6715 1.50268 28
110 Parametric  4.2345 1.49497 28
111       True 12.0060 1.72349 28
112 Structural  6.8670 1.20233 28
113         GP  9.6930 1.48321 29
114 Parametric 10.0815 1.79953 29
115       True 12.0060 1.72349 29
116 Structural  6.6330 1.40851 29
117         GP  8.8080 1.50953 30
118 Parametric  5.1855 1.60429 30
119       True 12.0060 1.72349 30
120 Structural  1.7145 0.64359 30
121         GP  8.2648 1.49954 31
122 Parametric  0.3405 0.29454 31
123       True 12.0060 1.72349 31
124 Structural  6.5148 0.18002 31
125         GP 10.0320 1.44746 32
126 Parametric  8.7105 1.79814 32
127       True 12.0060 1.72349 32
128 Structural  7.2540 1.54221 32
```




NULL

