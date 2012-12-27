










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
OptTime = 50  # stationarity with unstable models is tricky thing
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
sigma_m <- 0.0
x_grid <- seq(0, 1.5 * K, length=101)
h_grid <- x_grid
```


With parameters 1.1, 10, 5. 



```r
seed <- round(runif(32) * 1e6)
seed[1] <- 1
seed
```

```
 [1]      1 316507 231725 255413 335007 381962 238977 352168 982046 303711
[11] 278403 683508 444264 549747 605184 361138 976276 745228 731537 733902
[21] 189338  25004 925384 389512 259166 301249 251023 696861 440810 347912
[31] 563003 458116
```

```r
yields <- lapply(seed, 
function(seed_i){
  Xo <- K # observations start from
  x0 <- Xo # simulation under policy starts from
  obs <- sim_obs(Xo, z_g, f, p, Tobs=35, nz= 1, 
                 harvest = sort(rep(seq(0, .5, length=7), 5)), seed = seed_i)
    
  alt <- par_est(obs,  init = c(r=1.5, K=mean(obs$x), s=1))
  est <- par_est_allee(obs, f, p,  init = c(p[1]+1, p[2]-1, p[3]-2, 
                                   s = sigma_g + abs(rnorm(1,0,.1))))
  
  
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

![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-27-e114ca6d34-unnamed-chunk-11.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-31-e114ca6d34-unnamed-chunk-12.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-32-e114ca6d34-unnamed-chunk-13.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-33-e114ca6d34-unnamed-chunk-14.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-34-e114ca6d34-unnamed-chunk-15.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-34-e114ca6d34-unnamed-chunk-16.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-35-e114ca6d34-unnamed-chunk-17.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-36-e114ca6d34-unnamed-chunk-18.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-37-e114ca6d34-unnamed-chunk-19.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-38-e114ca6d34-unnamed-chunk-110.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-38-e114ca6d34-unnamed-chunk-111.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-39-e114ca6d34-unnamed-chunk-112.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-40-e114ca6d34-unnamed-chunk-113.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-41-e114ca6d34-unnamed-chunk-114.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-41-e114ca6d34-unnamed-chunk-115.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-42-e114ca6d34-unnamed-chunk-116.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-43-e114ca6d34-unnamed-chunk-117.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-44-e114ca6d34-unnamed-chunk-118.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-45-e114ca6d34-unnamed-chunk-119.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-46-e114ca6d34-unnamed-chunk-120.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-46-e114ca6d34-unnamed-chunk-121.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-47-e114ca6d34-unnamed-chunk-122.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-48-e114ca6d34-unnamed-chunk-123.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-49-e114ca6d34-unnamed-chunk-124.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-49-e114ca6d34-unnamed-chunk-125.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-50-e114ca6d34-unnamed-chunk-126.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-51-e114ca6d34-unnamed-chunk-127.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-52-e114ca6d34-unnamed-chunk-128.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-52-e114ca6d34-unnamed-chunk-129.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-53-e114ca6d34-unnamed-chunk-130.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-54-e114ca6d34-unnamed-chunk-131.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-55-e114ca6d34-unnamed-chunk-132.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-56-e114ca6d34-unnamed-chunk-133.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-56-e114ca6d34-unnamed-chunk-134.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-57-e114ca6d34-unnamed-chunk-135.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-58-e114ca6d34-unnamed-chunk-136.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-55-59-e114ca6d34-unnamed-chunk-137.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-00-e114ca6d34-unnamed-chunk-138.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-00-e114ca6d34-unnamed-chunk-139.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-01-e114ca6d34-unnamed-chunk-140.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-02-e114ca6d34-unnamed-chunk-141.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-03-e114ca6d34-unnamed-chunk-142.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-03-e114ca6d34-unnamed-chunk-143.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-04-e114ca6d34-unnamed-chunk-144.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-05-e114ca6d34-unnamed-chunk-145.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-06-e114ca6d34-unnamed-chunk-146.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-07-e114ca6d34-unnamed-chunk-147.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-07-e114ca6d34-unnamed-chunk-148.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-08-e114ca6d34-unnamed-chunk-149.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-09-e114ca6d34-unnamed-chunk-150.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-10-e114ca6d34-unnamed-chunk-151.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-11-e114ca6d34-unnamed-chunk-152.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-11-e114ca6d34-unnamed-chunk-153.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-12-e114ca6d34-unnamed-chunk-154.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-13-e114ca6d34-unnamed-chunk-155.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-14-e114ca6d34-unnamed-chunk-156.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-14-e114ca6d34-unnamed-chunk-157.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-15-e114ca6d34-unnamed-chunk-158.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-16-e114ca6d34-unnamed-chunk-159.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-17-e114ca6d34-unnamed-chunk-160.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-18-e114ca6d34-unnamed-chunk-161.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-18-e114ca6d34-unnamed-chunk-162.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-19-e114ca6d34-unnamed-chunk-163.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-20-e114ca6d34-unnamed-chunk-164.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-21-e114ca6d34-unnamed-chunk-165.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-21-e114ca6d34-unnamed-chunk-166.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-22-e114ca6d34-unnamed-chunk-167.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-23-e114ca6d34-unnamed-chunk-168.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-24-e114ca6d34-unnamed-chunk-169.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-25-e114ca6d34-unnamed-chunk-170.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-25-e114ca6d34-unnamed-chunk-171.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-26-e114ca6d34-unnamed-chunk-172.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-27-e114ca6d34-unnamed-chunk-173.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-28-e114ca6d34-unnamed-chunk-174.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-28-e114ca6d34-unnamed-chunk-175.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-29-e114ca6d34-unnamed-chunk-176.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-30-e114ca6d34-unnamed-chunk-177.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-31-e114ca6d34-unnamed-chunk-178.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-31-e114ca6d34-unnamed-chunk-179.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-32-e114ca6d34-unnamed-chunk-180.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-33-e114ca6d34-unnamed-chunk-181.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-34-e114ca6d34-unnamed-chunk-182.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-35-e114ca6d34-unnamed-chunk-183.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-36-e114ca6d34-unnamed-chunk-184.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-36-e114ca6d34-unnamed-chunk-185.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-37-e114ca6d34-unnamed-chunk-186.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-38-e114ca6d34-unnamed-chunk-187.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-39-e114ca6d34-unnamed-chunk-188.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-39-e114ca6d34-unnamed-chunk-189.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-40-e114ca6d34-unnamed-chunk-190.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-41-e114ca6d34-unnamed-chunk-191.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-42-e114ca6d34-unnamed-chunk-192.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-43-e114ca6d34-unnamed-chunk-193.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-44-e114ca6d34-unnamed-chunk-194.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-44-e114ca6d34-unnamed-chunk-195.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-45-e114ca6d34-unnamed-chunk-196.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-46-e114ca6d34-unnamed-chunk-197.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-47-e114ca6d34-unnamed-chunk-198.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-48-e114ca6d34-unnamed-chunk-199.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-48-e114ca6d34-unnamed-chunk-1100.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-49-e114ca6d34-unnamed-chunk-1101.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-50-e114ca6d34-unnamed-chunk-1102.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-51-e114ca6d34-unnamed-chunk-1103.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-52-e114ca6d34-unnamed-chunk-1104.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-52-e114ca6d34-unnamed-chunk-1105.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-53-e114ca6d34-unnamed-chunk-1106.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-54-e114ca6d34-unnamed-chunk-1107.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-55-e114ca6d34-unnamed-chunk-1108.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-56-e114ca6d34-unnamed-chunk-1109.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-56-e114ca6d34-unnamed-chunk-1110.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-57-e114ca6d34-unnamed-chunk-1111.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-58-e114ca6d34-unnamed-chunk-1112.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-56-59-e114ca6d34-unnamed-chunk-1113.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-57-00-e114ca6d34-unnamed-chunk-1114.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-57-00-e114ca6d34-unnamed-chunk-1115.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-57-01-e114ca6d34-unnamed-chunk-1116.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-57-02-e114ca6d34-unnamed-chunk-1117.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-57-03-e114ca6d34-unnamed-chunk-1118.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-57-03-e114ca6d34-unnamed-chunk-1119.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-57-04-e114ca6d34-unnamed-chunk-1120.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-57-05-e114ca6d34-unnamed-chunk-1121.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-57-06-e114ca6d34-unnamed-chunk-1122.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-57-06-e114ca6d34-unnamed-chunk-1123.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-57-07-e114ca6d34-unnamed-chunk-1124.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-57-08-e114ca6d34-unnamed-chunk-1125.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-57-09-e114ca6d34-unnamed-chunk-1126.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-57-10-e114ca6d34-unnamed-chunk-1127.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-14-57-11-e114ca6d34-unnamed-chunk-1128.png) 



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

![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-27-14-57-12-e114ca6d34-unnamed-chunk-2.png) 

```r
yields
```

```
        method     V1       sd L1
1           GP 27.102 3.216427  1
2   Parametric 23.235 2.949384  1
3         True 28.773 3.099881  1
4   Structural  6.598 0.003606  1
5           GP 28.608 3.102644  2
6   Parametric 21.296 2.631665  2
7         True 28.773 3.099881  2
8   Structural  6.005 0.025717  2
9           GP 27.841 3.175109  3
10  Parametric 15.474 3.058871  3
11        True 28.773 3.099881  3
12  Structural  6.728 0.029762  3
13          GP 25.944 2.886719  4
14  Parametric 27.197 2.952171  4
15        True 28.773 3.099881  4
16  Structural  6.413 0.124446  4
17          GP 27.442 3.513024  5
18  Parametric 27.312 2.863777  5
19        True 28.773 3.099881  5
20  Structural  6.057 0.145821  5
21          GP 28.431 3.141107  6
22  Parametric 28.092 3.020316  6
23        True 28.773 3.099881  6
24  Structural  6.312 0.109211  6
25          GP 27.102 3.216427  7
26  Parametric 25.647 2.897803  7
27        True 28.773 3.099881  7
28  Structural  6.139 0.059479  7
29          GP 24.820 3.266041  8
30  Parametric 28.461 3.007462  8
31        True 28.773 3.099881  8
32  Structural  6.821 0.069591  8
33          GP 28.782 2.998581  9
34  Parametric 23.343 3.256586  9
35        True 28.773 3.099881  9
36  Structural  6.364 0.670677  9
37          GP 23.634 7.467731 10
38  Parametric 28.462 3.126818 10
39        True 28.773 3.099881 10
40  Structural  0.000 0.000000 10
41          GP 25.360 6.561896 11
42  Parametric 28.151 2.940474 11
43        True 28.773 3.099881 11
44  Structural  6.158 0.125858 11
45          GP 27.184 3.200700 12
46  Parametric 28.083 3.025391 12
47        True 28.773 3.099881 12
48  Structural  6.642 0.068372 12
49          GP 27.930 3.169600 13
50  Parametric  0.000 0.000000 13
51        True 28.773 3.099881 13
52  Structural  6.495 0.130382 13
53          GP 28.431 3.141107 14
54  Parametric 26.490 2.794638 14
55        True 28.773 3.099881 14
56  Structural  6.652 0.136188 14
57          GP 28.393 3.498655 15
58  Parametric 28.488 2.953786 15
59        True 28.773 3.099881 15
60  Structural  6.489 0.098961 15
61          GP 28.431 3.141107 16
62  Parametric 28.433 3.006310 16
63        True 28.773 3.099881 16
64  Structural  6.371 0.954187 16
65          GP 28.090 3.131719 17
66  Parametric  0.000 0.000000 17
67        True 28.773 3.099881 17
68  Structural  6.582 0.028389 17
69          GP 28.644 3.056879 18
70  Parametric 25.647 2.897803 18
71        True 28.773 3.099881 18
72  Structural  6.585 0.117327 18
73          GP 27.183 3.200467 19
74  Parametric 26.490 2.794638 19
75        True 28.773 3.099881 19
76  Structural  6.384 0.100725 19
77          GP 27.885 3.188937 20
78  Parametric 27.195 2.844652 20
79        True 28.773 3.099881 20
80  Structural  6.501 0.108095 20
81          GP 27.616 3.153706 21
82  Parametric 28.486 2.955895 21
83        True 28.773 3.099881 21
84  Structural  6.489 0.098961 21
85          GP 27.886 3.188765 22
86  Parametric 27.197 2.952171 22
87        True 28.773 3.099881 22
88  Structural  6.393 0.141005 22
89          GP 26.079 3.247658 23
90  Parametric 27.622 2.897677 23
91        True 28.773 3.099881 23
92  Structural  6.273 0.150286 23
93          GP 26.079 3.247658 24
94  Parametric 28.607 3.068853 24
95        True 28.773 3.099881 24
96  Structural  0.000 0.000000 24
97          GP 23.185 7.772252 25
98  Parametric 23.988 3.252459 25
99        True 28.773 3.099881 25
100 Structural  7.950 0.000000 25
101         GP 24.843 3.284314 26
102 Parametric 27.165 3.166766 26
103       True 28.773 3.099881 26
104 Structural  6.823 0.103859 26
105         GP 26.710 3.154380 27
106 Parametric 24.866 2.927600 27
107       True 28.773 3.099881 27
108 Structural 28.548 2.969457 27
109         GP 27.885 3.188937 28
110 Parametric 28.083 3.025391 28
111       True 28.773 3.099881 28
112 Structural 26.490 2.794638 28
113         GP 27.387 3.164446 29
114 Parametric 27.427 3.212927 29
115       True 28.773 3.099881 29
116 Structural  6.213 0.166761 29
117         GP 21.231 7.913535 30
118 Parametric 25.485 2.738156 30
119       True 28.773 3.099881 30
120 Structural  0.000 0.000000 30
121         GP 24.820 3.266041 31
122 Parametric 28.192 2.945673 31
123       True 28.773 3.099881 31
124 Structural  6.547 0.121039 31
125         GP 27.102 3.216427 32
126 Parametric 27.418 3.174135 32
127       True 28.773 3.099881 32
128 Structural  6.507 0.073175 32
```




NULL

