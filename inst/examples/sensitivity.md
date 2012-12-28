










Fixed priors on hyperparameters, fixed model type.


```r
#inv gamma has mean b / (a - 1) (assuming a>1) and variance b ^ 2 / ((a - 2) * (a - 1) ^ 2) (assuming a>2)
s2.p <- c(5,5)  
tau2.p <- c(5,1)
d.p = c(10, 1/0.1, 10, 1/0.1)
nug.p = c(10, 1/0.1, 10, 1/0.1) # gamma mean
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
```




```r
f <- RickerAllee
p <- c(1.5, 10, 5) 
K <- 10
allee <- 5
```




```r
sigma_g <- 0.1
sigma_m <- 0.0
z_g = function() rlnorm(1, 0, sigma_g)
z_m = function() 1+(2*runif(1, 0,  1)-1) * sigma_m
x_grid <- seq(0, 1.5 * K, length=101)
h_grid <- x_grid
```


With parameters 1.5, 10, 5. 



```r
seed <- round(runif(32) * 1e6)
seed[1] <- 1
seed
```

```
 [1]      1 136332 327822 431050  59473 695943 408795 279337 328374 301268
[11] 314435  47331 423943 880512 846295 934813  52123 953486 318780  89776
[21]  51710 282494 153598 553654 176856 839022  27643 600168 827952 731727
[31] 359673 525685
```

```r
yields <- lapply(seed, 
function(seed_i){
  Xo <- K # observations start from
  x0 <- Xo # simulation under policy starts from
  obs <- sim_obs(Xo, z_g, f, p, Tobs=35, nz= 10, 
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

![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-51-58-c4cdce9c3b-unnamed-chunk-11.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-51-59-c4cdce9c3b-unnamed-chunk-12.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-00-c4cdce9c3b-unnamed-chunk-13.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-01-c4cdce9c3b-unnamed-chunk-14.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-01-c4cdce9c3b-unnamed-chunk-15.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-02-c4cdce9c3b-unnamed-chunk-16.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-03-c4cdce9c3b-unnamed-chunk-17.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-03-c4cdce9c3b-unnamed-chunk-18.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-04-c4cdce9c3b-unnamed-chunk-19.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-05-c4cdce9c3b-unnamed-chunk-110.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-05-c4cdce9c3b-unnamed-chunk-111.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-06-c4cdce9c3b-unnamed-chunk-112.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-07-c4cdce9c3b-unnamed-chunk-113.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-08-c4cdce9c3b-unnamed-chunk-114.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-08-c4cdce9c3b-unnamed-chunk-115.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-09-c4cdce9c3b-unnamed-chunk-116.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-10-c4cdce9c3b-unnamed-chunk-117.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-10-c4cdce9c3b-unnamed-chunk-118.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-11-c4cdce9c3b-unnamed-chunk-119.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-12-c4cdce9c3b-unnamed-chunk-120.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-13-c4cdce9c3b-unnamed-chunk-121.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-13-c4cdce9c3b-unnamed-chunk-122.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-14-c4cdce9c3b-unnamed-chunk-123.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-15-c4cdce9c3b-unnamed-chunk-124.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-15-c4cdce9c3b-unnamed-chunk-125.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-16-c4cdce9c3b-unnamed-chunk-126.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-17-c4cdce9c3b-unnamed-chunk-127.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-18-c4cdce9c3b-unnamed-chunk-128.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-18-c4cdce9c3b-unnamed-chunk-129.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-19-c4cdce9c3b-unnamed-chunk-130.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-20-c4cdce9c3b-unnamed-chunk-131.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-21-c4cdce9c3b-unnamed-chunk-132.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-21-c4cdce9c3b-unnamed-chunk-133.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-22-c4cdce9c3b-unnamed-chunk-134.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-22-c4cdce9c3b-unnamed-chunk-135.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-23-c4cdce9c3b-unnamed-chunk-136.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-24-c4cdce9c3b-unnamed-chunk-137.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-25-c4cdce9c3b-unnamed-chunk-138.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-25-c4cdce9c3b-unnamed-chunk-139.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-26-c4cdce9c3b-unnamed-chunk-140.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-27-c4cdce9c3b-unnamed-chunk-141.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-28-c4cdce9c3b-unnamed-chunk-142.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-28-c4cdce9c3b-unnamed-chunk-143.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-29-c4cdce9c3b-unnamed-chunk-144.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-30-c4cdce9c3b-unnamed-chunk-145.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-30-c4cdce9c3b-unnamed-chunk-146.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-31-c4cdce9c3b-unnamed-chunk-147.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-32-c4cdce9c3b-unnamed-chunk-148.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-33-c4cdce9c3b-unnamed-chunk-149.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-33-c4cdce9c3b-unnamed-chunk-150.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-34-c4cdce9c3b-unnamed-chunk-151.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-35-c4cdce9c3b-unnamed-chunk-152.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-36-c4cdce9c3b-unnamed-chunk-153.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-37-c4cdce9c3b-unnamed-chunk-154.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-37-c4cdce9c3b-unnamed-chunk-155.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-38-c4cdce9c3b-unnamed-chunk-156.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-39-c4cdce9c3b-unnamed-chunk-157.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-39-c4cdce9c3b-unnamed-chunk-158.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-40-c4cdce9c3b-unnamed-chunk-159.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-41-c4cdce9c3b-unnamed-chunk-160.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-42-c4cdce9c3b-unnamed-chunk-161.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-42-c4cdce9c3b-unnamed-chunk-162.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-43-c4cdce9c3b-unnamed-chunk-163.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-44-c4cdce9c3b-unnamed-chunk-164.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-45-c4cdce9c3b-unnamed-chunk-165.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-45-c4cdce9c3b-unnamed-chunk-166.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-46-c4cdce9c3b-unnamed-chunk-167.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-47-c4cdce9c3b-unnamed-chunk-168.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-47-c4cdce9c3b-unnamed-chunk-169.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-48-c4cdce9c3b-unnamed-chunk-170.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-49-c4cdce9c3b-unnamed-chunk-171.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-50-c4cdce9c3b-unnamed-chunk-172.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-50-c4cdce9c3b-unnamed-chunk-173.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-51-c4cdce9c3b-unnamed-chunk-174.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-52-c4cdce9c3b-unnamed-chunk-175.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-52-c4cdce9c3b-unnamed-chunk-176.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-53-c4cdce9c3b-unnamed-chunk-177.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-54-c4cdce9c3b-unnamed-chunk-178.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-54-c4cdce9c3b-unnamed-chunk-179.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-55-c4cdce9c3b-unnamed-chunk-180.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-56-c4cdce9c3b-unnamed-chunk-181.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-57-c4cdce9c3b-unnamed-chunk-182.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-57-c4cdce9c3b-unnamed-chunk-183.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-58-c4cdce9c3b-unnamed-chunk-184.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-59-c4cdce9c3b-unnamed-chunk-185.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-52-59-c4cdce9c3b-unnamed-chunk-186.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-00-c4cdce9c3b-unnamed-chunk-187.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-01-c4cdce9c3b-unnamed-chunk-188.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-02-c4cdce9c3b-unnamed-chunk-189.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-02-c4cdce9c3b-unnamed-chunk-190.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-03-c4cdce9c3b-unnamed-chunk-191.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-04-c4cdce9c3b-unnamed-chunk-192.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-05-c4cdce9c3b-unnamed-chunk-193.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-05-c4cdce9c3b-unnamed-chunk-194.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-06-c4cdce9c3b-unnamed-chunk-195.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-07-c4cdce9c3b-unnamed-chunk-196.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-07-c4cdce9c3b-unnamed-chunk-197.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-08-c4cdce9c3b-unnamed-chunk-198.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-09-c4cdce9c3b-unnamed-chunk-199.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-10-c4cdce9c3b-unnamed-chunk-1100.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-10-c4cdce9c3b-unnamed-chunk-1101.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-11-c4cdce9c3b-unnamed-chunk-1102.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-12-c4cdce9c3b-unnamed-chunk-1103.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-12-c4cdce9c3b-unnamed-chunk-1104.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-13-c4cdce9c3b-unnamed-chunk-1105.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-14-c4cdce9c3b-unnamed-chunk-1106.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-14-c4cdce9c3b-unnamed-chunk-1107.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-15-c4cdce9c3b-unnamed-chunk-1108.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-16-c4cdce9c3b-unnamed-chunk-1109.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-17-c4cdce9c3b-unnamed-chunk-1110.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-17-c4cdce9c3b-unnamed-chunk-1111.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-18-c4cdce9c3b-unnamed-chunk-1112.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-19-c4cdce9c3b-unnamed-chunk-1113.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-20-c4cdce9c3b-unnamed-chunk-1114.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-20-c4cdce9c3b-unnamed-chunk-1115.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-21-c4cdce9c3b-unnamed-chunk-1116.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-22-c4cdce9c3b-unnamed-chunk-1117.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-22-c4cdce9c3b-unnamed-chunk-1118.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-23-c4cdce9c3b-unnamed-chunk-1119.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-24-c4cdce9c3b-unnamed-chunk-1120.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-25-c4cdce9c3b-unnamed-chunk-1121.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-25-c4cdce9c3b-unnamed-chunk-1122.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-26-c4cdce9c3b-unnamed-chunk-1123.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-27-c4cdce9c3b-unnamed-chunk-1124.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-27-c4cdce9c3b-unnamed-chunk-1125.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-28-c4cdce9c3b-unnamed-chunk-1126.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-29-c4cdce9c3b-unnamed-chunk-1127.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-27-16-53-30-c4cdce9c3b-unnamed-chunk-1128.png) 



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

![plot of chunk unnamed-chunk-2](http://carlboettiger.info/assets/figures/2012-12-27-16-53-31-c4cdce9c3b-unnamed-chunk-2.png) 

```r
yields
```

```
        method     V1      sd L1
1           GP 38.838 7.08073  1
2   Parametric 39.669 6.56923  1
3         True 39.669 6.56923  1
4   Structural  6.430 0.12192  1
5           GP 32.143 8.06603  2
6   Parametric 38.157 6.92022  2
7         True 39.669 6.56923  2
8   Structural 10.405 5.19610  2
9           GP 39.886 6.52492  3
10  Parametric 38.421 6.59104  3
11        True 39.669 6.56923  3
12  Structural  6.903 0.10410  3
13          GP 38.846 7.08771  4
14  Parametric 39.154 6.58619  4
15        True 39.669 6.56923  4
16  Structural  6.829 0.09727  4
17          GP 36.139 8.34198  5
18  Parametric 38.411 7.01027  5
19        True 39.669 6.56923  5
20  Structural  6.765 0.24963  5
21          GP 39.159 7.15022  6
22  Parametric 37.019 6.90501  6
23        True 39.669 6.56923  6
24  Structural  6.843 1.25255  6
25          GP 38.486 7.07333  7
26  Parametric 36.001 7.28714  7
27        True 39.669 6.56923  7
28  Structural  6.472 0.13923  7
29          GP 39.669 6.56923  8
30  Parametric 36.768 7.64312  8
31        True 39.669 6.56923  8
32  Structural  6.793 0.08866  8
33          GP 39.051 6.58684  9
34  Parametric 39.339 6.49418  9
35        True 39.669 6.56923  9
36  Structural  6.598 0.83902  9
37          GP 35.763 6.53735 10
38  Parametric 34.166 7.34005 10
39        True 39.669 6.56923 10
40  Structural  6.479 0.28512 10
41          GP 37.230 7.81937 11
42  Parametric 32.681 7.82089 11
43        True 39.669 6.56923 11
44  Structural  6.513 0.07441 11
45          GP 26.537 9.38809 12
46  Parametric 37.785 7.38231 12
47        True 39.669 6.56923 12
48  Structural  6.587 0.26338 12
49          GP 32.185 8.29100 13
50  Parametric 36.780 7.59566 13
51        True 39.669 6.56923 13
52  Structural  6.617 1.15091 13
53          GP 37.207 6.57849 14
54  Parametric 39.154 6.58619 14
55        True 39.669 6.56923 14
56  Structural  6.159 0.53734 14
57          GP 38.547 7.05936 15
58  Parametric 38.577 6.57756 15
59        True 39.669 6.56923 15
60  Structural  6.276 0.54672 15
61          GP 39.154 6.58619 16
62  Parametric 35.692 8.07688 16
63        True 39.669 6.56923 16
64  Structural  6.901 0.04340 16
65          GP 39.148 6.57939 17
66  Parametric  0.219 0.43350 17
67        True 39.669 6.56923 17
68  Structural  6.112 0.54788 17
69          GP 39.885 6.52190 18
70  Parametric 37.930 6.97607 18
71        True 39.669 6.56923 18
72  Structural  6.511 0.56195 18
73          GP 38.334 6.59949 19
74  Parametric 35.478 7.27250 19
75        True 39.669 6.56923 19
76  Structural  6.433 0.53889 19
77          GP 35.593 8.35485 20
78  Parametric  0.000 0.00000 20
79        True 39.669 6.56923 20
80  Structural  6.910 1.01935 20
81          GP 37.480 7.83792 21
82  Parametric 35.692 8.07688 21
83        True 39.669 6.56923 21
84  Structural  6.371 0.41498 21
85          GP 39.773 6.49764 22
86  Parametric 37.930 6.97607 22
87        True 39.669 6.56923 22
88  Structural  6.237 0.87235 22
89          GP 34.087 8.17287 23
90  Parametric 35.779 7.21839 23
91        True 39.669 6.56923 23
92  Structural  6.516 0.38807 23
93          GP 38.332 6.59553 24
94  Parametric 32.115 8.60811 24
95        True 39.669 6.56923 24
96  Structural  7.039 0.11615 24
97          GP 37.702 7.83639 25
98  Parametric 39.733 6.48536 25
99        True 39.669 6.56923 25
100 Structural  6.750 0.23762 25
101         GP 32.143 8.06603 26
102 Parametric 19.462 9.26487 26
103       True 39.669 6.56923 26
104 Structural  6.823 0.07225 26
105         GP 39.886 6.52492 27
106 Parametric 39.669 6.56923 27
107       True 39.669 6.56923 27
108 Structural  7.087 2.04830 27
109         GP 35.593 8.35485 28
110 Parametric 35.212 6.23004 28
111       True 39.669 6.56923 28
112 Structural  6.593 0.68838 28
113         GP 39.289 7.14690 29
114 Parametric 38.516 6.98267 29
115       True 39.669 6.56923 29
116 Structural  7.037 0.11508 29
117         GP 38.486 7.07333 30
118 Parametric 39.115 7.03491 30
119       True 39.669 6.56923 30
120 Structural  6.521 0.48346 30
121         GP 38.486 7.07333 31
122 Parametric 38.434 7.02288 31
123       True 39.669 6.56923 31
124 Structural  6.569 0.09078 31
125         GP 39.498 6.58566 32
126 Parametric 38.504 6.99314 32
127       True 39.669 6.56923 32
128 Structural  7.715 1.95599 32
```




NULL

