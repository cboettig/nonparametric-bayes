










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
sigma_m <- 0.0
x_grid <- seq(0, 1.5 * K, length=101)
h_grid <- x_grid
```


With parameters 1.1, 10, 5. 




```r
seed_i <- 1

  Xo <- K # observations start from
  x0 <- Xo # simulation under policy starts from
  obs <- sim_obs(Xo, z_g, f, p, Tobs=40, nz=1, 
                 harvest = sort(rep(seq(0, 1.8, length=8), 5)), seed = seed_i)
```

![plot of chunk unnamed-chunk-1](figure/2012-12-27-11-01-18-a49d5666d1-unnamed-chunk-11.png) 

```r
  alt <- par_est(obs)
  est <- par_est_allee(obs, f, p)
  gp <- bgp(X=obs$x, XX=x_grid, Z=obs$y, verb=0,
          meanfn="constant", bprior="b0", BTE=c(2000,16000,2),
          m0r1=FALSE, corr="exp", trace=FALSE, 
          beta = beta, s2.p = s2.p, d.p = d.p, nug.p = nug.p, tau2.p = tau2.p,
          s2.lam = "fixed", d.lam = "fixed", nug.lam = "fixed", tau2.lam = "fixed")      
  gp_plot(gp, f, p, est$f, est$p, alt$f, alt$p, x_grid, obs, seed_i)
```

![plot of chunk unnamed-chunk-1](figure/2012-12-27-11-01-18-a49d5666d1-unnamed-chunk-12.png) 

```r
#  posteriors_plot(gp, priors) # needs trace=TRUE!
  OPT <- optimal_policy(gp, f, est$f, alt$f,
                        p, est$p, alt$p,
                        x_grid, h_grid, sigma_g, 
                        est$sigma_g, alt$sigma_g, 
                        delta, xT, profit, reward, OptTime)
  plot_policies(x_grid, OPT$gp_D, OPT$est_D, OPT$true_D, OPT$alt_D)
```

![plot of chunk unnamed-chunk-1](figure/2012-12-27-11-01-18-a49d5666d1-unnamed-chunk-13.png) 

```r
  dt <- simulate_opt(OPT, f, p, x_grid, h_grid, x0, z_g, profit)
  sim_plots(dt, seed=seed_i)
```

![plot of chunk unnamed-chunk-1](figure/2012-12-27-11-01-18-a49d5666d1-unnamed-chunk-14.png) 

```r
  profits_stats(dt)
```

```
       method     V1    sd
1:         GP 10.719 1.472
2: Parametric  5.185 1.604
3:       True 12.006 1.723
4: Structural  8.154 1.929
```

  



```r
est
```

```
$f
function (x, h, p) 
{
    sapply(x, function(x) {
        x <- max(0, x - h)
        x * exp(p[1] * (1 - x/p[2]) * (x - p[3])/p[2])
    })
}
<environment: namespace:pdgControl>

$p
[1] 0.3448 9.2312 7.0138

$sigma_g
[1] 0.001
```

```r
alt
```

```
$f_alt
function (x, h, p) 
{
    sapply(x, function(x) {
        x <- max(0, x - h)
        max(0, x * exp(p[1] * (1 - x/p[2])))
    })
}
<environment: namespace:pdgControl>

$p_alt
      r       K 
 0.4935 12.0461 

$sigma_g_alt
    s 
0.001 
```




NULL

