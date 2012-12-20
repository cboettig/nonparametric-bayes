









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


With parameters `1, 2, 4.5`. 


```r
require(snowfall)
```

```
Loading required package: snowfall
```

```
Loading required package: snow
```

```r
sfInit(cpu=2)
```

```
snowfall 1.84 initialized: sequential execution, one CPU.
```

```r
sfExportAll()
```

```
sfExportAll() ignored in sequential mode.
```

```r
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
yields <- 
  sfLapply(1:20, function(j){
    
x_0_observed <- allee + x_grid[5]
seed <- round(runif(1) * 1e6)
set.seed(seed)
print(seed)


## @knitr sim-obs
Tobs <- 50
x <- numeric(Tobs)
x[1] <- x_0_observed
for(t in 1:(Tobs-1))
  x[t+1] = z_g() * f(x[t], h=0, p=p)
qplot(1:Tobs, x)

## @knitr lag-data
obs <- data.frame(x=c(0,x[1:(Tobs-1)]),y=c(0,x[2:Tobs]))

## @knitr par-est
estf <- function(p){
  mu <- log(obs$x) + p["r"]*(1-obs$x/p["K"])
  -sum(dlnorm(obs$y, mu, p["s"]), log=TRUE)
}
o <- optim(par = c(r=1,K=mean(x),s=1), estf, method="L", lower=c(1e-3,1e-3,1e-3))
f_alt <- Ricker
p_alt <- c(o$par['r'], o$par['K'])
sigma_g_alt <- o$par['s']


## @knitr gp-fit
gp <- bgp(X=obs$x, XX=x_grid, Z=obs$y, verb=0,
          meanfn="constant", bprior="b0", BTE=c(2000,16000,2),
          m0r1=FALSE, corr="exp", trace=TRUE, 
          beta = beta, s2.p = s2.p, d.p = d.p, nug.p = nug.p, tau2.p = tau2.p,
          s2.lam = "fixed", d.lam = "fixed", nug.lam = "fixed", tau2.lam = "fixed")


## @knitr gp-data
V <- gp$ZZ.ks2
Ef = gp$ZZ.km
tgp_dat <- data.frame(x   = gp$XX[[1]], 
                      y   = gp$ZZ.km, 
                      ymin = gp$ZZ.km - 1.96 * sqrt(gp$ZZ.ks2), 
                      ymax = gp$ZZ.km + 1.96 * sqrt(gp$ZZ.ks2))


## @knitr gp-plot
true <- sapply(x_grid, f, 0, p)
est <- sapply(x_grid, f_alt, 0, p_alt)
models <- data.frame(x=x_grid, GP=tgp_dat$y, Parametric=est, True=true)
models <- melt(models, id="x")
names(models) <- c("x", "method", "value")
gp_plot <- ggplot(tgp_dat)  + geom_ribbon(aes(x,y,ymin=ymin,ymax=ymax), fill="gray80") +
  geom_line(data=models, aes(x, value, col=method), lwd=2, alpha=0.8) + 
  geom_point(data=obs, aes(x,y), alpha=0.8) + 
  xlab(expression(X[t])) + ylab(expression(X[t+1])) +
  scale_colour_manual(values=cbPalette)
print(gp_plot)


## @knitr gp-posteriors
hyperparameters <- c("index", "s2", "tau2", "beta0", "nug", "d", "ldetK")
posteriors <- melt(gp$trace$XX[[1]][,hyperparameters], id="index")
priors <- list(s2 = s2_prior, tau2 = tau2_prior, beta0 = dnorm, nug = nug_prior, d = d_prior, ldetK = function(x) 0)
prior_curves <- ddply(posteriors, "variable", function(dd){
  grid <- seq(min(dd$value), max(dd$value), length = 100)
  data.frame(value = grid, density = priors[[dd$variable[1]]](grid))
})
posteriors_plot <- ggplot(posteriors) + 
  geom_histogram(aes(x=value, y=..density..), lwd=2) +
  geom_line(data=prior_curves, aes(x=value, y=density), col="red", lwd=2) +
  facet_wrap(~ variable, scale="free")
#print(posteriors_plot)

## @knitr gp-opt
matrices_gp <- gp_transition_matrix(Ef, V, x_grid, h_grid)
opt_gp <- find_dp_optim(matrices_gp, x_grid, h_grid, OptTime, xT, profit, delta, reward=reward)

## @knitr true-opt
matrices_true <- f_transition_matrix(f, p, x_grid, h_grid, sigma_g)
opt_true <- find_dp_optim(matrices_true, x_grid, h_grid, OptTime, xT, profit, delta=delta, reward = reward)

## @knitr est-opt
matrices_estimated <- f_transition_matrix(f_alt, p_alt, x_grid, h_grid, sigma_g_alt)
opt_estimated <- find_dp_optim(matrices_estimated, x_grid, h_grid, OptTime, xT, profit, delta=delta, reward = reward)

## @knitr policy_plot
policies <- melt(data.frame(stock=x_grid, 
                            GP = x_grid[opt_gp$D[,1]], 
                            Parametric = x_grid[opt_estimated$D[,1]],
                            True = x_grid[opt_true$D[,1]]),
                 id="stock")
names(policies) <- c("stock", "method", "value")
policy_plot <- ggplot(policies, aes(stock, stock - value, color=method)) +
  geom_line(lwd=2, alpha=0.8) + 
  xlab("stock size") + ylab("escapement")  +
  scale_colour_manual(values=cbPalette)
print(policy_plot)


## @knitr stationary_policy_only
m <- sapply(1:OptTime, function(i) opt_gp$D[,1])
opt_gp$D <- m
mm <- sapply(1:OptTime, function(i) opt_true$D[,1])
opt_true$D <- mm
mmm <- sapply(1:OptTime, function(i) opt_estimated$D[,1])
opt_estimated$D <- mmm

## @knitr simulate
set.seed(1)
sim_gp <- lapply(1:100, function(i) ForwardSimulate(f, p, x_grid, h_grid, K, opt_gp$D, z_g, profit=profit))
set.seed(1)
sim_true <- lapply(1:100, function(i) ForwardSimulate(f, p, x_grid, h_grid, K, opt_true$D, z_g, profit=profit))
set.seed(1)
sim_est <- lapply(1:100, function(i) ForwardSimulate(f, p, x_grid, h_grid, K, opt_estimated$D, z_g, profit=profit))


## @knitr tidy
dat <- list(GP = sim_gp, Parametric = sim_est, True = sim_true)
dat <- melt(dat, id=names(dat[[1]][[1]]))
dt <- data.table(dat)
setnames(dt, c("L1", "L2"), c("method", "reps")) 

## @knitr sim-fish
fish_plot <- ggplot(dt) + 
  geom_line(aes(time, fishstock, group=interaction(reps,method), color=method), alpha=.1) +
  scale_colour_manual(values=cbPalette, guide = guide_legend(override.aes = list(alpha = 1)))
print(fish_plot)

## @knitr sim-harvest
harvest_plot <- ggplot(dt) +
  geom_line(aes(time, harvest, group=interaction(reps,method), color=method), alpha=.1) +
  scale_colour_manual(values=cbPalette, guide = guide_legend(override.aes = list(alpha = 1))) +
  ggtitle(paste("rep", j, "seed", seed))
print(harvest_plot)

## @knitr costs
profits <- dt[, sum(profit), by = c("reps", "method")]
means <- profits[, mean(V1), by = method]
sds <- profits[, sd(V1), by = method]
yield <- cbind(means, sd = sds$V1)
yield

  yield
})
```

```
[1] 578227
```

![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-52-33-fd14d91468-unnamed-chunk-11.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-52-35-fd14d91468-unnamed-chunk-12.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-52-38-fd14d91468-unnamed-chunk-13.png) 

```
[1] 284002
```

![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-52-42-fd14d91468-unnamed-chunk-14.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-52-45-fd14d91468-unnamed-chunk-15.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-52-47-fd14d91468-unnamed-chunk-16.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-52-49-fd14d91468-unnamed-chunk-17.png) 

```
[1] 284002
```

![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-52-53-fd14d91468-unnamed-chunk-18.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-52-56-fd14d91468-unnamed-chunk-19.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-52-59-fd14d91468-unnamed-chunk-110.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-53-01-fd14d91468-unnamed-chunk-111.png) 

```
[1] 284002
```

![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-53-05-fd14d91468-unnamed-chunk-112.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-53-08-fd14d91468-unnamed-chunk-113.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-53-10-fd14d91468-unnamed-chunk-114.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-53-12-fd14d91468-unnamed-chunk-115.png) 

```
[1] 284002
```

![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-53-16-fd14d91468-unnamed-chunk-116.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-53-19-fd14d91468-unnamed-chunk-117.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-53-21-fd14d91468-unnamed-chunk-118.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-53-24-fd14d91468-unnamed-chunk-119.png) 

```
[1] 284002
```

![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-53-28-fd14d91468-unnamed-chunk-120.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-53-32-fd14d91468-unnamed-chunk-121.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-53-34-fd14d91468-unnamed-chunk-122.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-53-36-fd14d91468-unnamed-chunk-123.png) 

```
[1] 284002
```

![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-53-40-fd14d91468-unnamed-chunk-124.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-53-44-fd14d91468-unnamed-chunk-125.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-53-51-fd14d91468-unnamed-chunk-126.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-53-53-fd14d91468-unnamed-chunk-127.png) 

```
[1] 284002
```

![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-53-57-fd14d91468-unnamed-chunk-128.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-54-00-fd14d91468-unnamed-chunk-129.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-54-02-fd14d91468-unnamed-chunk-130.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-54-04-fd14d91468-unnamed-chunk-131.png) 

```
[1] 284002
```

![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-54-11-fd14d91468-unnamed-chunk-132.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-54-14-fd14d91468-unnamed-chunk-133.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-54-16-fd14d91468-unnamed-chunk-134.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-54-18-fd14d91468-unnamed-chunk-135.png) 

```
[1] 284002
```

![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-54-22-fd14d91468-unnamed-chunk-136.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-54-25-fd14d91468-unnamed-chunk-137.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-54-27-fd14d91468-unnamed-chunk-138.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-54-29-fd14d91468-unnamed-chunk-139.png) 

```
[1] 284002
```

![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-54-33-fd14d91468-unnamed-chunk-140.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-54-36-fd14d91468-unnamed-chunk-141.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-54-39-fd14d91468-unnamed-chunk-142.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-54-41-fd14d91468-unnamed-chunk-143.png) 

```
[1] 284002
```

![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-54-45-fd14d91468-unnamed-chunk-144.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-54-48-fd14d91468-unnamed-chunk-145.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-54-50-fd14d91468-unnamed-chunk-146.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-54-52-fd14d91468-unnamed-chunk-147.png) 

```
[1] 284002
```

![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-55-01-fd14d91468-unnamed-chunk-148.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-55-04-fd14d91468-unnamed-chunk-149.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-55-06-fd14d91468-unnamed-chunk-150.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-55-09-fd14d91468-unnamed-chunk-151.png) 

```
[1] 284002
```

![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-55-13-fd14d91468-unnamed-chunk-152.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-55-22-fd14d91468-unnamed-chunk-153.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-55-24-fd14d91468-unnamed-chunk-154.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-55-26-fd14d91468-unnamed-chunk-155.png) 

```
[1] 284002
```

![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-55-31-fd14d91468-unnamed-chunk-156.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-55-34-fd14d91468-unnamed-chunk-157.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-55-36-fd14d91468-unnamed-chunk-158.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-55-38-fd14d91468-unnamed-chunk-159.png) 

```
[1] 284002
```

![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-55-42-fd14d91468-unnamed-chunk-160.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-55-45-fd14d91468-unnamed-chunk-161.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-55-47-fd14d91468-unnamed-chunk-162.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-55-50-fd14d91468-unnamed-chunk-163.png) 

```
[1] 284002
```

![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-55-54-fd14d91468-unnamed-chunk-164.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-55-57-fd14d91468-unnamed-chunk-165.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-55-59-fd14d91468-unnamed-chunk-166.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-56-01-fd14d91468-unnamed-chunk-167.png) 

```
[1] 284002
```

![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-56-05-fd14d91468-unnamed-chunk-168.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-56-09-fd14d91468-unnamed-chunk-169.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-56-11-fd14d91468-unnamed-chunk-170.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-56-13-fd14d91468-unnamed-chunk-171.png) 

```
[1] 284002
```

![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-56-18-fd14d91468-unnamed-chunk-172.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-56-21-fd14d91468-unnamed-chunk-173.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-56-23-fd14d91468-unnamed-chunk-174.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-56-25-fd14d91468-unnamed-chunk-175.png) 

```
[1] 284002
```

![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-56-29-fd14d91468-unnamed-chunk-176.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-56-32-fd14d91468-unnamed-chunk-177.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-56-34-fd14d91468-unnamed-chunk-178.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-56-36-fd14d91468-unnamed-chunk-179.png) ![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-20-13-56-40-fd14d91468-unnamed-chunk-180.png) 



```r
yields <- melt(yields, id=c("method", "V1", "sd"))
yields         
```

```
       method    V1     sd L1
1          GP 2.819 0.4465  1
2  Parametric 1.868 0.1102  1
3        True 3.114 0.4534  1
4          GP 3.149 0.4574  2
5  Parametric 1.755 0.1180  2
6        True 3.114 0.4534  2
7          GP 3.149 0.4574  3
8  Parametric 1.755 0.1180  3
9        True 3.114 0.4534  3
10         GP 3.149 0.4574  4
11 Parametric 1.755 0.1180  4
12       True 3.114 0.4534  4
13         GP 3.149 0.4574  5
14 Parametric 1.755 0.1180  5
15       True 3.114 0.4534  5
16         GP 3.149 0.4574  6
17 Parametric 1.755 0.1180  6
18       True 3.114 0.4534  6
19         GP 3.149 0.4574  7
20 Parametric 1.755 0.1180  7
21       True 3.114 0.4534  7
22         GP 3.149 0.4574  8
23 Parametric 1.755 0.1180  8
24       True 3.114 0.4534  8
25         GP 3.149 0.4574  9
26 Parametric 1.755 0.1180  9
27       True 3.114 0.4534  9
28         GP 3.149 0.4574 10
29 Parametric 1.755 0.1180 10
30       True 3.114 0.4534 10
31         GP 3.149 0.4574 11
32 Parametric 1.755 0.1180 11
33       True 3.114 0.4534 11
34         GP 3.149 0.4574 12
35 Parametric 1.755 0.1180 12
36       True 3.114 0.4534 12
37         GP 3.149 0.4574 13
38 Parametric 1.755 0.1180 13
39       True 3.114 0.4534 13
40         GP 3.149 0.4574 14
41 Parametric 1.755 0.1180 14
42       True 3.114 0.4534 14
43         GP 3.149 0.4574 15
44 Parametric 1.755 0.1180 15
45       True 3.114 0.4534 15
46         GP 3.149 0.4574 16
47 Parametric 1.755 0.1180 16
48       True 3.114 0.4534 16
49         GP 3.149 0.4574 17
50 Parametric 1.755 0.1180 17
51       True 3.114 0.4534 17
52         GP 3.149 0.4574 18
53 Parametric 1.755 0.1180 18
54       True 3.114 0.4534 18
55         GP 3.149 0.4574 19
56 Parametric 1.755 0.1180 19
57       True 3.114 0.4534 19
58         GP 3.149 0.4574 20
59 Parametric 1.755 0.1180 20
60       True 3.114 0.4534 20
```




<p>Myers RA, Barrowman NJ, Hutchings JA and Rosenberg AA (1995).
&ldquo;Population Dynamics of Exploited Fish Stocks at Low Population Levels.&rdquo;
<EM>Science</EM>, <B>269</B>.
ISSN 0036-8075, <a href="http://dx.doi.org/10.1126/science.269.5227.1106">http://dx.doi.org/10.1126/science.269.5227.1106</a>.

