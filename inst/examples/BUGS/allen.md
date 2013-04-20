










```r
f <- RickerAllee
p <- c(2, 10, 5)
K <- 10
allee <- 5
```



```r
sigma_g <- 0.05
sigma_m <- 0
z_g <- function() rlnorm(1, 0, sigma_g)
z_m <- function() 1 + (2 * runif(1, 0, 1) - 1) * sigma_m
x_grid <- seq(0, 1.5 * K, length = 101)
h_grid <- x_grid
profit <- function(x, h) pmin(x, h)
delta <- 0.01
OptTime <- 20  # stationarity with unstable models is tricky thing
reward <- 0
xT <- 0
seed_i <- 1
Xo <- K  # observations start from
x0 <- Xo  # simulation under policy starts from
Tobs <- 35
```



```r
obs <- sim_obs(Xo, z_g, f, p, Tobs = Tobs, nz = 1, harvest = sort(rep(seq(0, 
    0.5, length = 7), 5)), seed = seed_i)
```



```r
alt <- par_est(obs, init = c(r = p[1], K = mean(obs$x[obs$x > 0]), 
    s = sigma_g))
est <- par_est_allee(obs, f, p, init = c(r = p[1] + 1, K = p[2] + 
    2, C = p[3] + 2, s = sigma_g))
```






```r
library(R2jags)
```

```
## Loading required package: R2WinBUGS
```

```
## Loading required package: boot
```

```
## Attaching package: 'boot'
```

```
## The following object(s) are masked from 'package:lattice':
## 
## melanoma
```

```
## Loading required package: rjags
```

```
## Linked to JAGS 3.3.0
```

```
## Loaded modules: basemod,bugs
```

```
## Loading required package: abind
```

```
## Loading required package: parallel
```

```
## Attaching package: 'R2jags'
```

```
## The following object(s) are masked from 'package:coda':
## 
## traceplot
```



```r
init_p = est$p
names(init_p) = c("r0", "K", "theta")
y <- obs$y[-1]
N = length(y)
```



```r
jags.data <- list("N", "y")
n.chains = 1
n.iter = 20000
n.burnin = floor(n.iter/2)
n.thin = max(1, floor(n.chains * (n.iter - n.burnin)/1000))
```





## Gamma priors on precision


```r
jags.params=c("K","logr0","logtheta","iR","iQ"); # Don't need to save "x"
jags.inits <- function(){
  list("K"=init_p["K"],"logr0"=log(init_p["r0"]),"logtheta"=log(init_p["theta"]),"iQ"=1/0.05,"iR"=1/0.1,"x"=y,.RNG.name="base::Wichmann-Hill", .RNG.seed=123)
}
set.seed(12345)
time<-system.time(       
  jagsfit <- jags(data=jags.data, inits=jags.inits, jags.params, n.chains=n.chains, 
                  n.iter=n.iter, n.thin=n.thin, n.burnin=n.burnin,model.file="bugmodel-GammaPrior.txt")
)         
```

```
## module glm loaded
```

```
## Compiling model graph
##    Resolving undeclared variables
##    Allocating nodes
##    Graph Size: 282
## 
## Initializing model
## 
##   |                                                          |                                                  |   0%  |                                                          |++                                                |   4%  |                                                          |++++                                              |   8%  |                                                          |++++++                                            |  12%  |                                                          |++++++++                                          |  16%  |                                                          |++++++++++                                        |  20%  |                                                          |++++++++++++                                      |  24%  |                                                          |++++++++++++++                                    |  28%  |                                                          |++++++++++++++++                                  |  32%  |                                                          |++++++++++++++++++                                |  36%  |                                                          |++++++++++++++++++++                              |  40%  |                                                          |++++++++++++++++++++++                            |  44%  |                                                          |++++++++++++++++++++++++                          |  48%  |                                                          |++++++++++++++++++++++++++                        |  52%  |                                                          |++++++++++++++++++++++++++++                      |  56%  |                                                          |++++++++++++++++++++++++++++++                    |  60%  |                                                          |++++++++++++++++++++++++++++++++                  |  64%  |                                                          |++++++++++++++++++++++++++++++++++                |  68%  |                                                          |++++++++++++++++++++++++++++++++++++              |  72%  |                                                          |++++++++++++++++++++++++++++++++++++++            |  76%  |                                                          |++++++++++++++++++++++++++++++++++++++++          |  80%  |                                                          |++++++++++++++++++++++++++++++++++++++++++        |  84%  |                                                          |++++++++++++++++++++++++++++++++++++++++++++      |  88%  |                                                          |++++++++++++++++++++++++++++++++++++++++++++++    |  92%  |                                                          |++++++++++++++++++++++++++++++++++++++++++++++++  |  96%  |                                                          |++++++++++++++++++++++++++++++++++++++++++++++++++| 100%
##   |                                                          |                                                  |   0%  |                                                          |**                                                |   4%  |                                                          |****                                              |   8%  |                                                          |******                                            |  12%  |                                                          |********                                          |  16%  |                                                          |**********                                        |  20%  |                                                          |************                                      |  24%  |                                                          |**************                                    |  28%  |                                                          |****************                                  |  32%  |                                                          |******************                                |  36%  |                                                          |********************                              |  40%  |                                                          |**********************                            |  44%  |                                                          |************************                          |  48%  |                                                          |**************************                        |  52%  |                                                          |****************************                      |  56%  |                                                          |******************************                    |  60%  |                                                          |********************************                  |  64%  |                                                          |**********************************                |  68%  |                                                          |************************************              |  72%  |                                                          |**************************************            |  76%  |                                                          |****************************************          |  80%  |                                                          |******************************************        |  84%  |                                                          |********************************************      |  88%  |                                                          |**********************************************    |  92%  |                                                          |************************************************  |  96%  |                                                          |**************************************************| 100%
```

```r
time <- unname(time["elapsed"]);
```




```r
library(emdbook)
library(coda)
tfit_jags_m <- as.mcmc.bugs(jagsfit$BUGSoutput)
print(xyplot(tfit_jags_m))
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-61.png) 

```r
print(densityplot(tfit_jags_m))
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-62.png) 





```r
mcmc <- as.mcmc(jagsfit)
mcmcall <- mcmc[, -2]
who <- colnames(mcmcall)
who
```

```
## [1] "K"        "iQ"       "iR"       "logr0"    "logtheta"
```

```r
mcmcall <- cbind(mcmcall[, 1], mcmcall[, 2], mcmcall[, 3], mcmcall[, 
    4], mcmcall[, 5])
colnames(mcmcall) <- who
```





```r
theta <- exp(mcmcall[, "logtheta"])
theta_dist <- hist(theta, freq = FALSE)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8.png) 

```r
theta_dist$mids
```

```
##  [1] 0.25 0.75 1.25 1.75 2.25 2.75 3.25 3.75 4.25 4.75 5.25 5.75 6.25 6.75
## [15] 7.25
```

```r
theta_dist$density
```

```
##  [1] 1.154 0.256 0.130 0.092 0.068 0.058 0.044 0.042 0.042 0.018 0.022
## [12] 0.026 0.022 0.014 0.012
```

```r
delta <- theta_dist$mids[2] - theta_dist$mids[1]
```


evaluating the value function given $f$ fixed at each of `theta_dist$mids` multiplied by  `theta_dist$density * delta` and summed over each $\theta$ is going to be slow for a single parameter but simply unrealistic for higher dimensions.  Instead we will rely on Monte Carlo sampling of parameter values from their posteriors,


```r
mc_n <- 100
thetas <- sample(theta, mc_n, replace = T)
```



