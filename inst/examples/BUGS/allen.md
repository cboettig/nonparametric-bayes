










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
## linking to JAGS 3.3.0
```

```
## module basemod loaded
```

```
## module bugs loaded
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

![plot of chunk unnamed-chunk-6](http://farm9.staticflickr.com/8531/8643861412_025a46058c_o.png) 

```r
print(densityplot(tfit_jags_m))
```

![plot of chunk unnamed-chunk-6](http://farm9.staticflickr.com/8261/8642765173_13aee46c8e_o.png) 





```r
mcmc <- as.mcmc(jagsfit)
mcmcall <- mcmc[, -2]
mcmcall <- cbind(mcmcall[, 1], mcmcall[, 2], mcmcall[, 3], mcmcall[, 
    4], mcmcall[, 5])
## get parameter estimates (and transform inverse variances back to SDs)
coef <- apply(mcmcall, 2, mean)
coef[2:3] = log(sqrt(1/coef[2:3]))
coef[4:5] = exp(coef[2:3])

coef.median <- apply(mcmcall, 2, median)
coef.median[2:3] = log(sqrt(1/coef.median[2:3]))
coef.median[4:5] = exp(coef.median[4:5])

dens <- apply(mcmcall, 2, density)
coef.mode <- sapply(dens, function(x) x$x[which.max(x$y)])
coef.mode[2:3] = log(sqrt(1/coef.mode[2:3]))
coef.mode[4:5] = exp(coef.mode[4:5])

ci.mc <- t(apply(mcmcall, 2, quantile, probs = c(0.025, 0.975)))
## OBS remember that upper and lower are reversed!
ci.mc[2:3, 2:1] = log(sqrt(1/ci.mc[2:3, ]))
ci.mc[4:5, ] = exp(ci.mc[4:5, ])

## reorder the stuff (K,Q,R,r0,theta)=>(theta r0 K Q R)
myorder = c(5, 4, 1, 2, 3)
coef = coef[myorder]
coef.median = coef.median[myorder]
coef.mode = coef.mode[myorder]
ci.mc = ci.mc[myorder, ]


convinfo = list(effsize = effectiveSize(mcmc), geweke = geweke.diag(mcmc))  ## MCMC convergence info
```

