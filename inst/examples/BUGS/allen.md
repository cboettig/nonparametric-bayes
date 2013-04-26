






```r
opts_chunk$set(tidy = FALSE, warning = FALSE, message = FALSE, cache = FALSE)
opts_knit$set(upload.fun = socialR::flickr.url)
theme_set(theme_bw(base_size = 10))
theme_update(panel.background = element_rect(fill = "transparent", 
    colour = NA), plot.background = element_rect(fill = "transparent", colour = NA))
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
    "#0072B2", "#D55E00", "#CC79A7")
```


### Model and parameters


```r
f <- RickerAllee
p <- c(2, 10, 5) 
K <- p[2]
allee <- p[3]
```



```r
sigma_g <- 0.05
sigma_m <- 0.0
z_g <- function() rlnorm(1, 0, sigma_g)
z_m <- function() 1+(2*runif(1, 0,  1)-1) * sigma_m
x_grid <- seq(0, 1.5 * K, length=101)
h_grid <- x_grid
profit <- function(x,h) pmin(x, h)
delta <- 0.01
OptTime <- 20  # stationarity with unstable models is tricky thing
reward <- 0
xT <- 0
seed_i <- 111
Xo <- K # observations start from
x0 <- Xo # simulation under policy starts from
Tobs <- 35
```


### Sample Data


```r
obs <- sim_obs(Xo, z_g, f, p, Tobs=Tobs, nz=15, 
                harvest = sort(rep(seq(0, .5, length=7), 5)), seed = seed_i)
```



## Maximum Likelihood


```r
alt <- par_est(obs,  init = c(r = p[1], 
                              K = mean(obs$x[obs$x>0]), 
                              s = sigma_g))
est <- par_est_allee(obs, f, p,  
                     init = c(r = p[1] + 1, 
                              K = p[2] + 2, 
                              C = p[3] + 2, 
                              s = sigma_g))
```



## Non-parametric Bayes



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
  gp <- bgp(X=obs$x, XX=x_grid, Z=obs$y, verb=0,
          meanfn="constant", bprior="b0", BTE=c(2000,16000,2),
          m0r1=FALSE, corr="exp", trace=TRUE, 
          beta = beta, s2.p = s2.p, d.p = d.p, nug.p = nug.p, tau2.p = tau2.p,
          s2.lam = "fixed", d.lam = "fixed", nug.lam = "fixed", tau2.lam = "fixed")      
  gp_plot(gp, f, p, est$f, est$p, alt$f, alt$p, x_grid, obs, seed_i)
```

![plot of chunk unnamed-chunk-2](http://farm9.staticflickr.com/8115/8682465554_b2015d36b9_o.png) 



```r
  posteriors_plot(gp, priors) # needs trace=TRUE!
```

![plot of chunk unnamed-chunk-3](http://farm9.staticflickr.com/8116/8682465736_0fe1557e5b_o.png) 




#### Parametric Bayes



```r
library(R2jags)
```



```r
# a bit unfair to start with the correct values, but anyhow...
init_p = p
names(init_p) = c("r0", "K", "theta")
y <- obs$y[-1] 
N=length(y);
```



```r
jags.data <- list("N","y")
n.chains = 1
n.iter = 20000
n.burnin = floor(n.iter/2)
n.thin = max(1, floor(n.chains * (n.iter - n.burnin)/1000))
```





## Gamma priors on precision


```r
jags.params=c("K","logr0","logtheta","iR","iQ") # Don't need to save "x"
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
## Compiling model graph
##    Resolving undeclared variables
##    Allocating nodes
##    Graph Size: 394
## 
## Initializing model
```

```r
time <- unname(time["elapsed"]);
```



#### Convergence diagnostics for parametric bayes


```r
library(emdbook)
library(coda)
tfit_jags_m <- as.mcmc.bugs(jagsfit$BUGSoutput)
print(xyplot(tfit_jags_m))
```

![plot of chunk unnamed-chunk-8](http://farm9.staticflickr.com/8524/8682466588_db50662d7d_o.png) 

```r
print(densityplot(tfit_jags_m))
```

![plot of chunk unnamed-chunk-8](http://farm9.staticflickr.com/8396/8681355375_39e73ba28b_o.png) 




```r
mcmc <- as.mcmc(jagsfit)
mcmcall <- mcmc[,-2]
who <- colnames(mcmcall)
who 
```

```
## [1] "K"        "iQ"       "iR"       "logr0"    "logtheta"
```

```r
mcmcall <- cbind(mcmcall[,1],mcmcall[,2],mcmcall[,3],mcmcall[,4],mcmcall[,5])
colnames(mcmcall) <- who
```





```r
MaxT = 1000
matrices_gp <- gp_transition_matrix(gp$ZZ.km, gp$ZZ.ks2, x_grid, h_grid)
opt_gp <- value_iteration(matrices_gp, x_grid, h_grid, MaxT, xT, profit, delta, reward)

matrices_true <- f_transition_matrix(f, p, x_grid, h_grid, sigma_g)
opt_true <- value_iteration(matrices_true, x_grid, h_grid, OptTime=MaxT, xT, profit, delta=delta)

matrices_estimated <- f_transition_matrix(est$f, est$p, x_grid, h_grid, est$sigma_g)
opt_estimated <- value_iteration(matrices_estimated, x_grid, h_grid, OptTime=MaxT, xT, profit, delta=delta)

matrices_alt <- f_transition_matrix(alt$f, alt$p, x_grid, h_grid, alt$sigma_g)
opt_alt <- value_iteration(matrices_alt, x_grid, h_grid, OptTime=MaxT, xT, profit, delta=delta)


pardist <- mcmcall
pardist[,4] = exp(pardist[,4])
pardist[,5] = exp(pardist[,5])


matrices_par_bayes <- parameter_uncertainty_SDP(f, p, x_grid, h_grid, pardist)
opt_par_bayes <- value_iteration(matrices_par_bayes, x_grid, h_grid, OptTime=MaxT, xT, profit, delta=delta)


OPT = list(gp_D = opt_gp$D, true_D = opt_true$D, est_D = opt_estimated$D, alt_D = opt_alt$D, par_bayes = opt_par_bayes$D)
```





```r
policies <- melt(data.frame(stock=x_grid, 
                            GP = x_grid[opt_gp$D], 
                            MLE = x_grid[opt_estimated$D],
                            True = x_grid[opt_true$D],
                            Ricker_MLE = x_grid[opt_alt$D],
                            Parametric.Bayes = x_grid[opt_par_bayes$D]),
                   id="stock")
names(policies) <- c("stock", "method", "value")
policy_plot <- ggplot(policies, aes(stock, stock - value, color=method)) +
  geom_line(lwd=1.2, alpha=0.8) + 
  xlab("stock size") + ylab("escapement")  +
  scale_colour_manual(values=cbPalette)
```





```r
sims <- lapply(OPT, function(D){
  set.seed(1)
  lapply(1:100, function(i) 
    ForwardSimulate(f, p, x_grid, h_grid, x0, D, z_g, profit=profit, OptTime=OptTime)
  )
})

  
dat <- melt(sims, id=names(sims[[1]][[1]]))
dt <- data.table(dat)
setnames(dt, c("L1", "L2"), c("method", "reps")) 
setkey(dt, method) # change the ordering


ggplot(dt) + 
    geom_line(aes(time, fishstock, group=interaction(reps,method), color=method), alpha=.1) +
    scale_colour_manual(values=cbPalette, guide = guide_legend(override.aes = list(alpha = 1)))
```

![plot of chunk unnamed-chunk-12](http://farm9.staticflickr.com/8400/8682482276_1a75ac7286_o.png) 

```r

ggplot(dt) + 
    geom_line(aes(time, harvest, group=interaction(reps,method), color=method), alpha=.1) +
    scale_colour_manual(values=cbPalette, guide = guide_legend(override.aes = list(alpha = 1)))
```

![plot of chunk unnamed-chunk-12](http://farm9.staticflickr.com/8382/8682482350_ef0a4ebfb3_o.png) 

```r

```


