# Comparison of Nonparametric Bayesian Gaussian Process estimates to standard the Parametric Bayesian approach

Plotting and knitr options, (can generally be ignored)






```r
posterior.mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
```





### Model and parameters

Uses the model derived in <a href="http://dx.doi.org/10.1080/10236190412331335373">Allen _et. al._ (2005)</a>, of a Ricker-like growth curve with an allee effect, defined in the pdgControl package,



```r
f <- Myers
p <- c(1, 2, 6)
K <- 5  # approx, a li'l' less
allee <- 1.2 # approx, a li'l' less
```



Various parameters defining noise dynamics, grid, and policy costs.  

Various parameters defining noise dynamics, grid, and policy costs.  


```r
sigma_g <- 0.05
sigma_m <- 0.0
z_g <- function() rlnorm(1, 0, sigma_g)
z_m <- function() 1+(2*runif(1, 0,  1)-1) * sigma_m
x_grid <- seq(0, 1.5 * K, length=50)
h_grid <- x_grid
profit <- function(x,h) pmin(x, h)
delta <- 0.01
OptTime <- 50  # stationarity with unstable models is tricky thing
reward <- 0
xT <- 0
Xo <- K # observations start from
x0 <- Xo # simulation under policy starts from
Tobs <- 40
```


### Sample Data


```r
  set.seed(1234)
  #harvest <- sort(rep(seq(0, .5, length=7), 5))
  x <- numeric(Tobs)
  x[1] <- 2
  nz <- 1
  for(t in 1:(Tobs-1))
    x[t+1] = z_g() * f(x[t], h=0, p=p)
  obs <- data.frame(x = c(rep(0,nz), 
                          pmax(rep(0,Tobs-1), x[1:(Tobs-1)])), 
                    y = c(rep(0,nz), 
                          x[2:Tobs]))
raw_plot <- ggplot(data.frame(time = 1:Tobs, x=x), aes(time,x)) + geom_line()
raw_plot
```

![plot of chunk obs](http://farm8.staticflickr.com/7346/8779148472_8920b6e922_o.png) 



## Maximum Likelihood


```r
est <- par_est_allee(obs, f, p,  
                     init = c(r = p[1] + 1, 
                              K = p[2] + 2, 
                              C = p[3] + 2, 
                              s = sigma_g))
```



Mean predictions


```r
true_means <- sapply(x_grid, f, 0, p)
est_means <- sapply(x_grid, est$f, 0, est$p)
```



## Non-parametric Bayes



```r
#inv gamma has mean b / (a - 1) (assuming a>1) and variance b ^ 2 / ((a - 2) * (a - 1) ^ 2) (assuming a>2)
s2.p <- c(5,5)  
d.p = c(10, 1/0.1)
```



Estimate the Gaussian Process (nonparametric Bayesian fit)


```r
gp <- gp_mcmc(obs$x, y=obs$y, n=1e5, s2.p = s2.p, d.p = d.p)
gp_dat <- gp_predict(gp, x_grid, burnin=1e4, thin=300)
```



Show traces and posteriors against priors


```r
plots <- summary_gp_mcmc(gp)
```

![plot of chunk gp_traces_densities](http://farm3.staticflickr.com/2892/8772628717_025a6c8797_o.png) ![plot of chunk gp_traces_densities](http://farm8.staticflickr.com/7420/8779192054_1936b81d6e_o.png) 




```r
# Summarize the GP model
tgp_dat <- 
    data.frame(  x = x_grid, 
                 y = gp_dat$E_Ef, 
                 ymin = gp_dat$E_Ef - 2 * sqrt(gp_dat$E_Vf), 
                 ymax = gp_dat$E_Ef + 2 * sqrt(gp_dat$E_Vf) )
```



Parametric Bayesian Models
---------------------------

We use the JAGS Gibbs sampler, a recent open source BUGS
implementation with an R interface that works on most platforms.
We initialize the usual MCMC parameters; see `?jags` for details.

All parametric Bayesian estimates use the following basic parameters for the JAGS MCMC:



```r
y <- obs$x[-1] 
N <- length(y);
jags.data <- list("N","y")
n.chains <- 1
n.iter <- 40000
n.burnin <- floor(10000)
n.thin <- max(1, floor(n.chains * (n.iter - n.burnin)/1000))
```



We will use the same priors for process and observation noise in each model, 


```r
stdQ_prior_p <- c(0.0001, 100)
stdR_prior_p <- c(0.0001, 100)
stdQ_prior  <- function(x) dunif(x, stdQ_prior_p[1], stdQ_prior_p[2])
stdR_prior  <- function(x) dunif(x, stdR_prior_p[1], stdR_prior_p[2])
```


### Parametric Bayes of correct (Allen) model

We initiate the MCMC chain (`init_p`) using the true values of the
parameters `p` from the simulation.  While impossible in real data, this
gives the parametric Bayesian approach the best chance at succeeding.
`y` is the timeseries (recall `obs` has the $x_t$, $x_{t+1}$ pairs)


```r
# a bit unfair to start with the correct values, but anyhow...
init_p = p
names(init_p) = c("r0", "K", "theta")
```


The actual model is defined in a `model.file` that contains an R function
that is automatically translated into BUGS code by *R2WinBUGS*.  The file
defines the priors and the model. We write the file from R as follows: 



```r
K_prior_p <- c(0.01, 40.0)
logr0_prior_p <- c(-6.0, 6.0)
logtheta_prior_p <- c(-6.0, 6.0)

bugs.model <- 
paste(sprintf(
"model{
  K     ~ dunif(%s, %s)
  logr0    ~ dunif(%s, %s)
  logtheta ~ dunif(%s, %s)
  stdQ ~ dunif(%s, %s)
  stdR ~ dunif(%s, %s)", 
  K_prior_p[1], K_prior_p[2],
  logr0_prior_p[1], logr0_prior_p[2],
  logtheta_prior_p[1], logtheta_prior_p[2],
  stdQ_prior_p[1], stdQ_prior_p[2],
  stdR_prior_p[1], stdR_prior_p[2]),

  "
  iQ <- 1 / (stdQ * stdQ);
  iR <- 1 / (stdR * stdR);
  r0 <- exp(logr0)
  theta <- exp(logtheta)
  x[1] ~ dunif(0, 10)
  for(t in 1:(N-1)){
    mu[t] <- x[t] * exp(r0 * (1 - x[t]/K)* (x[t] - theta) / K )
    x[t+1] ~ dnorm(mu[t], iQ) 
  }
  for(t in 1:(N)){
    y[t] ~ dnorm(x[t], iR)
  }
}")
writeLines(bugs.model, "allen.bugs")
```


Write the priors into a list for later reference


```r
K_prior     <- function(x) dunif(x, K_prior_p[1], K_prior_p[2])
logr0_prior <- function(x) dunif(x, logr0_prior_p[1], logr0_prior_p[2])
logtheta_prior <- function(x) dunif(x, -6, 6)
par_priors  <- list(K = K_prior, deviance = function(x) 0 * x, 
                    logr0 = logr0_prior, logtheta = logtheta_prior,
                    stdQ = stdQ_prior, stdR = stdR_prior)
```


We define which parameters to keep track of, and set the initial values of
parameters in the transformed space used by the MCMC.  We use logarithms
to maintain strictly positive values of parameters where appropriate.



```r
jags.params=c("K","logr0","logtheta","stdQ", "stdR")
jags.inits <- function(){
  list("K"=init_p["K"],"logr0"=log(init_p["r0"]),"logtheta"=log(init_p["theta"]), 
       "stdQ"=sqrt(0.05),"stdR"=sqrt(0.1),"x"=y,
       .RNG.name="base::Wichmann-Hill", .RNG.seed=123)
}

set.seed(12345)
jagsfit <- jags(data=jags.data, inits=jags.inits, jags.params, n.chains=n.chains, 
                n.iter=n.iter, n.thin=n.thin, n.burnin=n.burnin, 
                model.file="allen.bugs", progress.bar="none")
```

```
Compiling model graph
   Resolving undeclared variables
   Allocating nodes
   Graph Size: 365

Initializing model
```



#### Convergence diagnostics for parametric bayes


```r
jags_matrix <- as.data.frame(as.mcmc.bugs(jagsfit$BUGSoutput))
allen_posteriors <- melt(cbind(index = 1:dim(jags_matrix)[1], jags_matrix), id = "index")
ggplot(allen_posteriors) + geom_line(aes(index, value)) + facet_wrap(~ variable, scale="free", ncol=1)
```

![plot of chunk allen-traces](http://farm8.staticflickr.com/7313/8772645113_1b1858446e_o.png) 




```r
allen_priors <- ddply(allen_posteriors, "variable", function(dd){
    grid <- seq(min(dd$value), max(dd$value), length = 100) 
    data.frame(value = grid, density = par_priors[[dd$variable[1]]](grid))
})

ggplot(allen_posteriors, aes(value)) + 
  stat_density(geom="path", position="identity", alpha=0.7) +
  geom_line(data=allen_priors, aes(x=value, y=density), col="red") + 
  facet_wrap(~ variable, scale="free", ncol=3)
```

![plot of chunk allen-posteriors](http://farm4.staticflickr.com/3707/8779208996_9caab81725_o.png) 



Reshape the posterior parameter distribution data, transform back into original space, and calculate the mean parameters and mean function


```r
pardist <- acast(allen_posteriors, index ~ variable, subset=.(variable != "deviance") )
pardist[,2] = exp(pardist[,2]) # transform model parameters back first
pardist[,3] = exp(pardist[,3])
bayes_coef <- apply(pardist,2, posterior.mode) 
bayes_pars <- unname(c(bayes_coef[2], bayes_coef[1], bayes_coef[3])) # parameters formatted for f
allen_means <- sapply(x_grid, f, 0, bayes_pars)
```







## Parametric Bayes based on the structurally wrong model (Ricker)




```r
K_prior_p <- c(0.01, 40.0)
logr0_prior_p <- c(-6.0, 6.0)

bugs.model <- 
paste(sprintf(
"model{
  K    ~ dunif(%s, %s)
  logr0    ~ dunif(%s, %s)
  stdQ ~ dunif(%s, %s)
  stdR ~ dunif(%s, %s)", 
  K_prior_p[1], K_prior_p[2],
  logr0_prior_p[1], logr0_prior_p[2],
  stdQ_prior_p[1], stdQ_prior_p[2],
  stdR_prior_p[1], stdR_prior_p[2]),

  "
  iQ <- 1 / (stdQ * stdQ);
  iR <- 1 / (stdR * stdR);
  r0 <- exp(logr0)
  x[1] ~ dunif(0, 10)
  for(t in 1:(N-1)){
    mu[t] <- x[t] * exp(r0 * (1 - x[t]/K) )
    x[t+1] ~ dnorm(mu[t], iQ) 
  }
  for(t in 1:(N)){
    y[t] ~ dnorm(x[t], iR)
  }
}")
writeLines(bugs.model, "ricker.bugs")
```



Compute prior curves


```r
K_prior     <- function(x) dunif(x, K_prior_p[1], K_prior_p[2])
logr0_prior <- function(x) dunif(x, logr0_prior_p[1], logr0_prior_p[2])
par_priors <- list(K = K_prior, deviance = function(x) 0 * x, 
                   logr0 = logr0_prior, stdQ = stdQ_prior, stdR = stdR_prior)
```


We define which parameters to keep track of, and set the initial values of
parameters in the transformed space used by the MCMC.  We use logarithms
to maintain strictly positive values of parameters where appropriate.



```r
# Uniform priors on standard deviation terms
jags.params=c("K","logr0", "stdQ", "stdR")
jags.inits <- function(){
  list("K"=init_p["K"], "logr0"=log(init_p["r0"]), 
       "stdQ"=sqrt(0.05), "stdR"=sqrt(0.1), "x"=y, 
       .RNG.name="base::Wichmann-Hill", .RNG.seed=123)
}
set.seed(12345) 
jagsfit <- jags(data=jags.data, inits=jags.inits, jags.params, n.chains=n.chains, 
                  n.iter=n.iter, n.thin=n.thin, n.burnin=n.burnin,
                  model.file="ricker.bugs", progress.bar="none")
```

```
Compiling model graph
   Resolving undeclared variables
   Allocating nodes
   Graph Size: 287

Initializing model
```



#### Convergence diagnostics for parametric bayes Ricker model


```r
jags_matrix <- as.data.frame(as.mcmc.bugs(jagsfit$BUGSoutput))
ricker_posteriors <- melt(cbind(index = 1:dim(jags_matrix)[1], jags_matrix), id = "index")
ggplot(ricker_posteriors) + geom_line(aes(index, value)) + 
  facet_wrap(~ variable, scale="free", ncol=1)
```

![plot of chunk ricker_traces](http://farm3.staticflickr.com/2873/8779218746_08f33c3e5b_o.png) 



```r
ricker_priors <- ddply(ricker_posteriors, "variable", function(dd){
    grid <- seq(min(dd$value), max(dd$value), length = 100) 
    data.frame(value = grid, density = par_priors[[dd$variable[1]]](grid))
})
# plot posterior distributions
ggplot(ricker_posteriors, aes(value)) + 
  stat_density(geom="path", position="identity", alpha=0.7) +
  geom_line(data=ricker_priors, aes(x=value, y=density), col="red") + 
  facet_wrap(~ variable, scale="free", ncol=2)
```

![plot of chunk ricker_posteriors](http://farm4.staticflickr.com/3675/8779220516_3f503bd4ed_o.png) 



Reshape posteriors data, transform back, calculate mode and corresponding function.  


```r
ricker_pardist <- acast(ricker_posteriors, index ~ variable, subset=.(variable != "deviance") )
ricker_pardist[,"logr0"] = exp(ricker_pardist[,"logr0"]) # transform model parameters back first
bayes_coef <- apply(ricker_pardist,2, posterior.mode) # much better estimates from mode then mean
ricker_bayes_pars <- unname(c(bayes_coef[2], bayes_coef[1]))
ricker_means <- sapply(x_grid, Ricker, 0, ricker_bayes_pars[c(1,2)])
```



## Myers Parametric Bayes


```r
logr0_prior_p <- c(-6.0, 6.0)
logtheta_prior_p <- c(-6.0, 6.0)
logK_prior_p <- c(-6.0, 6.0)

bugs.model <- 
paste(sprintf(
"model{
  logr0    ~ dunif(%s, %s)
  logtheta    ~ dunif(%s, %s)
  logK    ~ dunif(%s, %s)
  stdQ ~ dunif(%s, %s)
  stdR ~ dunif(%s, %s)", 
  logr0_prior_p[1], logr0_prior_p[2],
  logtheta_prior_p[1], logtheta_prior_p[2],
  logK_prior_p[1], logK_prior_p[2],
  stdQ_prior_p[1], stdQ_prior_p[2],
  stdR_prior_p[1], stdR_prior_p[2]),

  "
  iQ <- 1 / (stdQ * stdQ);
  iR <- 1 / (stdR * stdR);

  r0 <- exp(logr0)
  theta <- exp(logtheta)
  K <- exp(logK)

  x[1] ~ dunif(0, 10)

  for(t in 1:(N-1)){
    mu[t] <- r0 * pow(abs(x[t]), theta) / (1 + pow(abs(x[t]), theta) / K)
    x[t+1] ~ dnorm(mu[t], iQ) 
  }

  for(t in 1:(N)){
    y[t] ~ dnorm(x[t], iR)
  }
}")
writeLines(bugs.model, "myers.bugs")
```




```r
logK_prior     <- function(x) dunif(x, logK_prior_p[1], logK_prior_p[2])
logr_prior     <- function(x) dunif(x, logr0_prior_p[1], logr0_prior_p[2])
logtheta_prior <- function(x) dunif(x, logtheta_prior_p[1], logtheta_prior_p[2])
par_priors <- list( deviance = function(x) 0 * x, logK = logK_prior,
                    logr0 = logr_prior, logtheta = logtheta_prior, 
                    stdQ = stdQ_prior, stdR = stdR_prior)
```




```r
init_p = c(1, 1, 1)
names(init_p) = c("r0", "theta", "K")
jags.params=c("logr0", "logtheta", "logK", "stdQ", "stdR")
jags.inits <- function(){
  list("logr0"=log(init_p["r0"]), "logtheta" = log(init_p["theta"]),  
       "logK"=log(init_p["K"]), "stdQ"=sqrt(0.05), 
       "stdR"=sqrt(0.1),"x"=y,.RNG.name="base::Wichmann-Hill", .RNG.seed=123)
}
set.seed(12345)
jagsfit <- jags(data=jags.data, inits=jags.inits, jags.params, 
                  n.chains=n.chains, n.iter=n.iter, n.thin=n.thin,
                  n.burnin=n.burnin, model.file="myers.bugs", progress.bar="none")
```

```
Compiling model graph
   Resolving undeclared variables
   Allocating nodes
   Graph Size: 326

Initializing model


Deleting model
```

```
Error: Error in node x[39] Failure to calculate log density
```


Convergence diagnostics for parametric bayes


```r
jags_matrix <- as.data.frame(as.mcmc.bugs(jagsfit$BUGSoutput))
myers_posteriors <- melt(cbind(index = 1:dim(jags_matrix)[1], jags_matrix), id = "index")
ggplot(myers_posteriors) + geom_line(aes(index, value)) +
  facet_wrap(~ variable, scale="free", ncol=1)
```

![plot of chunk myers-traces](http://farm8.staticflickr.com/7445/8779222440_1c36f295c6_o.png) 





```r
par_prior_curves <- ddply(myers_posteriors, "variable", function(dd){
    grid <- seq(min(dd$value), max(dd$value), length = 100) 
    data.frame(value = grid, density = par_priors[[dd$variable[1]]](grid))
})

ggplot(myers_posteriors, aes(value)) + 
  stat_density(geom="path", position="identity", alpha=0.7) +
  geom_line(data=par_prior_curves, aes(x=value, y=density), col="red") + 
  facet_wrap(~ variable, scale="free", ncol=3)
```

![plot of chunk myers-posteriors](http://farm6.staticflickr.com/5446/8779224062_05746c11d7_o.png) 




```r
myers_pardist <- as.matrix(jags_matrix[2:6])
```

```
Error: undefined columns selected
```

```r
myers_pardist[,1] = exp(myers_pardist[,1]) # transform model parameters back first
myers_pardist[,2] = exp(myers_pardist[,2]) # transform model parameters back first
myers_pardist[,3] = exp(myers_pardist[,3]) # transform model parameters back first
colnames(myers_pardist) = c("K", "r0", "theta", "stdQ", "stdR")
bayes_coef <- apply(myers_pardist,2, posterior.mode) # much better estimates
myers_bayes_pars <- unname(c(bayes_coef[2], bayes_coef[3], bayes_coef[1]))
myers_means <- sapply(x_grid, Myer_harvest, 0, myers_bayes_pars)
```




### Phase-space diagram of the expected dynamics


```r
models <- data.frame(x=x_grid, GP=tgp_dat$y, True=true_means, 
                     MLE=est_means, Ricker=ricker_means, 
                     Allen = allen_means,
                     Myers = myers_means)

models <- melt(models, id="x")
names(models) <- c("x", "method", "value")
```




```r
plot_gp <- ggplot(tgp_dat) + geom_ribbon(aes(x,y,ymin=ymin,ymax=ymax), fill="gray80") +
    geom_line(data=models, aes(x, value, col=method), lwd=1, alpha=0.8) + 
    geom_point(data=obs, aes(x,y), alpha=0.8) + 
    xlab(expression(X[t])) + ylab(expression(X[t+1])) +
    scale_colour_manual(values=cbPalette) 
print(plot_gp)
```

![plot of chunk Figure1](http://farm3.staticflickr.com/2822/8772663395_cff4ed50cc_o.png) 





## Optimal policies by value iteration

Compute the optimal policy under each model using stochastic dynamic programming. We begin with the policy based on the GP model,



```r
MaxT = 1000
# uses expected values from GP, instead of integrating over posterior
#matrices_gp <- gp_transition_matrix(gp_dat$E_Ef, gp_dat$E_Vf, x_grid, h_grid)

# Integrate over posteriors 
matrices_gp <- gp_transition_matrix(gp_dat$Ef_posterior, gp_dat$Vf_posterior, x_grid, h_grid) 

# Solve the SDP using the GP-derived transition matrix
opt_gp <- value_iteration(matrices_gp, x_grid, h_grid, MaxT, xT, profit, delta, reward)
```



Determine the optimal policy based on the allen and MLE models


```r
matrices_true <- f_transition_matrix(f, p, x_grid, h_grid, sigma_g)
opt_true <- value_iteration(matrices_true, x_grid, h_grid, OptTime=MaxT, xT, profit, delta=delta)

matrices_estimated <- f_transition_matrix(est$f, est$p, x_grid, h_grid, est$sigma_g)
opt_estimated <- value_iteration(matrices_estimated, x_grid, h_grid, OptTime=MaxT, xT, profit, delta=delta)
```


Determine the optimal policy based on Bayesian Allen model


```r
allen_f <- function(x,h,p) unname(f(x,h,p[c(2, 1, 3)]))
matrices_allen <- parameter_uncertainty_SDP(allen_f, x_grid, h_grid, pardist, 4)
opt_allen <- value_iteration(matrices_allen, x_grid, h_grid, OptTime=MaxT, xT, profit, delta=delta)
```


Bayesian Ricker


```r
ricker_f <- function(x, h, p) Ricker(x, h, unname(p[c(2, 1)])) # defined by pdgControl 
matrices_ricker <- parameter_uncertainty_SDP(ricker_f, x_grid, h_grid, as.matrix(ricker_pardist), 3)
opt_ricker <- value_iteration(matrices_ricker, x_grid, h_grid, OptTime=MaxT, xT, profit, delta=delta)
```



Bayesian Myers model


```r
myers_f <- function(x,h,p) Myer_harvest(x, h, p[c(2, 3, 1)])
matrices_myers <- parameter_uncertainty_SDP(myers_f, x_grid, h_grid, as.matrix(myers_pardist), 4)
myers_alt <- value_iteration(matrices_myers, x_grid, h_grid, OptTime=MaxT, xT, profit, delta=delta)
```



Assemble the data


```r
OPT = data.frame(GP = opt_gp$D, True = opt_true$D, MLE = opt_estimated$D, Ricker = opt_ricker$D, Allen = opt_allen$D, Myers = myers_alt$D)
colorkey=cbPalette
names(colorkey) = names(OPT) 
```




## Graph of the optimal policies


```r
policies <- melt(data.frame(stock=x_grid, sapply(OPT, function(x) x_grid[x])), id="stock")
names(policies) <- c("stock", "method", "value")

ggplot(policies, aes(stock, stock - value, color=method)) +
  geom_line(lwd=1.2, alpha=0.8) + xlab("stock size") + ylab("escapement")  +
  scale_colour_manual(values=colorkey)
```

![plot of chunk Figure2](http://farm9.staticflickr.com/8119/8779384440_33f33207b2_o.png) 




## Simulate 100 realizations managed under each of the policies



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
# Legend in original ordering please, not alphabetical: 
dt$method = factor(dt$method, ordered=TRUE, levels=names(OPT))
```



```r
ggplot(dt) + 
  geom_line(aes(time, fishstock, group=interaction(reps,method), color=method), alpha=.1) +
  scale_colour_manual(values=colorkey, guide = guide_legend(override.aes = list(alpha = 1)))
```

![plot of chunk Figure3](http://farm8.staticflickr.com/7334/8779389508_7429c8a9dd_o.png) 




```r
Profit <- dt[, sum(profit), by=c("reps", "method")]
Profit[, mean(V1), by="method"]
```

```
   method    V1
1:     GP 31.06
2:   True 31.14
3:    MLE 14.94
4: Ricker 25.58
5:  Allen  5.00
6:  Myers  5.00
```



```r
ggplot(Profit, aes(V1)) + geom_histogram() + 
  facet_wrap(~method, scales = "free_y") + guides(legend.position = "none") + xlab("Total profit by replicate")
```

![plot of chunk totalprofits](http://farm8.staticflickr.com/7419/8772826363_db47abe5fe_o.png) 

