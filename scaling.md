---
layout: page
---



```r
f <- RickerAllee
p <- c(.5, .8, .5)
K <- .8  # approx, a li'l' less
allee <- .5 # approx, a li'l' less
```



```r
sigma_g <- 0.005
sigma_m <- 0.0
z_g <- function() rlnorm(1, 0, sigma_g)
z_m <- function() 1
x_grid <- seq(0, 1.5 * K, length=50)
h_grid <- x_grid
profit <- function(x,h) pmin(x, h)
delta <- 0.01
OptTime <- 50  # stationarity with unstable models is tricky thing
reward <- 0
xT <- 0
Xo <-  allee+.2# observations start from
x0 <- K # simulation under policy starts from
Tobs <- 40
MaxT = 1000 # timeout for value iteration convergence
```



```r
set.seed(1234)
x <- numeric(Tobs)
x[1] <- Xo
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

![plot of chunk unnamed-chunk-4](http://farm6.staticflickr.com/5513/9488645423_86e44c3679_o.png) 



```r
set.seed(12345)
estf <- function(p){ 
  mu <- f(obs$x,0,p)
  -sum(dlnorm(obs$y, log(mu), p[4]), log=TRUE)
}
par <- c(p[1]*rlnorm(1,0,.1), 
         p[2]*rlnorm(1,0,.1), 
         p[3]*rlnorm(1,0, .1), 
         sigma_g * rlnorm(1,0,.1))
o <- optim(par, estf, method="L", lower=c(1e-5,1e-5,1e-5,1e-5))
f_alt <- f
p_alt <- c(as.numeric(o$par[1]), as.numeric(o$par[2]), as.numeric(o$par[3]))
sigma_g_alt <- as.numeric(o$par[4])

est <- list(f = f_alt, p = p_alt, sigma_g = sigma_g_alt, mloglik=o$value)


## @knitr mle-output
true_means <- sapply(x_grid, f, 0, p)
est_means <- sapply(x_grid, est$f, 0, est$p)


## @knitr gp-priors
#inv gamma has mean b / (a - 1) (assuming a>1) and variance b ^ 2 / ((a - 2) * (a - 1) ^ 2) (assuming a>2)
s2.p <- c(5,5)  
d.p <- c(1,1) # = c(10, 1/0.1)
```



```r
gp <- gp_mcmc(obs$x, y=obs$y, n=1e5, s2.p = s2.p, d.p = d.p)
gp_dat <- gp_predict(gp, x_grid, burnin=1e4, thin=300)
```





```r
## @knitr gp_traces_densities
gp_assessment_plots <- summary_gp_mcmc(gp, burnin=1e4, thin=300)
gp_assessment_plots
```

```
$traces_plot
```

![plot of chunk unnamed-chunk-7](http://farm4.staticflickr.com/3802/9491452730_34ef90a140_o.png) 

```

$posteriors_plot
```

![plot of chunk unnamed-chunk-7](http://farm3.staticflickr.com/2811/9488657493_ffd898069f_o.png) 



```r
## @knitr gp-output 
tgp_dat <- data.frame(x = x_grid, 
               y = gp_dat$E_Ef, 
               ymin = pmax(gp_dat$E_Ef - 2 * sqrt(gp_dat$E_Vf),0), 
               ymax = gp_dat$E_Ef + 2 * sqrt(gp_dat$E_Vf) )


## @knitr assemble-models
models <- data.frame(x=x_grid, 
                     GP=tgp_dat$y, 
                     True=true_means)
models <- melt(models, id="x")
names(models) <- c("x", "method", "value")


## @knitr Figure1
plot_gp <- ggplot(tgp_dat) + geom_ribbon(aes(x,y,ymin=ymin,ymax=ymax), fill="gray80") +
  geom_line(data=models, aes(x, value, col=method), lwd=1, alpha=0.8) + 
  geom_point(data=obs, aes(x,y), alpha=0.8) + 
  xlab(expression(X[t])) + ylab(expression(X[t+1])) 
print(plot_gp)
```

![plot of chunk unnamed-chunk-8](http://farm4.staticflickr.com/3681/9488658009_efdb474fe9_o.png) 



```r
dat_fig1c <- melt(lapply(sample(100), 
  function(i){
    data.frame(stock = x_grid, 
                GP = pmax( # non-negative curves only
                          mvrnorm(1, 
                                  gp_dat$Ef_posterior[,i], 
                                  gp_dat$E_Cf # could draw from Cf_posterior after assuring is pos-def
                                  )
                          ,0)
               )
    
  }), id="stock")

ggplot(dat_fig1c, aes(stock, value, group=L1)) + geom_line(alpha=.3)
```

![plot of chunk unnamed-chunk-9](http://farm3.staticflickr.com/2848/9491454504_4b8d995bd3_o.png) 

```r

# a <- dcast(dat_fig1c, stock ~ L1, value.var='value')
# dat_fig1c <- data.table(dat_fig1c)
# to_grid <- function(x, grid) sapply(x, function(x) which.min(abs(grid - x)))  
# dat_fig1c$stock <- to_grid(dat_fig1c$stock, x_grid)
# a <- dat_fig1c[, to_grid(value, x_grid), by=c('stock')]
```




```r
matrices_gp <- gp_transition_matrix(gp_dat$Ef_posterior, gp_dat$Vf_posterior, x_grid, h_grid) 
opt_gp <- value_iteration(matrices_gp, x_grid, h_grid, MaxT, xT, profit, delta)

opt_gp <- find_dp_optim(matrices_gp, x_grid, h_grid, 5, xT=0.5, profit, delta, reward=1)
opt_gp$D <- opt_gp$D[,1]

## @knitr mle-opt
matrices_true <- f_transition_matrix(f, p, x_grid, h_grid, sigma_g)
opt_true <- value_iteration(matrices_true, x_grid, h_grid, OptTime=MaxT, xT, profit, delta=.01)
```





```r
OPT = data.frame(GP = opt_gp$D, True = opt_true$D)

## @knitr Figure2
policies <- melt(data.frame(stock=x_grid, sapply(OPT, function(x) x_grid[x])), id="stock")
names(policies) <- c("stock", "method", "value")

ggplot(policies, aes(stock, stock - value, color=method)) +
  geom_line(lwd=1.2, alpha=0.8) + xlab("stock size") + ylab("escapement") 
```

![plot of chunk unnamed-chunk-11](http://farm4.staticflickr.com/3683/9488670189_2e7168e731_o.png) 





```r
## @knitr sims
sims <- lapply(OPT, function(D){
  set.seed(1)
  lapply(1:100, function(i) 
    ForwardSimulate(f, p, x_grid, h_grid, x0, D, z_g, profit=profit, OptTime=OptTime)
  )
})

dat <- melt(sims, id=names(sims[[1]][[1]]))
sims_data <- data.table(dat)
setnames(sims_data, c("L1", "L2"), c("method", "reps")) 
# Legend in original ordering please, not alphabetical: 
sims_data$method = factor(sims_data$method, ordered=TRUE, levels=names(OPT))


## @knitr Figure3
ggplot(sims_data) + 
  geom_line(aes(time, fishstock, group=interaction(reps,method), color=method), alpha=.1) 
```

![plot of chunk unnamed-chunk-12](http://farm4.staticflickr.com/3728/9488671349_5612a1de51_o.png) 




```r
## @knitr profits
Profit <- sims_data[, sum(profit), by=c("reps", "method")]
tmp <- dcast(Profit, reps ~ method)
tmp <- tmp / tmp[,"True"]
tmp <- melt(tmp[2:dim(tmp)[2]])
actual_over_optimal <-subset(tmp, variable != "True")


## @knitr Figure4 
ggplot(actual_over_optimal, aes(value)) + geom_histogram(aes(fill=variable)) + 
  facet_wrap(~variable, scales = "free_y")  + guides(legend.position = "none") +
  xlab("Total profit by replicate")
```

![plot of chunk unnamed-chunk-13](http://farm4.staticflickr.com/3761/9488672199_2355ecaa47_o.png) 











# Forecasts


## One-step ahead predictor plots and forecasting plots



Now we show how each of the fitted models performs on the training data (e.g. plot of the step-ahead predictors). For the GP, we need to predict on the training data first:



```r
gp_f_at_obs <- gp_predict(gp, x, burnin=1e4, thin=300)
```


For the parametric models, prediction is a matter of feeding in each of the observed data points to the fitted parameters, like so 



```r
step_ahead <- function(x, f, p){
  h = 0
  x_predict <- sapply(x, f, h, p)
  n <- length(x_predict) - 1
  y <- c(x[1], x_predict[1:n])
  y
}
```


Which we apply over each of the fitted models, including the GP, organizing the "expected" transition points (given the previous point) into a data frame.  



```r
df <- melt(data.frame(time = 1:length(x), stock = x, 
                GP = gp_f_at_obs$E_Ef,
                True = step_ahead(x,f,p)
                 ), id=c("time", "stock"))

ggplot(df) + geom_point(aes(time, stock)) + 
  geom_line(aes(time, value, col=variable)) 
```

![plot of chunk unnamed-chunk-16](http://farm8.staticflickr.com/7407/9488672745_58e2364a6a_o.png) 



## Posterior predictive curves

This shows only the mean predictions.  For the Bayesian cases, we can instead loop over the posteriors of the parameters (or samples from the GP posterior) to get the distribution of such curves in each case.  

We will need a vector version (`pmin` in place of `min`) of this function that can operate on the posteriors, others are vectorized already.


```r
ricker_f <- function(x,h,p){
  sapply(x, function(x){ 
    x <- pmax(0, x-h) 
    pmax(0, x * exp(p[2] * (1 - x / p[1] )) )
  })
}
```



Then we proceed as before, now looping over 100 random samples from the posterior for each Bayesian estimate.  We write this as a function for easy reuse.  


```r
require(MASS)

step_ahead_posteriors <- function(x){
gp_f_at_obs <- gp_predict(gp, x, burnin=1e4, thin=300)
df_post <- melt(lapply(sample(100), 
  function(i){
    data.frame(time = 1:length(x), stock = x, 
                GP = pmax(mvrnorm(1, 
                                  gp_f_at_obs$Ef_posterior[,i], 
                                  pmax(gp_f_at_obs$E_Cf, 1e-4)
                                  #pmax(gp_f_at_obs$Cf_posterior[[i]],1e-5)
                                  ),0),
                True = step_ahead(x,f,p))
  }), id=c("time", "stock"))

}


df_post <- step_ahead_posteriors(x)

ggplot(df_post) + geom_point(aes(time, stock)) + 
  geom_line(aes(time, value, col=variable, group=interaction(L1,variable)), alpha=.1) 
```

![plot of chunk unnamed-chunk-18](http://farm3.staticflickr.com/2828/9491469452_02e30cee68_o.png) 



## Performance on data outside of observations

Of course it is hardly suprising that all models do reasonably well on the data on which they were trained.  A crux of the problem is the model performance on data outside the observed range.  (Though we might also wish to repeat the above plot on data in the observed range but a different sequence from the observed data).  

First we generate some data from the underlying model coming from below the tipping point: 



```r
Tobs <- 30
y <- numeric(Tobs)
y[1] = .45
for(t in 1:(Tobs-1))
    y[t+1] = z_g() * f(y[t], h=0, p=p)
```


Proceed as before on this data: 


```r
crash_data <- step_ahead_posteriors(y)
```

```
Error: 'Sigma' is not positive definite
```

```r

ggplot(crash_data) + geom_point(aes(time, stock)) + 
  geom_line(aes(time, value, col=variable, group=interaction(L1,variable)), alpha=.1) + 
  facet_wrap(~variable) 
```

```
Error: object 'crash_data' not found
```






## Forecast Distributions

Another way to visualize this is to look directly at the distribution predicted under each model, one step and several steps into the future (say, at a fixed harvest level).  



```r
library(expm)
get_forecasts <- function(i, Tobs, h_i){

df <- data.frame(
  x = x_grid,  
  GP = matrices_gp[[h_i]][i,], 
  True = matrices_true[[h_i]][i,])

df2 <- data.frame(
  x = x_grid,
  GP = (matrices_gp[[h_i]] %^% Tobs)[i,], 
  True = (matrices_true[[h_i]] %^% Tobs)[i,])

T2 <- 2 * Tobs

df4 <- data.frame(
  x = x_grid,
  GP = (matrices_gp[[h_i]] %^% T2)[i,], 
  True = (matrices_true[[h_i]] %^% T2)[i,])

forecasts <- melt(list(T_start = df, T_mid = df2, T_end = df4), id="x")

}

i = 15
```


This takes `i`, an index to `x_grid` value (e.g. for `i`=15 corresponds to a starting postion `x` = 0.3429)



```r
forecasts <- get_forecasts(i = 15, Tobs = 15, h_i = 1)
ggplot(forecasts) +
  geom_line(aes(x, value, group=interaction(variable, L1), col=variable, lty=L1)) +
  facet_wrap(~ variable, scale="free_y", ncol=2) 
```

![plot of chunk unnamed-chunk-22](http://farm3.staticflickr.com/2834/9488674537_253891e90f_o.png) 



We can compare to a better starting stock,


```r
i<-45
```


where `i`=45 corresponds to a starting postion `x` = 1.0776



```r
forecasts <- get_forecasts(i = i, Tobs = 10, h_i = 1)

ggplot(forecasts) +
  geom_line(aes(x, value, group=interaction(variable, L1), col=variable, lty=L1)) +
  facet_wrap(~ variable, scale="free_y", ncol=2) 
```

![plot of chunk unnamed-chunk-24](http://farm6.staticflickr.com/5519/9488674991_bd7fd57cfc_o.png) 




and a middle value showing the population decline:


```r
i<-35
```


where `i`=35 corresponds to a starting postion `x` = 0.8327



```r
forecasts <- get_forecasts(i = i, Tobs = 10, h_i = 1)

ggplot(forecasts) +
  geom_line(aes(x, value, group=interaction(variable, L1), col=variable, lty=L1)) +
  facet_wrap(~ variable, scale="free_y", ncol=2) 
```

![plot of chunk unnamed-chunk-26](http://farm8.staticflickr.com/7290/9488675377_3b9fb86078_o.png) 



