Reed (1979) with Gaussian processes 
========================================================







```r
knit("~/Documents/code/nonparametric-bayes/inst/examples/beverton_holt_data.Rmd")
```

```


processing file:
~/Documents/code/nonparametric-bayes/inst/examples/beverton_holt_data.Rmd
```

```
output file:
/home/cboettig/Documents/code/nonparametric-bayes/inst/examples/beverton_holt_data.md
```

```r
true <- data.frame(x = X, y = sapply(X, f, 0, p))
```


Parameters


```r
x_grid <- seq(0, 1.5 * K, length=101)
sigma_g <- 0.01
delta <- 0.01
OptTime = 20
reward = numeric(length(x_grid))
xT = x_grid[4]
reward = profit(x_grid[101],x_grid[101]) + 1/(1-delta)^OptTime # reward for more than 
```






Noise function, profit function

```r
z_g <- function(sigma_g) rlnorm(1, 0, sigma_g) #1+(2*runif(1, 0,  1)-1)*sigma_g #
#profit <- profit_harvest(1,0,0)
profit = function(x,h) pmin(x, h)
```




Simulation 


```r
Tobs <- 100
x <- numeric(Tobs)
x[1] <- x_grid[2]
for(t in 1:(Tobs-1))
  x[t+1] = z_g(sigma_g) * f(x[t], h=0, p=p)
plot(x)
```

![plot of chunk unnamed-chunk-4](http://carlboettiger.info/assets/figures/2012-12-12-80a5ac4d29-unnamed-chunk-4.png) 



Very crude estimate of parameter values


```r
K = mean(x)
r = 0.5 * sigma_g^2 / var(x)

#cor <- acf(x)$acf[,,1]
#m <- which(cor < 0)[1] - 1
#l <- lm(y ~ x, data=data.frame(x = 1:6, y = log(cor[1:6])))
#r <- - as.numeric(l$coefficients[2] / 2)

f_alt <- Ricker
p_alt <- p
```



Predict the function over the target grid


```r
obs <- data.frame(x=x[1:(T-1)],y=x[2:T])
X <- x_grid
```


Gaussian Process inference from this model




```r
gp <- bgp(X=obs$x, XX=X, Z=obs$y, verb=0,
          meanfn="linear", bprior="b0", BTE=c(2000,6000,2), m0r1=FALSE, 
          corr="exp", trace=TRUE, beta = c(0,0),
          s2.p = c(50,50), d.p = c(10, 1/0.01, 10, 1/0.01), nug.p = c(10, 1/0.01, 10, 1/0.01),
          s2.lam = "fixed", d.lam = "fixed", nug.lam = "fixed", 
          tau2.lam = "fixed", tau2.p = c(50,1))
```



Extract the posterior Gaussian process mean and the $\pm 2$ standard deviations over the predicted grid from the fit:


```r
V <- gp$ZZ.ks2
Ef = gp$ZZ.km
tgp_dat <- data.frame(x   = gp$XX[[1]], 
                  y   = gp$ZZ.km, 
                 ymin = gp$ZZ.km - 1.96 * sqrt(gp$ZZ.ks2), 
                 ymax = gp$ZZ.km + 1.96 * sqrt(gp$ZZ.ks2))
```





```r
true <- data.frame(x=X, y=sapply(X,f, 0, p))
ggplot(tgp_dat)  + geom_ribbon(aes(x,y,ymin=ymin,ymax=ymax), fill="gray80") +
  geom_line(aes(x,y)) + geom_point(data=obs, aes(x,y)) +
  geom_line(data=true, aes(x,y), col='red', lty=2)
```

![plot of chunk unnamed-chunk-9](http://carlboettiger.info/assets/figures/2012-12-12-80a5ac4d29-unnamed-chunk-9.png) 




## Stochastic Dynamic programming solution based on the posterior Gaussian process:


```r
h_grid <- x_grid
```


## Understanding the GP Solution


```r
h <- x_grid[1]
X <- numeric(length(x_grid))
X[18] = 1
S <- Ef-h
 F_ <- t(sapply(1:length(S), function(i){
      if(S[i]>0) {
        out <- dnorm(x_grid, S[i], sqrt(V[i]))
#        negatives <- sum(dnorm(seq(-K,0, length=50), S[i], sqrt(V[i])))
#        out[1] <- out[1] + negatives
        out <- out/sum(out)
      } else {
        out <- numeric(length(x_grid))
        out[1] <- 1
        out
      }
    }))
xt1 <- X %*% F_

xt10 <- xt1
for(s in 1:20)
  xt10 <- xt10 %*% F_

qplot(x_grid, xt10[1,]) + geom_point(aes(y=xt1[1,]), col="grey")
```

![plot of chunk unnamed-chunk-11](http://carlboettiger.info/assets/figures/2012-12-12-80a5ac4d29-unnamed-chunk-11.png) 



whereas 



```r
rownorm <- function(M) t(apply(M, 1, function(x) x/sum(x)))

mu <- sapply(x_grid, f, h, p)
F_true <- t(sapply(mu, function(m){
      if(m>0) {
        out <- dlnorm(x_grid/m, 0, sigma_g)
                
      } else {
        out <- numeric(length(x_grid))
        out[1] <- 1
        out
      }
    }))
    F_true <- rownorm(F_true)
F_true <- rownorm(F_true)
yt1 <- X %*% F_true


yt10 <- yt1
for(s in 1:20)
  yt10 <- yt10 %*% F_true
qplot(x_grid, yt10[1,]) + geom_point(aes(y=yt1[1,]), col="grey")
```

![plot of chunk unnamed-chunk-12](http://carlboettiger.info/assets/figures/2012-12-12-80a5ac4d29-unnamed-chunk-12.png) 



```r
transition <- melt(data.frame(x = x_grid, gp = xt1[1,], parametric = yt1[1,]), id="x")
ggplot(transition) + geom_point(aes(x,value, col=variable))
```

![plot of chunk unnamed-chunk-13](http://carlboettiger.info/assets/figures/2012-12-12-80a5ac4d29-unnamed-chunk-13.png) 






```r
matrices_gp <- gp_transition_matrix(Ef, V, x_grid, h_grid)
opt_gp <- find_dp_optim(matrices_gp, x_grid, h_grid, OptTime, xT, profit, delta, reward=1e6)
```







```r
F <- matrices_gp
Tmax <- 20
    n_x <- length(x_grid)
    n_h <- length(h_grid)
    D <- matrix(NA, nrow=n_x, ncol=Tmax)  
    P <-  outer(x_grid, h_grid, profit)   
    V <- P
    for(t in 1:Tmax){
      D[,(Tmax-t+1)] <- apply(V, 1, which.max) 
      v_t <- apply(V, 1, max) # vector of current values 
      V <- sapply(1:n_h, function(j) # updated value matrix
        P[,j] + (1-delta) * F[[j]] %*% v_t)
    }
```






Optimal solution


```r
matrices_true <- f_transition_matrix(f, p, x_grid, h_grid, sigma_g)
opt_true <- find_dp_optim(matrices_true, x_grid, h_grid, OptTime, xT, profit, delta=delta, reward = reward)
```


Optimal parametric solution in which a parameter is estimated with error


```r
matrices_estimated <- f_transition_matrix(f_alt, p_alt, x_grid, h_grid, sigma_g)
opt_estimated <- find_dp_optim(matrices_estimated, x_grid, h_grid, OptTime, xT, profit, delta=delta, reward = reward)
```




```r
policies <- melt(data.frame(stock=x_grid, 
                            GP = x_grid[opt_gp$D[,1]], 
                            Exact = x_grid[opt_true$D[,1]],
                            Approx = x_grid[opt_estimated$D[,1]]),
                  id="stock")

policy_plot <- ggplot(policies, aes(stock, stock - value, color=variable)) +
  geom_point() + xlab("stock size") + ylab("escapement") 
policy_plot
```

![plot of chunk policy_plot](http://carlboettiger.info/assets/figures/2012-12-12-80a5ac4d29-policy_plot.png) 


We can see what happens when we attempt to manage a stock using this:


```r
z_g <- function() rlnorm(1,0, sigma_g)
```



```r
set.seed(1)
sim_gp <- ForwardSimulate(f, p, x_grid, h_grid, K, opt_gp$D, z_g, profit=profit)
set.seed(1)
sim_true <- ForwardSimulate(f, p, x_grid, h_grid, K, opt_true$D, z_g, profit=profit)
set.seed(1)
sim_est <- ForwardSimulate(f, p, x_grid, h_grid, K, opt_estimated$D, z_g, profit=profit)
```




```r
dat <- list(est = sim_est, gp = sim_gp, true = sim_true)
dat <- melt(dat, id=names(dat[[1]]))
dt <- data.table(dat)
setnames(dt, "L1", "method") 
```



```r
ggplot(dt) + geom_line(aes(time,fishstock, color=method))
```

![plot of chunk sim-fish](http://carlboettiger.info/assets/figures/2012-12-12-80a5ac4d29-sim-fish.png) 




```r
ggplot(dt) + geom_line(aes(time,harvest, color=method))
```

![plot of chunk sim-harvest](http://carlboettiger.info/assets/figures/2012-12-12-80a5ac4d29-sim-harvest.png) 


Total Profits


```r
c( gp = sum(sim_gp$profit), true = sum(sim_true$profit), est = sum(sim_est$profit))
```

```
    gp   true    est 
26.400 27.000  8.958 
```

