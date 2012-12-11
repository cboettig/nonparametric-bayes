GP Example using the May (1979) bistable model
==============================================


Use my custom libraries and the `ggplot` plotting tools


```r
library(pdgControl)
library(nonparametricbayes)
library(reshape2)
library(ggplot2)
library(data.table)
library(kernlab)
library(tgp)
```






State-space model


```r
f <- May
p <- c(r = .75, k = 10, a=1.3, H=1, Q = 3)
K <- 8 # approx
```


Parameters


```r
x_grid = seq(0, 1.5 * K, length=101)
T <- 40
sigma_g <- 0.3
x <- numeric(T)
x[1] <- 1
```



Model dynamics look like this:


```r
birth <- function(x) p["r"] * (1-  x / p["k"])
death <- function(x) p["a"] * x ^ (p["Q"] - 1) / (x ^ p["Q"] + p["H"])
df <- data.frame(x = x_grid, b = sapply(x_grid, birth), d = sapply(x_grid, death))
ggplot(df) + geom_line(aes(x, b), col = "blue") + geom_line(aes(x,d), col = "red")
```

![plot of chunk unnamed-chunk-3](http://carlboettiger.info/assets/figures/2012-12-11-9b2d792442-unnamed-chunk-3.png) 




Noise function, profit function

```r
z_g <- function(sigma_g) rlnorm(1, 0, sigma_g) #1+(2*runif(1, 0,  1)-1)*sigma_g #
profit <- profit_harvest(1,0,0)
```




Simulation 


```r
x[1] = 2.5
for(t in 1:(T-1))
  x[t+1] = z_g(sigma_g) * f(x[t], h=0, p=p)
plot(x)
```

![plot of chunk unnamed-chunk-5](http://carlboettiger.info/assets/figures/2012-12-11-9b2d792442-unnamed-chunk-5.png) 


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

![plot of chunk unnamed-chunk-9](http://carlboettiger.info/assets/figures/2012-12-11-9b2d792442-unnamed-chunk-9.png) 




## Stochastic Dynamic programming solution based on the posterior Gaussian process:


```r
h_grid <- x_grid
```


GP Solution


```r
matrices_gp <- gp_transition_matrix(Ef, sqrt(V), x_grid, h_grid)
opt_gp <- find_dp_optim(matrices_gp, x_grid, h_grid, 20, 0, profit, delta=.01)
```


Optimal solution


```r
matrices_true <- f_transition_matrix(f, p, x_grid, h_grid, sigma_g)
opt_true <- find_dp_optim(matrices_true, x_grid, h_grid, 20, 0, profit, delta=.01)
```


Optimal parametric solution in which a parameter is estimated with error


```r
matrices_estimated <- f_transition_matrix(f, p, x_grid, h_grid, sigma_g=0.01)
opt_estimated <- find_dp_optim(matrices_estimated, x_grid, h_grid, 20, 0, profit, delta=.01)
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

![plot of chunk policy_plot](http://carlboettiger.info/assets/figures/2012-12-11-9b2d792442-policy_plot.png) 


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

![plot of chunk simplot](http://carlboettiger.info/assets/figures/2012-12-11-9b2d792442-simplot1.png) 

```r
ggplot(dt) + geom_line(aes(time,harvest, color=method))
```

![plot of chunk simplot](http://carlboettiger.info/assets/figures/2012-12-11-9b2d792442-simplot2.png) 


Total Profits


```r
c( gp = sum(sim_gp$profit), true = sum(sim_true$profit), est = sum(sim_est$profit))
```

```
   gp  true   est 
 8.00 25.44 25.14 
```


