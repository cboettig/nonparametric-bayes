









```r
knit("~/Documents/code/nonparametric-bayes/inst/examples/beverton_holt_data.Rmd")
```

```


processing file:
~/Documents/code/nonparametric-bayes/inst/examples/beverton_holt_data.Rmd
```

![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-04-046b509f8d-unnamed-chunk-1.png) 

```
output file:
/home/cboettig/Documents/code/nonparametric-bayes/inst/examples/beverton_holt_data.md
```

```r
true <- data.frame(x = x_grid, y = sapply(x_grid, f, 0, p))
XX <- as.matrix(x_grid)
X <- as.matrix(obs$x)
Z <- as.matrix(obs$y)
```




## Call to `tgp` function `bgp()`

Tweak the hyperparameter prior's parameters (`s2.p`) to influence the GP inferred.  


```r
gp <- bgp(X=X, XX=XX, Z=Z, meanfn="constant", bprior="b0", BTE=c(1000,6000,2), m0r1=TRUE, verb=4, corr="exp", trace=TRUE, s2.p=c(5,10), tau2.p=c(5,10), s2.lam="fixed", tau2.lam="fixed")
```

```

n=39, d=1, nn=101
BTE=(1000,6000,2), R=1, linburn=0
RNG state RK using rk_seed
preds: data krige
T[alpha,beta,nmin,smin,bmax]=[0,0,10,1,1]
mean function: constant
beta prior: b0 hierarchical
s2[a0,g0]=[5,10]
s2 prior fixed
tau2[a0,g0]=[5,10]
tau2 prior fixed
corr prior: isotropic power
nug[a,b][0,1]=[1,1],[1,1]
nug prior fixed
gamlin=[0,0.2,0.7]
d[a,b][0,1]=[1,20],[10,10]
d prior fixed

burn in:
r=1000 d=0.922051; n=39

Sampling @ nn=101 pred locs: [with traces]
r=1000 d=1.36989; mh=1 n=39
r=2000 d=0.688873; mh=1 n=39
r=3000 d=0.746951; mh=1 n=39
r=4000 d=0.852561; mh=1 n=39
r=5000 d=1.0243; mh=1 n=39

Gathering traces
  XX 63% done     XX 63% done     XX 63% done     XX 63% done     XX 63% done     XX 63% done     XX 63% done     XX 63% done     XX 63% done     XX 63% done     XX 63% done     XX 63% done     XX 64% done     XX 64% done     XX 64% done     XX 64% done     XX 64% done     XX 64% done     XX 64% done     XX 64% done     XX 64% done     XX 64% done     XX 64% done     XX 64% done     XX 64% done     XX 64% done     XX 64% done     XX 65% done     XX 65% done     XX 65% done     XX 65% done     XX 65% done     XX 65% done     XX 65% done     XX 65% done     XX 65% done     XX 65% done     XX 65% done     XX 65% done     XX 65% done     XX 66% done     XX 66% done     XX 66% done     XX 66% done     XX 66% done     XX 66% done     XX 66% done     XX 66% done     XX 66% done     XX 66% done     XX 67% done     XX 67% done     XX 67% done     XX 67% done     XX 67% done     XX 67% done     XX 67% done     XX 67% done     XX 68% done     XX 68% done     XX 68% done     XX 68% done     XX 68% done     XX 68% done     XX 68% done     XX 69% done     XX 69% done     XX 69% done     XX 69% done     XX 69% done     XX 69% done     XX 70% done     XX 70% done     XX 70% done     XX 70% done     XX 71% done     XX 71% done     XX 71% done     XX 71% done     XX 72% done     XX 72% done     XX 72% done     XX 73% done     XX 73% done     XX 73% done     XX 74% done     XX 74% done     XX 75% done     XX 75% done     XX 76% done     XX 77% done     XX 77% done     XX 78% done     XX 79% done     XX 80% done     XX 81% done     XX 83% done     XX 85% done     XX 88% done     XX 92% done     XX 100% done  
  hier-params done
  posts done
  Zp done
  Zp.km done
  Zp.ks2 done
  ZZ done
  ZZ.km done
  ZZ.ks2 done

```


## parsing the function call

* `X` is the observed X values (or matrix of appropriate dimension)
* `XX` the desired predicted grid. (will also predict on X)
* `Z` is observed Y values (or matrix of appropriate dimension)
* `m0r1` means scale data to mean zero and unit range
* `bprior` is the prior distribution, defaulting to `bfalt`. `b0` is a Gaussian.  This prior is for the Gaussian process and is not to be confused with the hyperpriors describing the distribution of various parameters.
* `BTE` Burn-in B steps, Terminate after T steps, sample Every E steps. 
* `verb` Very verbose output
* Other options (for adaptive/query learning, etc, igored for now.  Focus on priors).  
* The priors for the hyperparameters are given by `s2.p`, `tau.p` etc (both are Inverse Gammas).  Their parameters can vary hiearchically in general, but we hold them fixed, `s2.lam = 'fixed'`, etc.)

## Hyperparameters

Hyperparameter distributions returned in the trace:


```r
names(gp$trace$XX[[1]][2500,])
```

```
[1] "index"  "lambda" "s2"     "tau2"   "beta0"  "nug"    "d"      "b"     
[9] "ldetK" 
```


* `s2` is $\sigma^2$ in Gramacy and $\tau$ in Rassmussen and Williams, is the scale of the noise. Comes from an inverse gamma prior.  Set by inverse-gamma parameters `a0` and `g0`.  

* `tau2` Also comes from an inverse gamma prior. It is also set by inverse-gamma parameters `a0` and `g0`.  

$$Z | \beta, \sigma^2, K \sim N_n(\mathbf{F} \beta, \sigma^2 \mathbf{K}) $$
$$ \beta | \sigma^2, \tau^2, \mathbf{W}, \beta_0 ~ N_m(\beta_0, \sigma^2 \tau^2 \mathbf{W} )  $$
$$\beta_0 | N_m(\mathbf{\mu}, \mathbf{B})$$

The parameters in this run are constant across the trace points (because hyperparameters are fixed).  


```r
gp$trace$XX[[1]][2500,]
```

```
     index lambda    s2  tau2  beta0      nug     d b  ldetK
2500  2500  89.83 2.064 3.616 0.6892 0.003601 1.024 1 -199.9
```

```r
gp$trace$XX[[3]][2500,]
```

```
     index lambda    s2  tau2  beta0      nug     d b  ldetK
2500  2500  89.83 2.064 3.616 0.6892 0.003601 1.024 1 -199.9
```




```r
pars <- melt(gp$trace$XX[[1]], id = "index")
ggplot(pars, aes(index, value)) + geom_line() + facet_wrap(~variable, scales="free_y")
```

![plot of chunk unnamed-chunk-5](http://carlboettiger.info/assets/figures/2012-12-04-046b509f8d-unnamed-chunk-5.png) 


It is possible to calculate summary statistics and check out the distribution of these fellows.  


```r
mean(gp$trace$XX[[1]][,"nug"])
```

```
[1] 0.009755
```

```r
ggplot(pars, aes(value)) + geom_histogram() + facet_wrap(~variable, scales="free")
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

![plot of chunk unnamed-chunk-6](http://carlboettiger.info/assets/figures/2012-12-04-046b509f8d-unnamed-chunk-6.png) 


* We keep the heirachical priors fixed in the `tau.lm`, `s2.lm`, etc. `nug` fixed by default.  Note this results in the the priors from which these are sampled are fixed: e.g. `gp$trace$hier$s2.a0`, `gp$trace$hier$s2.g0`

* `nug.p` This is the "nugget," the measurement error.  Comes from a mixture of gamma prior parameter (initial values) for the nugget parameter `c(a1,g1,a2,g2)` where `g1` and `g2` are scale (`1/rate`) parameters. _Default reduces to simple exponential prior._ Specifying `nug.p = 0` fixes the nugget parameter to the "starting" value in `gd[1]`, i.e., it is excluded from the MCMC


* `index` is the sampling point (from `BTE`, we see we sample starting at step 1000 and ending at step 6000, recording every 2 steps, so there are $(T-B)/E = 2500$ index points).  
* `lambda` possibly a mixing parameter in the importance sampler temperature??
* `beta0` No mention in the manual
* `d` The parameter for priors using the hierachical exponential distribution for the parameters `a1`, `a2`, `g1`, `g2`.  
* `b` ??

* Where are the hyperparameters for the correlation function (e.g. the length-scale for the Gaussian??)  


## Extracting the estimated process 

(The **MAP**, maximum a posteriori expected mean and variance).

Extract the posterior Gaussian process mean and the $\pm 2$ standard deviations over the predicted grid from the fit:


```r
V <- gp$ZZ.ks2
dat <- data.frame(x   = gp$XX[[1]], 
                 y    = gp$ZZ.km, 
                 ymin = gp$ZZ.km - 1.96 * sqrt(gp$ZZ.ks2), 
                 ymax = gp$ZZ.km + 1.96 * sqrt(gp$ZZ.ks2))
```


Plot the posterior Gaussian Process:


```r
  p1 <- ggplot(dat) +
    geom_ribbon(aes(x=x,y=y, ymin=ymin, ymax=ymax), fill="grey80") + # Var
    geom_line(aes(x=x,y=y), size=1) + #MEAN
    geom_point(data=gp$obs,aes(x=x,y=y)) + 
    labs(title=paste("llik =", prettyNum(gp$llik)))
  if(!is.null(true))
    p1 <- p1 + geom_line(data = true, aes(x, y), col = "red", lty = 2) 
p1 + theme_notebook
```

![plot of chunk gp-plot](http://carlboettiger.info/assets/figures/2012-12-04-046b509f8d-gp-plot.png) 




Find the actual optimum policy from the true parametric model


```r
x_grid <- dat$x
h_grid <- x_grid
matrices_true <- f_transition_matrix(f, p, x_grid, h_grid)
opt_true <- find_dp_optim(matrices_true, x_grid, h_grid, 20, 0, profit, delta=.01)
```



Estimate a policy from the Gaussian Process:


```r
rownorm <- function(M) t(apply(M, 1, function(x) x/sum(x)))
matrices_gp <- lapply(h_grid, function(h){
    S <- dat$y - h
    F_ <- t(sapply(1:length(S), function(i){
      if(S[i]>0) {
        out <- dlnorm(x_grid/S[i], 0, V[i])
      } else {
        out <- numeric(length(x_grid))
        out[1] <- 1
        out
      }
    }))
    F_ <- rownorm(F_)
  })
opt_gp <- find_dp_optim(matrices_gp, x_grid, h_grid, 20, 0, profit, delta=.01)
```




```r
matrices_true[[1]][1:10,1:10]
```

```
      [,1]       [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
 [1,]    1  0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00
 [2,]    0  5.361e-02 9.464e-01 1.078e-09 2.276e-20 3.993e-31 1.614e-41
 [3,]    0  9.446e-26 7.322e-04 9.915e-01 7.731e-03 5.996e-07 6.519e-12
 [4,]    0  8.377e-48 6.260e-14 8.680e-04 6.372e-01 3.564e-01 5.502e-03
 [5,]    0  1.484e-67 3.055e-25 3.665e-10 8.567e-04 2.498e-01 6.400e-01
 [6,]    0  4.313e-85 2.801e-36 2.128e-17 2.477e-08 8.930e-04 1.172e-01
 [7,]    0 8.411e-101 1.021e-46 9.403e-25 1.687e-13 3.028e-07 9.679e-04
 [8,]    0 5.121e-115 1.654e-56 5.908e-32 7.277e-19 3.471e-11 1.618e-06
 [9,]    0 5.398e-128 1.115e-65 6.708e-39 3.138e-24 2.515e-15 1.175e-09
[10,]    0 6.254e-140 2.784e-74 1.492e-45 1.688e-29 1.601e-19 5.633e-13
           [,8]      [,9]     [,10]
 [1,] 0.000e+00 0.000e+00 0.000e+00
 [2,] 1.969e-51 7.453e-61 8.297e-70
 [3,] 3.099e-17 1.112e-22 3.983e-28
 [4,] 1.210e-05 8.834e-09 3.442e-12
 [5,] 1.060e-01 3.272e-03 3.465e-05
 [6,] 5.412e-01 2.984e-01 4.018e-02
 [7,] 6.648e-02 3.799e-01 4.023e-01
 [8,] 1.072e-03 4.360e-02 2.608e-01
 [9,] 5.465e-06 1.203e-03 3.190e-02
[10,] 1.444e-08 1.394e-05 1.362e-03
```

```r
matrices_gp[[1]][1:10,1:10]
```

```
      [,1]    [,2]    [,3]    [,4]    [,5]    [,6]    [,7]    [,8]    [,9]
 [1,]    0 0.21390 0.12563 0.08674 0.06495 0.05111 0.04161 0.03474 0.02955
 [2,]    0 0.19036 0.11927 0.08504 0.06499 0.05189 0.04271 0.03596 0.03080
 [3,]    0 0.16918 0.11295 0.08312 0.06479 0.05246 0.04363 0.03703 0.03193
 [4,]    0 0.14976 0.10655 0.08091 0.06433 0.05280 0.04437 0.03797 0.03296
 [5,]    0 0.13188 0.10001 0.07838 0.06357 0.05291 0.04493 0.03876 0.03386
 [6,]    0 0.11551 0.09339 0.07555 0.06251 0.05277 0.04529 0.03939 0.03464
 [7,]    0 0.10066 0.08675 0.07245 0.06116 0.05237 0.04542 0.03984 0.03528
 [8,]    0 0.08735 0.08020 0.06913 0.05954 0.05171 0.04534 0.04010 0.03575
 [9,]    0 0.07556 0.07382 0.06566 0.05769 0.05082 0.04503 0.04017 0.03606
[10,]    0 0.06525 0.06772 0.06211 0.05565 0.04971 0.04453 0.04006 0.03621
        [,10]
 [1,] 0.02553
 [2,] 0.02676
 [3,] 0.02789
 [4,] 0.02894
 [5,] 0.02990
 [6,] 0.03075
 [7,] 0.03149
 [8,] 0.03210
 [9,] 0.03257
[10,] 0.03289
```




```r
policies <- melt(data.frame(stock=x_grid, GP = x_grid[opt_gp$D[,1]], Exact = x_grid[opt_true$D[,1]]), id="stock")
policy_plot <- ggplot(policies, aes(stock, stock - value, color=variable)) +
  geom_point() + xlab("stock size") + ylab("escapement") 
policy_plot + theme_notebook
```

![plot of chunk policy_plot](http://carlboettiger.info/assets/figures/2012-12-04-046b509f8d-policy_plot.png) 


We can see what happens when we attempt to manage a stock using this:


```r
z_g <- function() rlnorm(1,0, sigma_g)
```



```r
set.seed(1)
sim_gp <- ForwardSimulate(f, p, x_grid, h_grid, K, opt_gp$D, z_g, profit=profit)
set.seed(1)
sim_true <- ForwardSimulate(f, p, x_grid, h_grid, K, opt_true$D, z_g, profit=profit)
```



```r
df <- data.frame(time = sim_gp$time, stock_gp = sim_gp$fishstock, stock_true = sim_true$fishstock, harvest_gp = sim_gp$harvest, havest_true = sim_true$harvest)
df <- melt(df, id="time")
simplot <- ggplot(df) + geom_line(aes(time,value, color=variable))
simplot + theme_notebook
```

![plot of chunk simplot](http://carlboettiger.info/assets/figures/2012-12-04-046b509f8d-simplot.png) 


Total Profits


```r
sum(sim_gp$profit)
```

```
[1] 16.35
```

```r
sum(sim_true$profit)
```

```
[1] 29.95
```





