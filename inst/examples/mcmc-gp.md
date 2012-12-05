









```r
knit("~/Documents/code/nonparametric-bayes/inst/examples/beverton_holt_data.Rmd")
```

```


processing file:
~/Documents/code/nonparametric-bayes/inst/examples/beverton_holt_data.Rmd
```

![plot of chunk unnamed-chunk-1](http://carlboettiger.info/assets/figures/2012-12-04-ec2ffbc711-unnamed-chunk-1.png) 

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


Normalize data to the unit (hyper)cube.  








```r
gp <- bgp(X=X, XX=XX, Z=Z, meanfn="constant", bprior="b0", BTE=c(1000,6000,2), m0r1=FALSE, verb=4, corr="exp", trace=TRUE, s2.p=c(5,10), tau2.p=c(5,10), s2.lam="fixed", tau2.lam="fixed")
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
r=1000 d=0.760499; n=39

Sampling @ nn=101 pred locs: [with traces]
r=1000 d=0.204246; mh=1 n=39
r=2000 d=0.0573078; mh=1 n=39
r=3000 d=0.0503173; mh=1 n=39
r=4000 d=0.155435; mh=1 n=39
r=5000 d=0.700194; mh=1 n=39

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







Take a look at the trace


```r
names(gp$trace$XX[[1]])
```

```
[1] "index"  "lambda" "s2"     "tau2"   "beta0"  "nug"    "d"      "b"     
[9] "ldetK" 
```

```r
dim(gp$trace$XX[[1]])
```

```
[1] 2500    9
```

```r
gp$trace$XX[[1]][1,]
```

```
  index lambda    s2   tau2 beta0     nug      d b  ldetK
1     1   1511 26.98 0.8792 3.154 0.02537 0.8052 1 -129.2
```

```r
gp$trace$XX[[1]][2500,]
```

```
     index lambda    s2  tau2 beta0     nug      d b  ldetK
2500  2500  686.5 14.39 3.942 7.617 0.05902 0.6202 1 -98.12
```


* `index` is the sampling point (from `BTE`, we see we sample starting at step 1000 and ending at step 6000, recording every 2 steps, so there are $(T-B)/E = 2500$ index points).  

* `lambda`
* `s2`
* `tau2`
* `beta0`
* `nug`
* `d`

Extract the posterior Gaussian process mean and the $\pm 2$ standard deviations over the predicted grid from the fit:


```r

dat <- data.frame(x    = gp$XX[[1]], 
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

![plot of chunk gp-plot](http://carlboettiger.info/assets/figures/2012-12-04-ec2ffbc711-gp-plot.png) 





is mean ZZ.mean or ZZ.km?  is the var I want ZZ.s2 or ZZ.ks2?)  I want the kreiging values (km, ks2).  The others slightly don't agree with the above.  We do this by comparing to the plotted option (which handles the rescaling(?) from `choose.center` first).



```r
center <- tgp:::tgp.choose.center(gp, "mean")
Z.mean <- center$Z
X_vals <- center$X[, 1]
o <- order(X_vals)
Zb.q1 <- Z.mean + 1.96 * sqrt(c(gp$Zp.ks2, gp$ZZ.ks2))
Zb.q2 <- Z.mean - 1.96 * sqrt(c(gp$Zp.ks2, gp$ZZ.ks2))
V <- c(gp$Zp.ks2, gp$ZZ.ks2)[o]
df <- data.frame(x=X_vals[o], y=Z.mean[o], ymin=Zb.q1[o], ymax=Zb.q2[o])



raw <- data.frame(x    = gp$XX[[1]], 
                 y    = gp$ZZ.mean, 
                 ymin = gp$ZZ.mean - 1.96 * sqrt(gp$ZZ.s2), 
                 ymax = gp$ZZ.mean + 1.96 * sqrt(gp$ZZ.s2))


p1 + geom_ribbon(aes(x,y, ymin=ymin, ymax=ymax), fill="red", alpha=.1, data=df) + geom_line(aes(x,y), col="red", data=df)  + geom_ribbon(aes(x,y, ymin=ymin, ymax=ymax), fill="blue", alpha=.1, data=raw) + geom_line(aes(x,y), col="blue", data=raw) 
```

![plot of chunk unnamed-chunk-6](http://carlboettiger.info/assets/figures/2012-12-04-ec2ffbc711-unnamed-chunk-6.png) 



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
 [1,]    0 0.11125 0.08108 0.06383 0.05258 0.04461 0.03865 0.03403 0.03033
 [2,]    0 0.10462 0.07925 0.06347 0.05279 0.04507 0.03921 0.03462 0.03092
 [3,]    0 0.09774 0.07718 0.06297 0.05293 0.04549 0.03976 0.03521 0.03152
 [4,]    0 0.09075 0.07486 0.06229 0.05295 0.04584 0.04027 0.03579 0.03212
 [5,]    0 0.08376 0.07231 0.06142 0.05284 0.04610 0.04071 0.03632 0.03269
 [6,]    0 0.07691 0.06956 0.06036 0.05258 0.04625 0.04108 0.03680 0.03321
 [7,]    0 0.07031 0.06665 0.05910 0.05215 0.04626 0.04133 0.03719 0.03368
 [8,]    0 0.02677 0.04296 0.04825 0.04876 0.04719 0.04473 0.04196 0.03916
 [9,]    0 0.06229 0.06236 0.05677 0.05092 0.04568 0.04115 0.03727 0.03392
[10,]    0 0.05640 0.05914 0.05505 0.05004 0.04531 0.04109 0.03739 0.03417
        [,10]
 [1,] 0.02730
 [2,] 0.02787
 [3,] 0.02846
 [4,] 0.02905
 [5,] 0.02963
 [6,] 0.03018
 [7,] 0.03068
 [8,] 0.03645
 [9,] 0.03102
[10,] 0.03135
```




```r
policies <- melt(data.frame(stock=x_grid, GP = x_grid[opt_gp$D[,1]], Exact = x_grid[opt_true$D[,1]]), id="stock")
policy_plot <- ggplot(policies, aes(stock, stock - value, color=variable)) +
  geom_point() + xlab("stock size") + ylab("escapement") 
policy_plot + theme_notebook
```

![plot of chunk policy_plot](http://carlboettiger.info/assets/figures/2012-12-04-ec2ffbc711-policy_plot.png) 


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

![plot of chunk simplot](http://carlboettiger.info/assets/figures/2012-12-04-ec2ffbc711-simplot.png) 


Total Profits


```r
sum(sim_gp$profit)
```

```
[1] 10
```

```r
sum(sim_true$profit)
```

```
[1] 29.95
```





