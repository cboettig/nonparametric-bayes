# Learning ADMB

Plotting and knitr options, (can generally be ignored)






### Model and parameters

Uses the model derived in ` citet("10.1080/10236190412331335373")`, of a Ricker-like growth curve with an allee effect, defined in the pdgControl package,



```r
f <- function(x,h,p)  x * exp(p[1] * (1 - x / p[2]) * (x - p[3]) / p[2] ) 
p <- c(1, 10, 5)
K <- 10  # approx, a li'l' less
allee <- 5 # approx, a li'l' less
```


  

Various parameters defining noise dynamics, grid, and policy costs.  


```r
sigma_g <- 0.05
z_g <- function() rlnorm(1,0, sigma_g)
x_grid <- seq(0, 1.5 * K, length=50)
Xo <-  allee + 0.5 # observations start from
Tobs <- 40
```


### Sample Data


```r
  set.seed(1234)
  #harvest <- sort(rep(seq(0, .5, length=7), 5))
  x <- numeric(Tobs)
  x[1] <- Xo
  nz <- 1
  for(t in 1:(Tobs-1))
    x[t+1] = z_g() * f(x[t], h=0, p=p)
  obs <- data.frame(x = c(rep(0,nz), 
                          pmax(rep(0,Tobs-1), x[1:(Tobs-1)])), 
                    y = c(rep(0,nz), 
                          x[2:Tobs]))
raw_plot <- ggplot(data.frame(time = 1:Tobs, x=x), aes(time,x)) + geom_line() + geom_point()
raw_plot
```

![plot of chunk obs](http://farm8.staticflickr.com/7300/8907300763_0a55663127_o.png) 



## Maximum Likelihood "by hand"


```r
set.seed(12345)
estf <- function(p){ 
    mu <- f(obs$x,0,p)
    -sum(dlnorm(obs$y, log(mu), p[4]), log=TRUE)
}
par <- c(p[1]+abs(rnorm(1,0,1)), 
         p[2]+rnorm(1,0, .1), 
         p[3]+rnorm(1,0, .1), 
         sigma_g + abs(rnorm(1,0, .1)))
o <- optim(par, estf, method="L", lower=c(1e-5,1e-5,1e-5,1e-5))
f_alt <- f
p_alt <- c(as.numeric(o$par[1]), as.numeric(o$par[2]), as.numeric(o$par[3]))
sigma_g_alt <- as.numeric(o$par[4])
est <- list(f = f_alt, p = p_alt, sigma_g = sigma_g_alt, mloglik=o$value)
```






```r
# install_github("R2admb", "bbolker", subdir="R2admb") # dev version
library(R2admb)
```



## ADMB definition


```r
model <- 
paste("
PARAMETER_SECTION
  vector mu(1,n) // per capita mort prob
      
PROCEDURE_SECTION
      
  mu = x + r * elem_prod((1 - x / k), (x - c) / k);
  f -= sum( (0.5 * log(2*M_PI) + log(s) + log(x)  + 0.5 * square(log(x) - mu) / square(s)) ); 
")
writeLines(model, "model.tpl")
```



Or try defining with a loop.  (Not used currently)


```r
model <- 
paste("
PARAMETER_SECTION
  vector mu(1,n) // per capita mort prob
      
PROCEDURE_SECTION
  for(int i = 0; i < n; i++)
  {
    mu[i] = x[i] + r * (1 - x[i] / k) * (x[i] - c) / k;
    f -= (0.5 * log(2*M_PI) + log(s) + log(x[i])  + 0.5 * square(log(x[i]) - mu) / square(s)) );
  }
  
")
writeLines(model, "loop_model.tpl")
```






```r
setup_admb("/var/admb")
```

```
[1] "/var/admb"
```

```r

df <- data.frame(x=x)
params <- list(r = 1, k = 1, c = 1, s = 1) ## starting parameters
bounds <- list(r = c(0, 100), k=c(0, 100), c=c(0, 100), s = c(0,100)) ## bounds
dat <- c(list(n = nrow(df)), df)
m1 <- do_admb("model",
              data = dat,
              params = params,
              bounds = bounds,
              run.opts = run.control(checkparam="write",
                                     checkdata="write"))
```

```
Warning: running command './model_gen > model_gen.out' had status 1
```

```
Warning: std file missing: some problem with fit, but retrieving parameter
estimates anyway
```



```r
unlink(c("model.tpl",
         "model_gen.tpl",
         "model")) ## clean up leftovers

m1
```

```
Model file: model_gen 
Negative log-likelihood: -1.764e+42 
Coefficients:
        r         k         c         s 
2.154e+00 1.169e-09 8.858e-01 3.370e-01 
```




Not so good. Not sure how to debug.  We can confirm ADMB works correctly on example from the Vignette:




## Ben Bolker Model 

from R2admb Vignette.  Works just fine:




```r
model <- 
paste("
PARAMETER_SECTION
  vector prob(1,nobs) // per capita mort prob     
      
PROCEDURE_SECTION

  dvariable fpen=0.0; // penalty variable
  // power-Ricker
  prob = c*pow(elem_prod(TBL/d,exp(1-TBL/d)),g);
  // penalties: constrain 0.001 <= prob <= 0.999
  prob = posfun(prob,0.001,fpen);
  f += 1000*fpen;
  prob = 1-posfun(1-prob,0.001,fpen);
  f += 1000*fpen;
  // binomial negative log-likelihood
  f -= sum( log_comb(nexposed,Kill)+
  elem_prod(Kill,log(prob))+
  elem_prod(nexposed-Kill,log(1-prob)));
")
writeLines(model, "frogs.tpl")
```




```r
setup_admb("/var/admb")
```

```
[1] "/var/admb"
```

```r
df <- data.frame(TBL = rep(c(9,12,21,25,37),each=3), 
                 Kill = c(0,2,1,3,4,5,0,0,0,0,1,0,0,0,0L))
rfs_params <- list(c = 0.45, d = 13, g=1) ## starting parameters
rfs_bounds <- list(c = c(0, 1), d=c(0, 50), g=c(-1, 25)) ## bounds
rfs_dat <- c(list(nobs = nrow(df),
                  nexposed=rep(10, nrow(df))), df)
m1 <- do_admb("frogs",
              data=rfs_dat,
              params=rfs_params,
              bounds=rfs_bounds,
              run.opts=run.control(checkparam="write",
                                   checkdata="write"))
```



```r
unlink(c("frogs.tpl",
         "frogs_gen.tpl",
         "frogs")) ## clean up leftovers

m1
```

```
Model file: frogs_gen 
Negative log-likelihood: 12.89 
Coefficients:
      c       d       g 
 0.4138 13.3508 18.2479 
```




