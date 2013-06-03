# Learning ADMB

Plotting and knitr options, (can generally be ignored)






### Model and parameters


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
set.seed(123)
```


### Sample Data


```r
x <- numeric(Tobs)
x[1] <- Xo
for(t in 1:(Tobs-1))
  x[t+1] = z_g() * f(x[t], h=0, p=p)
```



## Maximum Likelihood "by hand"


```r
n = length(x)
loglik <- function(p){ 
  r = p[1]; k = p[2]; c = p[3]; s = p[4];
  mu = x + r * (1 - x / k) * (x - c) / k;
  0.5 * n * log(2 * pi) + n * log(s) + 0.5 * sum(x - mu)^2/ s^2;
  #-sum(dlnorm(x, mu, s), log=TRUE) # above equation is equivalent to results from built-in function
  }
init <- c(1, 1, 1, 1)
o <- optim(init, loglik, method="L", lower=1e-5, upper=100)
o$par
```

```
[1] 1.000e-05 1.000e+02 1.000e+02 9.993e-04
```


Okay, now lets try admb.  We use R2admb which is just a convenient way to write our data and parameters into an admb file.  



```r
# install_github("R2admb", "bbolker", subdir="R2admb") # dev version
library(R2admb)
```



## ADMB definition

We still need to define the model using ADMB notation in the procedure section.  This is mostly like R or C++, with the exception of special functions like `square` in place of `^2`, `norm2` for the sum of squares, and `elem_prod` istead of `*` for the element-wise product of two arrays. The constant `pi` is given as `M_PI`, as typical of C/C++ libraries.  Where these other functions are defined I'm not sure, but some useful guides to [ADMB vector/matrix operations](http://fish.washington.edu/research/MPAM/resources/ADMB_Minte-Vera.pdf) or an (undefined) list of [keywords](http://www.admb-project.org/developers/contribute-documentation/functions/keywords.txt/view)...

The equivalent model 


```r
model <- 
paste("
PARAMETER_SECTION
  vector mu(1,n) // per capita mort prob
      
PROCEDURE_SECTION
  mu = x + r * elem_prod((1 - x / k), (x - c) / k);
  f = 0.5 * n * log(2 * M_PI) + n * log(s) + 0.5 * norm2(x - mu) / square(s);
")
writeLines(model, "model.tpl")
```







```r
setup_admb("/var/admb")
```

```
[1] "/var/admb"
```

```r

df <- data.frame(x=x)
params <- list(r = 1, k = 10, c = 2, s = .1) ## starting parameters
bounds <- list(r = c(1e-5, 100), k=c(1e-5, 100), c=c(1e-5, 100), s = c(1e-5,100)) ## bounds
dat <- c(list(n = nrow(df)), df)
m1 <- do_admb("model",
              data = dat,
              params = params,
              bounds = bounds,
              run.opts = run.control(checkparam="write",
                                     checkdata="write", clean=FALSE))
m1
```

```
Model file: model_gen 
Negative log-likelihood: -423.8 
Coefficients:
       r        k        c        s 
 0.00001 10.08153  6.38403  0.00001 
```



Which now matches the results from `optim`, now that we avoid the k=0 bound.  

Hans suggests adding an error term in the function definitions rather than in limiting the bounds or log transforming the variables:

> The most common plase where
 this goes wrong is 1/0, log(0), sqrt(0), pow(0,1) etc.
> Your suggestion is OK, but usually I prefer to put
 in log(1e-10+my_expression), sqrt(1e-10+my_expression), pow(1e-10+my_expression,1)

