



### Model and parameters



Various parameters defining noise dynamics, grid, and policy costs.  


```r
sigma_g <- 0.05
sigma_m <- 0.0
z_g <- function() rlnorm(1, 0, sigma_g)
z_m <- function() 1
x_grid <- seq(0, 1.5 * 8, length=50)
h_grid <- x_grid
profit <- function(x,h) pmin(x, h)
delta <- 0.01
OptTime <- 50  # stationarity with unstable models is tricky thing
reward <- 0
xT <- 0
Xo <-  5.7 # observations start from
x0 <- 8 # simulation under policy starts from
Tobs <- 40
MaxT = 1000 # timeout for value iteration convergence
```





```r
f <- RickerAllee
```




```r
set.seed(1234)
seed <- round(runif(12) * 1e6)
seed
```

```
 [1] 113703 622299 609275 623379 860915 640311   9496 232551 666084 514251
[11] 693591 544975
```



```r
require(snowfall)
sfInit(parallel=TRUE, cpu=8)
```

```
R Version:  R version 3.0.1 (2013-05-16) 
```

```r
sfLibrary(nonparametricbayes)
```

```
Library nonparametricbayes loaded.
```

```r
sfExportAll()
```


32 data-sets from the same model


```r
yields <- sfSapply(seed, 
function(seed_i){
  set.seed(seed_i)
  
  p <- c(2 + rnorm(1, 0, .1), 8 + rnorm(1, 0, .1), 5 + rnorm(1, 0, .2)) 
  
  initial_data <- function(){
    x <- numeric(Tobs)
    x[1] <- Xo
    nz <- 1
    for(t in 1:(Tobs-1))
      x[t+1] = z_g() * f(x[t], h=0, p=p)
  
    X = c(rep(0,nz), pmax(rep(0,Tobs-1), x[1:(Tobs-1)]))
    Y = c(rep(0,nz), x[2:Tobs])
  
    out = list(X=X,Y=Y)
    if(x[Tobs-1] < x[1]){
      warning("Training population crashed, trying again")
      out <- initial_data()
    } 
    out
  }
  
  attach(initial_data())
  ## @knitr gp-priors
  s2.p <- c(5,5)  
  d.p = c(10, 1/0.1)

  ## @knitr gp
  gp <- gp_mcmc(X, y=Y, n=1e5, s2.p = s2.p, d.p = d.p)
  gp_dat <- gp_predict(gp, x_grid, burnin=1e4, thin=300)
  
  matrices_gp <- gp_transition_matrix(gp_dat$Ef_posterior, gp_dat$Vf_posterior, x_grid, h_grid) 
  opt_gp <- value_iteration(matrices_gp, x_grid, h_grid, MaxT, xT, profit, delta, reward)

  ## @knitr mle-opt
  matrices_true <- f_transition_matrix(f, p, x_grid, h_grid, sigma_g)
  opt_true <- value_iteration(matrices_true, x_grid, h_grid, OptTime=MaxT, xT, profit, delta=delta)

  OPT = data.frame(GP = opt_gp$D, True = opt_true$D)
  policies <- melt(data.frame(stock=x_grid, sapply(OPT, function(x) x_grid[x])), id="stock")
  names(policies) <- c("stock", "method", "value")

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


  
  Profit <- sims_data[, sum(profit), by=c("reps", "method")]
  tmp <- dcast(Profit, reps ~ method)
  tmp <- tmp / tmp[,"True"]
  tmp <- melt(tmp[2:dim(tmp)[2]])
  tmp$value[tmp$variable == "GP"]
  
})
```

```
The following object is masked from initial_data() (position 3):

    X, Y
```

```
The following object is masked from initial_data() (position 3):

    X, Y
The following object is masked from initial_data() (position 4):

    X, Y
```

```
The following object is masked from initial_data() (position 3):

    X, Y
The following object is masked from initial_data() (position 4):

    X, Y
The following object is masked from initial_data() (position 5):

    X, Y
```

```
The following object is masked from initial_data() (position 3):

    X, Y
The following object is masked from initial_data() (position 4):

    X, Y
The following object is masked from initial_data() (position 5):

    X, Y
The following object is masked from initial_data() (position 6):

    X, Y
```

```
The following object is masked from initial_data() (position 3):

    X, Y
The following object is masked from initial_data() (position 4):

    X, Y
The following object is masked from initial_data() (position 5):

    X, Y
The following object is masked from initial_data() (position 6):

    X, Y
The following object is masked from initial_data() (position 7):

    X, Y
```

```
The following object is masked from initial_data() (position 3):

    X, Y
The following object is masked from initial_data() (position 4):

    X, Y
The following object is masked from initial_data() (position 5):

    X, Y
The following object is masked from initial_data() (position 6):

    X, Y
The following object is masked from initial_data() (position 7):

    X, Y
The following object is masked from initial_data() (position 8):

    X, Y
```

```
The following object is masked from initial_data() (position 3):

    X, Y
The following object is masked from initial_data() (position 4):

    X, Y
The following object is masked from initial_data() (position 5):

    X, Y
The following object is masked from initial_data() (position 6):

    X, Y
The following object is masked from initial_data() (position 7):

    X, Y
The following object is masked from initial_data() (position 8):

    X, Y
The following object is masked from initial_data() (position 9):

    X, Y
```

```
The following object is masked from initial_data() (position 3):

    X, Y
The following object is masked from initial_data() (position 4):

    X, Y
The following object is masked from initial_data() (position 5):

    X, Y
The following object is masked from initial_data() (position 6):

    X, Y
The following object is masked from initial_data() (position 7):

    X, Y
The following object is masked from initial_data() (position 8):

    X, Y
The following object is masked from initial_data() (position 9):

    X, Y
The following object is masked from initial_data() (position 10):

    X, Y
```

```
The following object is masked from initial_data() (position 3):

    X, Y
The following object is masked from initial_data() (position 4):

    X, Y
The following object is masked from initial_data() (position 5):

    X, Y
The following object is masked from initial_data() (position 6):

    X, Y
The following object is masked from initial_data() (position 7):

    X, Y
The following object is masked from initial_data() (position 8):

    X, Y
The following object is masked from initial_data() (position 9):

    X, Y
The following object is masked from initial_data() (position 10):

    X, Y
The following object is masked from initial_data() (position 11):

    X, Y
```

```
The following object is masked from initial_data() (position 3):

    X, Y
The following object is masked from initial_data() (position 4):

    X, Y
The following object is masked from initial_data() (position 5):

    X, Y
The following object is masked from initial_data() (position 6):

    X, Y
The following object is masked from initial_data() (position 7):

    X, Y
The following object is masked from initial_data() (position 8):

    X, Y
The following object is masked from initial_data() (position 9):

    X, Y
The following object is masked from initial_data() (position 10):

    X, Y
The following object is masked from initial_data() (position 11):

    X, Y
The following object is masked from initial_data() (position 12):

    X, Y
```

```
The following object is masked from initial_data() (position 3):

    X, Y
The following object is masked from initial_data() (position 4):

    X, Y
The following object is masked from initial_data() (position 5):

    X, Y
The following object is masked from initial_data() (position 6):

    X, Y
The following object is masked from initial_data() (position 7):

    X, Y
The following object is masked from initial_data() (position 8):

    X, Y
The following object is masked from initial_data() (position 9):

    X, Y
The following object is masked from initial_data() (position 10):

    X, Y
The following object is masked from initial_data() (position 11):

    X, Y
The following object is masked from initial_data() (position 12):

    X, Y
The following object is masked from initial_data() (position 13):

    X, Y
```



```r
yields_dat <- melt(yields)
names(yields_dat) <- c("replicate", "simulation", "value")
#class(yields_dat$simulation) = "factor"

col_metadata <- c(replicate = "id of individual replicates using identical settings, relative to the theoretical optimum",
  simulation = "id number for the simulation, which varies the simulated data used",
  value = "total profit achieved by specified run")
unit_metadata <- list(replicate = list(1:dim(yields)[1]), simulation = list(1:dim(yields)[2]), value = "currency (arbitrary)")
```



```r
head(yields_dat)
```

```
  replicate simulation  value
1         1          1 0.3670
2         2          1 0.3368
3         3          1 0.4667
4         4          1 0.3798
5         5          1 0.3936
6         6          1 0.3889
```



```r
ggplot(yields_dat) + geom_density(aes(value)) 
```

![plot of chunk unnamed-chunk-5](http://farm6.staticflickr.com/5455/9387766618_9772fca7e1_o.png) 



Compare over 12 simulation configurations with different parameters and different training data, the GP does consistently well.  Only one case has relatively poor outcomes.  



```r
ggplot(yields_dat) + geom_density(aes(value, color=as.factor(simulation), fill=as.factor(simulation)), alpha=.5)
```

![plot of chunk unnamed-chunk-6](http://farm4.staticflickr.com/3834/9385132341_ea3cffbb8d_o.png) 



Within a replicate number, outcomes not that varied over parameters.  


```r
ggplot(yields_dat) + geom_density(aes(value, group=as.factor(replicate)), alpha=.7)
```

![plot of chunk unnamed-chunk-7](http://farm4.staticflickr.com/3733/9385133429_af9ff61c2a_o.png) 



```r
save(list=ls(), file="sensitivity.rda")
```

