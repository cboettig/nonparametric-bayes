




To provide a minimal working example for comparison: If I were to do this manually (following [Rasmussen & Williams (2006)](http://www.GaussianProcess.org/gpml) Chapter 2), I would do:


Consider we have the observed `x,y` points and `x` values where we desire predicted `y` values:


```r
obs <- data.frame(x = c(-4, -3, -1,  0,  2),
                  y = c(-2,  0,  1,  2, -1))
X <- seq(-5,5,len=50)
```




```r
knit("beverton_holt_data.Rmd")
```

![plot of chunk unnamed-chunk-3](http://carlboettiger.info/assets/figures/2012-11-16-73396d797b-unnamed-chunk-3.png) 




Use a radial basis kernel:


```r
SE <- function(Xi,Xj, l=1) exp(-0.5 * (Xi - Xj) ^ 2 / l ^ 2)
cov <- function(X, Y) outer(X, Y, SE)
```


Calculate mean and variance: 


```r
sigma_n <- 0.3
cov_xx_inv <- solve(cov(obs$x, obs$x) + sigma_n^2 * diag(1, length(obs$x)))
Ef <- cov(X, obs$x) %*% cov_xx_inv %*% obs$y
system.time(Cf <- cov(X, X) - cov(X, obs$x)  %*% cov_xx_inv %*% cov(obs$x, X))
```

```
   user  system elapsed 
  0.004   0.000   0.005 
```


Plot


```r
require(ggplot2)
dat <- data.frame(x=X, y=(Ef), ymin=(Ef-2*sqrt(diag(Cf))), ymax=(Ef+2*sqrt(diag(Cf))))    
ggplot(dat) +
  geom_ribbon(aes(x=x,y=y, ymin=ymin, ymax=ymax), fill="grey80") + # Var
  geom_line(aes(x=x,y=y), size=1) + #MEAN
  geom_point(data=obs,aes(x=x,y=y)) +  #OBSERVED DATA
  scale_y_continuous(name="output, f(x)") + xlab("input, x") 
```

![plot of chunk unnamed-chunk-6](http://carlboettiger.info/assets/figures/2012-11-16-73396d797b-unnamed-chunk-6.png) 




## kernlab

I think I see how I get the equivalent expected values in kernlab:

rename hyperparameters


```r
l = 1
lengthscale = 1/(2*l^2)
var = sigma_n^2
```



```r
library(kernlab)

gp <- gausspr(obs$x, obs$y, kernel="rbfdot", kpar=list(sigma=0.5), fit=FALSE, scaled=FALSE, var=.09)
Ef_k <- predict(gp, X)
```


Manually get the covariance.  First we compute $K(x, x)$


```r
K <- kernelMatrix(kernelf(gp), xmatrix(gp))
```


Then compute $(K(x,x) + \sigma_n^2 \mathbb{I})^{-1}$,


```r
Inv <- solve(K + diag(rep(gp@kcall$var, length = length(xmatrix(gp)))))
```


The covariance is $K(x_*, x_*) - K(x_*, x)(K(x, x) - \sigma_n^2\mathbb{I})^{-1}K(x,x_*)$. 


```r
system.time(Cf_k <- kernelMatrix(kernelf(gp), as.matrix(X)) - kernelMatrix(kernelf(gp), as.matrix(X), xmatrix(gp)) %*% Inv %*% kernelMatrix(kernelf(gp), xmatrix(gp), as.matrix(X)))
```

```
   user  system elapsed 
  0.040   0.000   0.039 
```




Identical to computing the $K(x_*, x)(K(x, x) - \sigma_n^2\mathbb{I})^{-1}$ part first, multiplying the off-diagonal matrix block by the inverse


```r
system.time(Cf_k2 <- kernelMatrix(kernelf(gp), as.matrix(X)) - kernelMult(kernelf(gp), as.matrix(X), xmatrix(gp),  Inv) %*% kernelMatrix(kernelf(gp), xmatrix(gp), as.matrix(X)))
```

```
   user  system elapsed 
  0.032   0.000   0.029 
```



which is also faster.  



compare the two in plots


```r
require(ggplot2)
dat_k <- data.frame(x=X, y=(Ef_k), ymin=(Ef_k-2*sqrt(diag(Cf_k))), ymax=(Ef_k+2*sqrt(diag(Cf_k))))    
ggplot(dat_k) +
  geom_ribbon(aes(x=x,y=y, ymin=ymin, ymax=ymax), fill=rgb(0,0,1,.5)) + # Var
  geom_line(aes(x=x,y=y), size=1, col="blue") + #MEAN
  geom_point(data=obs,aes(x=x,y=y)) +  #OBSERVED DATA
  geom_ribbon(data=dat, aes(x=x,y=y, ymin=ymin, ymax=ymax), fill=rgb(1,0,0,.5)) + # Var
  geom_line(data=dat, aes(x=x,y=y), size=1, col="red") + #MEAN
  scale_y_continuous(name="output, f(x)") + xlab("input, x") 
```

![plot of chunk unnamed-chunk-13](http://carlboettiger.info/assets/figures/2012-11-16-73396d797b-unnamed-chunk-13.png) 


