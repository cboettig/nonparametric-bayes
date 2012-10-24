




Basic regression in Gaussian processes  
---------------------------------------

Working out coded examples for basic Gaussian process regression using R.  I've just read through the first few chapters of [Rasmussen & Williams (2006)](http://www.GaussianProcess.org/gpml), this implements the examples discussed in Chapter 2.1-2.5.  

<!-- Also motivated by [James Keirstead](http://www.jameskeirstead.ca/r/gaussian-process-regression-with-r/), with cleaner code.  -->

Required R libraries (for multivariate normal, also for plotting):


```r
require(MASS)
require(reshape2)
require(ggplot2)
```



Set a seed for repeatable plots


```r
set.seed(12345)
```



Define the points at which we want to compute the function values (x values of the prediction points or test points), and the scale parameter for the covariance function $\ell=1$


```r
x_predict <- seq(-5,5,len=50)
l <- 1
```


We will use the squared exponential (also called radial basis or Gaussian, though it is not this that gives Gaussian process it's name; any covariance function would do) as the covariance function, $\operatorname{cov}(X_i, X_j) = \exp\left(-\frac{(X_i - X_j)^2}{2 \ell ^ 2}\right)$


```r
SE <- function(Xi,Xj, l) exp(-0.5 * (Xi - Xj) ^ 2 / l ^ 2)
cov <- function(X, Y) outer(X, Y, SE, l)
```



```r
COV <- cov(x_predict, x_predict)
```

Generate a number of functions from the process


```r
values <- mvrnorm(3, rep(0, length=length(x_predict)), COV)
```



Reshape the data into long (tidy) form, listing x value, y value, and sample number

```r
dat <- data.frame(x=x_predict, t(values))
dat <- melt(dat, id="x")
head(dat)
```

```
       x variable   value
1 -5.000       X1 -0.6450
2 -4.796       X1 -0.9227
3 -4.592       X1 -1.1587
4 -4.388       X1 -1.3277
5 -4.184       X1 -1.4139
6 -3.980       X1 -1.4103
```



Plot the result


```r
fig2a <- ggplot(dat,aes(x=x,y=value)) +
  geom_rect(xmin=-Inf, xmax=Inf, ymin=-2, ymax=2, fill="grey80") +
  geom_line(aes(group=variable)) +   theme_bw() +
  scale_y_continuous(lim=c(-3,3), name="output, f(x)") +
  xlab("input, x")
fig2a
```

![plot of chunk unnamed-chunk-9](../../../assets/figures/2012-10-23-cad2271750-unnamed-chunk-9.png) 


### Posterior distribution given the data

Imagine we have some data,


```r
obs <- data.frame(x = c(-4, -3, -1,  0,  2),
                  y = c(-2,  0,  1,  2, -1))
```



In general we aren't interested in drawing from the prior, but want to include information from the data as well.  We want the joint distribution of the observed values and the prior is:

$$\begin{pmatrix} y_{\textrm{obs}} \\ y_{\textrm{pred}} \end{pmatrix} \sim \mathcal{N}\left( \mathbf{0}, \begin{bmatrix} cov(X_o,X_o) & cov(X_o, X_p) \\ cov(X_p,X_o) & cov(X_p, X_p) \end{bmatrix} \right)$$


## No observation noise

Assuming the data are known without error and conditioning on the data, and given $x \sim \mathcal{N}(0, cov(X_o, X_o))$, then the conditional probability of observing our data is easily solved by exploiting the nice properties of Gaussians,

$$x|y \sim \mathcal{N}(E,C)$$
$$E = cov(X_p, X_o) cov(X_o,X_o)^{-1} y$$
$$C= cov(X_p, X_p) - cov(X_p, X_o) cov(X_o,X_o)^{-1} cov(X_o, X_p)$$

(We use `solve(M)` which with no second argument will simply invert the matrix `M`, but should use the Cholsky decomposition instead for better numerical stability)


```r
cov_xx_inv <- solve(cov(obs$x, obs$x))
Ef <- cov(x_predict, obs$x) %*% cov_xx_inv %*% obs$y
Cf <- cov(x_predict, x_predict) - cov(x_predict, obs$x)  %*% cov_xx_inv %*% cov(obs$x, x_predict)
```



Now lets take 3 random samples from the posterior distribution,

```r
values <- mvrnorm(3, Ef, Cf)
```


and plot our solution (mean, 2 standard deviations, and the ranom samples.)


```r
dat <- data.frame(x=x_predict, t(values))
dat <- melt(dat, id="x")

fig2b <- ggplot(dat,aes(x=x,y=value)) +
  geom_ribbon(data=NULL, 
              aes(x=x_predict, y=Ef, ymin=(Ef-2*sqrt(diag(Cf))), ymax=(Ef+2*sqrt(diag(Cf)))),
              fill="grey80") +
  geom_line(aes(color=variable)) + #REPLICATES
  geom_line(data=NULL,aes(x=x_predict,y=Ef), size=1) + #MEAN
  geom_point(data=obs,aes(x=x,y=y)) +  #OBSERVED DATA
  scale_y_continuous(lim=c(-3,3), name="output, f(x)") +
  xlab("input, x")
fig2b
```

![plot of chunk unnamed-chunk-13](../../../assets/figures/2012-10-23-cad2271750-unnamed-chunk-13.png) 


Additive noise
--------------

In general the model may have process error, and rather than observe the deterministic map $f(x)$ we only observe $y = f(x) + \varepsilon$.  Let us assume for the moment that $\varepsilon$ are independent, normally distributed random variables with variance $\sigma_n^2$.  


```r
sigma.n <- 0.8
cov_xx_inv <- solve(cov(obs$x, obs$x) + sigma.n^2 * diag(1, length(obs$x)))
Ef <- cov(x_predict, obs$x) %*% cov_xx_inv %*% obs$y
Cf <- cov(x_predict, x_predict) - cov(x_predict, obs$x)  %*% cov_xx_inv %*% cov(obs$x, x_predict)
```




Now lets take 3 random samples from the posterior distribution,

```r
values <- mvrnorm(3, Ef, Cf)
```


and plot 


```r
dat <- data.frame(x=x_predict, t(values))
dat <- melt(dat, id="x")
Ey <- data.frame(x=x_predict, y=(Ef), ymin=(Ef-2*sqrt(diag(Cf))), ymax=(Ef+2*sqrt(diag(Cf))))

      
fig2c <- ggplot(dat,aes(x=x,y=value)) +
  geom_ribbon(data=Ey, aes(x=x,y=y, ymin=ymin, ymax=ymax), fill="grey80") + # Var
  geom_line(aes(color=variable)) + #REPLICATES
  geom_line(data=Ey,aes(x=x,y=y), size=1) + #MEAN
  geom_point(data=obs,aes(x=x,y=y)) +  #OBSERVED DATA
  scale_y_continuous(lim=c(-3,3), name="output, f(x)") +
  xlab("input, x") 
fig2c + opts(panel.background = theme_rect(fill = "transparent"), 
        plot.background = theme_rect(fill = "transparent"))
```

![plot of chunk unnamed-chunk-16](../../../assets/figures/2012-10-23-cad2271750-unnamed-chunk-16.png) 


Note that unlike the previous case, the posterior no longer collapses completely around the neighborhood of the test points.  


We can also compute the likelihood (and marginal likelihood over the prior) of the data directly from the inferred multivariate normal distribution, which can allow us to tune the hyperparameters such as the characteristic length scale $\ell$ and the observation noise $\sigma_n$.  The most obvious approach would be to do so by maximum likelihood, giving point estimates of the hyper-parameters, though presumably we could be Bayesian about these as well.  

<!-- Continuing... 2012-10-22 -->

The likelihood is given by 
$$\log(p(y | X)) = -\tfrac{1}{2} \mathbf{y}^T (K + \sigma_n^2 \mathbf{I})^{-1} y  - \tfrac{1}{2} \log\left| K + \sigma_n^2 \mathbf{I} \right| - \tfrac{n}{2}\log 2 \pi$$

Which is:


```r

#' observation pairs, obs$x, obs$x, 
#' SE covariance between two observations 
minusloglik <- function(pars){ 

  l <- pars[1]
  sigma.n <- pars[2]
  
  cov <- function(X, Y) outer(X, Y, SE, l)
  I <- diag(1, length(obs$x))
  K <- cov(obs$x, obs$x) 
  
  0.5 * t(obs$y) %*% solve(K + sigma.n^2 * I) %*% obs$y +
    log(det(K + sigma.n^2*I)) +
    length(obs$y) * log(2 * pi) / 2
  }
```




```r
pars <- c(l=1, sigma.n=1)
o <- optim(pars, minusloglik)
o$par
```

```
       l  sigma.n 
0.703891 0.002685 
```


## A comparison to built-in R utilities

Some basic Gaussian process tools are provided in the `kernlab` package,


```r
library(kernlab)
```


through the `gausspr` function,


```r
gp <- gausspr(obs$x, obs$y, kernel="rbfdot", kpar=list(sigma=1/(2*l^2)), fit=FALSE, scaled=FALSE, var=0.8)
```


Where `rbfdot` indicates the Gaussian (Radial basis function) kernel we use above, with parameter `sigma=1/(2l^2)` (note that gausspr defines the kernel parameter differently, see:


```r
gp@kernelf@.Data
```

```
function (x, y = NULL) 
{
    if (!is(x, "vector")) 
        stop("x must be a vector")
    if (!is(y, "vector") && !is.null(y)) 
        stop("y must a vector")
    if (is(x, "vector") && is.null(y)) {
        return(1)
    }
    if (is(x, "vector") && is(y, "vector")) {
        if (!length(x) == length(y)) 
            stop("number of dimension must be the same on both data points")
        return(exp(sigma * (2 * crossprod(x, y) - crossprod(x) - 
            crossprod(y))))
    }
}
<environment: 0x9415ce0>
```

note this takes vectors `x`, `y` and returns a scalar, not the covaraiance matrix).  The function will try to estimate hyperparameters such as `sigma` if they are not given unless `fit=FALSE`.  The data is also scaled to 0 mean, unit variance unless `scaled=FALSE` is set.  

We can then use the estimated process to infer values on a specified grid through `predict`, e.g. for the example above,


```r
y_p <- predict(gp, x_predict)
dat2 <- data.frame(x = x_predict, variable = "Y_p", value = y_p)
```




```r
ggplot(dat2, aes(x, value)) + geom_point() + geom_line(data=Ey,aes(x=x,y=y), size=1) 
```

![plot of chunk unnamed-chunk-23](../../../assets/figures/2012-10-23-cad2271750-unnamed-chunk-23.png) 


### Fitting hyperparameters by maximum likelihood

If we do not specify `fit=false` then `gausspr` will estimate the hyperparamters by maximum likelihood:



```r
gp_fit <- gausspr(obs$x, obs$y, kernel="rbfdot", kpar=list(sigma=1/(2*l^2)), scaled=FALSE, var=0.8)
```

