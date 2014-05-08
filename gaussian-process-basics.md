---
layout: page
---



```r
opts_chunk$set(tidy = FALSE, warning = FALSE, message = FALSE, cache = FALSE,
    comment = NA, dev = "CairoPDF", fig.width = 6, fig.height = 4)
library(ggplot2)  # plotting
opts_knit$set(upload.fun = socialR::flickr.url)
theme_set(theme_bw(base_size = 10))
theme_update(panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA))
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
    "#D55E00", "#CC79A7")
```



Basic regression in Gaussian processes
---------------------------------------

Working out coded examples for basic Gaussian process regression using R.  I've just read through the first few chapters of [Rasmussen & Williams (2006)](http://www.GaussianProcess.org/gpml), this implements the examples discussed in Chapter 2.1-2.5.

<!-- Also motivated by [James Keirstead](http://www.jameskeirstead.ca/r/gaussian-process-regression-with-r/), with cleaner code.  -->

Required R libraries (for multivariate normal, also for plotting):


```r
require(MASS)
require(reshape2)
```



Set a seed for repeatable plots


```r
set.seed(12345)
```



Define the points at which we want to compute the function values (x values of the prediction points or test points), and the scale parameter for the covariance function \\( \ell=1 \\)


```r
x_predict <- seq(-5,5,len=50)
l <- 1
```


We will use the squared exponential (also called radial basis or Gaussian, though it is not this that gives Gaussian process it's name; any covariance function would do) as the covariance function, <span> \( \operatorname{cov}(X_i, X_j) = \exp\left(-\frac{(X_i - X_j)^2}{2 \ell ^ 2}\right) \) </span>


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
1 -5.000       X1 -0.4055
2 -4.796       X1 -0.6389
3 -4.592       X1 -0.8300
4 -4.388       X1 -0.9550
5 -4.184       X1 -0.9996
6 -3.980       X1 -0.9581
```



Plot the result


```r
fig2a <- ggplot(dat,aes(x=x,y=value)) +
  geom_rect(xmin=-Inf, xmax=Inf, ymin=-2, ymax=2, fill="grey80") +
  geom_line(aes(group=variable)) +
  scale_y_continuous(lim=c(-3,3), name="output, f(x)") +
  xlab("input, x")
fig2a
```

![plot of chunk fig2a](http://farm8.staticflickr.com/7283/8738379037_e1314a3c34_o.jpg)


### Posterior distribution given the data

Imagine we have some data,


```r
obs <- data.frame(x = c(-4, -3, -1,  0,  2),
                  y = c(-2,  0,  1,  2, -1))
```



In general we aren't interested in drawing from the prior, but want to include information from the data as well.  We want the joint distribution of the observed values and the prior is:

<div>
$$\begin{pmatrix} y_{\textrm{obs}} \\ y_{\textrm{pred}} \end{pmatrix} \sim \mathcal{N}\left( \mathbf{0}, \begin{bmatrix} cov(X_o,X_o) & cov(X_o, X_p) \\ cov(X_p,X_o) & cov(X_p, X_p) \end{bmatrix} \right)$$
</div>

## No observation noise

Assuming the data are known without error and conditioning on the data, and given <span>\(x \sim \mathcal{N}(0, cov(X_o, X_o)) \)</span>, then the conditional probability of observing our data is easily solved by exploiting the nice properties of Gaussians,

<div>

$$x|y \sim \mathcal{N}(E,C)$$

$$E = cov(X_p, X_o) cov(X_o,X_o)^{-1} y$$

$$C= cov(X_p, X_p) - cov(X_p, X_o) cov(X_o,X_o)^{-1} cov(X_o, X_p)$$

</div>

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

![plot of chunk fig2b](http://farm8.staticflickr.com/7288/8738379253_b298f9e11a_o.jpg)






```r
ggplot(data=obs,aes(x=x,y=y)) +
  geom_point()   #OBSERVED DATA
```

![plot of chunk fig2_data](http://farm8.staticflickr.com/7281/8739499454_0b63134241_o.jpg)



Additive noise
--------------

In general the model may have process error, and rather than observe the deterministic map \\(f(x)\\) we only observe \\(y = f(x) + \varepsilon\\).  Let us assume for the moment that \\(\varepsilon\\) are independent, normally distributed random variables with variance \\(\sigma_n^2\\).


```r
sigma.n <- 0.3
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
fig2c
```

![plot of chunk fig2c](http://farm8.staticflickr.com/7284/8739499586_dfe825b0e4_o.jpg)


Note that unlike the previous case, the posterior no longer collapses completely around the neighborhood of the test points.


We can also compute the likelihood (and marginal likelihood over the prior) of the data directly from the inferred multivariate normal distribution, which can allow us to tune the hyperparameters such as the characteristic length scale $\ell$ and the observation noise $\sigma_n$.  The most obvious approach would be to do so by maximum likelihood, giving point estimates of the hyper-parameters, though presumably we could be Bayesian about these as well.

<!-- Continuing... 2012-10-22 -->

The likelihood is given by
<div>
$$\log(p(y | X)) = -\tfrac{1}{2} \mathbf{y}^T (K + \sigma_n^2 \mathbf{I})^{-1} y  - \tfrac{1}{2} \log\left| K + \sigma_n^2 \mathbf{I} \right| - \tfrac{n}{2}\log 2 \pi$$
</div>

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



```r
lpriors <- function(pars){
   d.p <- c(5, 5)
  s2.p <- c(5, 5)

  prior <- unname(
    dgamma(exp(pars[1]), d.p[1], scale = d.p[2]) *
    dgamma(exp(pars[2]), s2.p[1], s2.p[2])
#    dunif(exp(pars[1]), 0, 100) *
#    dunif(exp(pars[2]), 0, 100)
  )

  log(prior)
}

posterior <- function(pars, x, y){

  l <- exp(pars[1])
  sigma.n <- exp(pars[2])

  cov <- function(X, Y) outer(X, Y, SE, l)
  I <- diag(1, length(x))
  K <- cov(x, x)

  loglik <- - 0.5 * t(y) %*% solve(K + sigma.n^2 * I) %*% y -
    log(det(K + sigma.n^2*I)) -
    length(y) * log(2 * pi) / 2

  loglik + lpriors(pars)
}

```



```r
posterior(log(pars), obs$x, obs$y)
```

```
       [,1]
[1,] -21.87
```

```r
posterior(log(o$par), obs$x, obs$y)
```

```
       [,1]
[1,] -40.78
```



```r
require(mcmc)
n <- 1e4
out <- metrop(posterior, log(pars), n, x = obs$x, y = obs$y)
out$accept
```

```
[1] 0.1509
```




```r
postdist <- cbind(index=1:n, as.data.frame(exp(out$batch)))
names(postdist) <- c("index", names(pars))
df <- melt(postdist, id="index")
# TRACES
ggplot(df) + geom_line(aes(index, value)) + facet_wrap(~ variable, scale="free", ncol=1)
```

![plot of chunk gp_mcmc](http://farm8.staticflickr.com/7285/8738380125_14bb7a3296_o.jpg)

```r

ggplot(df) + geom_line(aes(index, log(value))) + facet_wrap(~ variable, scale="free", ncol=1)
```

![plot of chunk gp_mcmc](http://farm8.staticflickr.com/7285/8739500886_70b018f8ab_o.jpg)



```r

d.p <- c(5, 5)
s2.p <- c(5, 5)

s2_prior <- function(x) dgamma(x, s2.p[1], s2.p[2])
d_prior <- function(x) dgamma(x, d.p[1], scale = d.p[2])
prior_fns <- list(l = d_prior, sigma.n = s2_prior)


require(plyr)
prior_curves <- ddply(df, "variable", function(dd){
    grid <- seq(min(dd$value), max(dd$value), length = 100)
    data.frame(value = grid, density = prior_fns[[dd$variable[1]]](grid))
})

# Posteriors (easier to read than histograms)
ggplot(df, aes(value)) +
  stat_density(geom="path", position="identity", alpha=0.7) +
  geom_line(data=prior_curves, aes(x=value, y=density), col="red") +
  facet_wrap(~ variable, scale="free", ncol=2)
```

![plot of chunk gp_posteriors](http://farm8.staticflickr.com/7281/8738381585_dc0bd82493_o.jpg)




