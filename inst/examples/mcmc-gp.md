




```r
library(MCMCpack)
library(tgp)
library(reshape2)
library(ggplot2)
```



I am trying to understand the interface for tgp method, but much of it is not well documented (yes, despite reading through the package R manual, two nice vignettes, and a the PhD Thesis describing it.)

Just to get some consistent notation down: consider the case of a Gaussian process with a Gaussian/radial basis function kernel, parameterized as:

$$K(X, X') = \sigma^2 e^{\frac{(X-X')^2}{d}}$$

And observations from $Z = f(X) + \varepsilon$, for which we seek to approximate $f$ as a (multivariate) Gaussian process, $Z | X ~ N(\mu, C)$. The covariance matrix $C$ is given by our kernel $K$, conditioned on our observations; the mean $\mu$ is given by a constant or linear model, likewise conditioned on our observations.  If the multivariate normal of observed and predicted points is,

$$\begin{pmatrix} y_{\textrm{obs}} \ y_{\textrm{pred}} \end{pmatrix} \sim \mathcal{N}\left( \mathbf{0}, \begin{bmatrix} cov(X_o,X_o) & cov(X_o, X_p) \ cov(X_p,X_o) & cov(X_p, X_p) \end{bmatrix} \right)$$

the conditional probability is

$$x|y \sim \mathcal{N}(E,C)$$ $$E = cov(X_p, X_o) (cov(X_o,X_o) + \varepsilon \mathbb{I}) ^{-1} y$$ $$C= cov(X_p, X_p) - cov(X_p, X_o) (cov(X_o,X_o)+ \varepsilon \mathbb{I} )^{-1} cov(X_o, X_p)$$


Our hyperparameters are $d$, $\sigma^2$, and $\varepsilon$, for which we need priors.  We'll use the inverse gamma,

$$P(x; a, g) = \frac{g^a}{\Gamma(a)} x^{-a - 1}\exp\left(-\frac{g}{x}\right)$$

The distribution arises as the marginal posterior distribution for the unknown variance of a normal distribution if an uninformative prior is used; and as a useful conjugate prior if a less uninformative prior is preferred. However, it is common among Bayesians to consider an alternative parametrization of the normal distribution in terms of the precision, defined as the reciprocal of the variance, which allows the gamma distribution to be used directly as a conjugate prior.

* $X$ is `X`, $Z$ is `Z` 

* $\sigma^2$ is the parameter `s2` in the `tgp` package, coming from inverse gamma prior with parameters $a = a_0$ and $g = g_0$, which are set by `s2.p = c(a_0, g_0)`.  Optionally, $a_0$ and $g_0$ can be taken from an exponential prior $\lambda e^{-\lambda x}$, in which the value of $\lambda$ is given by `s2.lam = c(lambda_a, lambda_g)`.  Setting `s2.lam = "fixed"` "turns off" the use of a hyperprior.  

* $d$ I believe is `d`, which Grammarcy calls the range parameter (other terms include the "lengthscale"). I am not sure what either `d` or `s2` correspond to if the covariance function `corr` is something other than the radial basis function.  I'm also not clear how to set the correlation function to a strictly Gaussian covariance -- `corr` is set to the `exp` family, which leaves the power as a parameter `p_0`:

$$K(X, X') =  e^{\frac{(X-X')^p_0}{d}}$$

(from equation 4 in the Vignette (1)).  I do not see a prior for `p_0` or any mention of it.  Nor can I tell what is the functional form for it's prior (`d.p`) (which is not an inverse gamma?)?  Not clear why it's prior gets defined with four parameters, `c(a1, g1, a2, g2)`, instead of two.   Form the manual:

> Mixture of gamma prior parameter (initial values) for the range parameter(s) c(a1,g1,a2,g2) where g1 and g2 are scale (1/rate) parameters.  

As before, these parameters  `c(a1, g1, a2, g2)` can be drawn from an exponential distribution hyperprior whose $\lambda$ values are given by `d.lam`, or turned off.  

* $\varepsilon$ is `nug`.  It seems to have same shape prior as $d$, specified by `c(a1,g1,a2,g2)`, and optional hyperprior as well.  Additionally it can be fixed to a constant:

> specifying nug.p = 0 fixes the nugget parameter to the starting value in gd[1], i.e., it is excluded from the MCMC

But what is `gd[1]` and where is it defined?

If the mean is not fixed at zero but is instead a linear model, we need hyperparameters for it as well.  The key parameter there is the slope $\beta$.  The functional form for it's prior is specified as `bprior`, which takes values `bflat`, `b0`, `bmzt`, or `bmznot`, but I can find no precise statement of the functional forms of these priors (e.g. I see that `b0` is Gaussian prior and `bflat` an improper or flat prior, but beyond that), nor do I see where their shape parameters are set.  

The manual suggests that we need to specify start values for the hyperparameters for the MCMC (presumably we could have just drawn them from the prior? we burn in anyhow...)  The are given in 

> `start =  c(0.5,0.1,1.0,1.0)` starting values for range `d`, nugget `g`, `s2`, and `tau2`

Actually specifiying start values throws an error though, it appears the documentation is out of date and the start parameter has been removed.  

* I'm still not clear on the role of `tau2`, whether it's unique to the "treed" partioning or whether I've misunderstood $sigma^2$.  Grammacy phrases the model as: 

$$Z | \beta, \sigma^2, K \sim N_n(\mathbf{F} \beta, \sigma^2 \mathbf{K}) $$
$$ \beta | \sigma^2, \tau^2, \mathbf{W}, \beta_0 ~ N_m(\beta_0, \sigma^2 \tau^2 \mathbf{W} )  $$
$$\beta_0 | N_m(\mathbf{\mu}, \mathbf{B})$$

(I guess I haven't parsed the second two equations)


## Example comparison

So, how now would I make a call to `tgp::bgp` that would give me the equivalent output of, say, a direct calculation of the above model (given fixed hyperparameters)?  Say, given some observed data pairs $x,y$ and predictions desired on a grid $X$,


```r
X = seq(-5, 5, length=50)
x = c(-4, -3, -1,  0,  2)
y = c(-2,  0,  1,  2, -1)
```


Manually I could just do:


```r
d = .5; epsilon = .1; sigma = 1; #fixed hyperparamaters
SE <- function(Xi,Xj, d) sigma * exp(-0.5 * (Xi - Xj) ^ 2 / d ^ 2)
cov <- function(X, Y) outer(X, Y, SE, d) 
cov_xx_inv <- solve(cov(x, x) + epsilon * diag(1, length(x)))
Ef <- cov(X, x) %*% cov_xx_inv %*% y
Cf <- cov(X, X) - cov(X, x)  %*% cov_xx_inv %*% cov(x, X)
```


And `Ef` and `Cf` are the estimated (krieging) mean and covariance.

## Priors


```r
s2.p <- c(50,50)
tau2.p <- c(20,1)
d.p = c(10, 1/0.01, 10, 1/0.01)
nug.p = c(10, 1/0.01, 10, 1/0.01)
```



Define as random samples (avoids having to specify range)


```r
s2_rprior <- function(x) rinvgamma(x, s2.p[1], s2.p[2])
tau2_rprior <- function(x) rinvgamma(x, tau2.p[1], tau2.p[2])
d_rprior <- function(x) rgamma(x, d.p[1], scale = d.p[2]) + rgamma(x, d.p[3], scale = d.p[4])
nug_rprior <- function(x) rgamma(x, nug.p[1], scale = nug.p[2]) + rgamma(x, nug.p[3], scale = nug.p[4])
beta0_rprior <- function(x, tau) rnorm(x, 0, tau)
n = 1000
rpriors <- data.frame(s2 = s2_rprior(n), tau2 = tau2_rprior(n), beta0 = beta0_rprior(n, 1), nug = nug_rprior(n), d = d_rprior(n))
rpriors <- melt(rpriors)
ggplot(rpriors) + geom_histogram(aes(value)) + facet_wrap(~variable, scale="free")
```

![plot of chunk unnamed-chunk-5](http://carlboettiger.info/assets/figures/2012-12-10-d0529e905d-unnamed-chunk-5.png) 




Define as curves


```r
s2_prior <- function(x) dinvgamma(x, s2.p[1], s2.p[2])
tau2_prior <- function(x) dinvgamma(x, tau2.p[1], tau2.p[2])
d_prior <- function(x) dgamma(x, d.p[1], scale = d.p[2]) + dgamma(x, d.p[3], scale = d.p[4])
nug_prior <- function(x) dgamma(x, nug.p[1], scale = nug.p[2]) + dgamma(x, nug.p[3], scale = nug.p[4])
beta0_prior <- function(x, tau) dnorm(x, 0, tau)

xx <- seq(.0001, 5, length.out=100)
priors <- data.frame(x = xx, s2 = s2_prior(xx), tau2 = tau2_prior(xx), beta0 = beta0_prior(xx, 1), nug = nug_prior(xx), d = d_prior(xx))
priors <- melt(priors, id="x")
ggplot(priors) + geom_line(aes(x, value)) + facet_wrap(~variable, scale="free")
```

![plot of chunk unnamed-chunk-6](http://carlboettiger.info/assets/figures/2012-12-10-d0529e905d-unnamed-chunk-6.png) 



```r
gp <- bgp(X=x, XX=X, Z=y, verb=0,
          meanfn="constant", bprior="b0", BTE=c(1,2000,1), m0r1=FALSE, 
          corr="exp", trace=TRUE, beta = 0,
          s2.p = s2.p, d.p = d.p, nug.p = nug.p,
          s2.lam = "fixed", d.lam = "fixed", nug.lam = "fixed", 
          tau2.lam = "fixed", tau2.p = tau2.p)
```




Extract the posterior Gaussian process mean and the $\pm 2$ standard deviations over the predicted grid from the fit:


```r
V <- gp$ZZ.ks2
tgp_dat <- data.frame(x   = gp$XX[[1]], 
                  y   = gp$ZZ.km, 
                 ymin = gp$ZZ.km - 1.96 * sqrt(gp$ZZ.ks2), 
                 ymax = gp$ZZ.km + 1.96 * sqrt(gp$ZZ.ks2))
```


and likewise for our manual data:


```r
manual_dat <- data.frame(x=X, y=Ef, ymin=(Ef - 1.96 * sqrt(diag(Cf))), ymax=(Ef + 2 * sqrt(diag(Cf))))
```



Compare the GP posteriors:


```r
ggplot(tgp_dat) +
    geom_ribbon(aes(x, y, ymin = ymin, ymax = ymax), fill="red", alpha = .1) + # Var
    geom_line(aes(x, y), col="red") + # mean
    geom_point(data = data.frame(x = x, y = y), aes(x, y)) + # raw data
    geom_ribbon(dat = manual_dat, aes(x, y, ymin = ymin, ymax = ymax), fill = "blue", alpha = .1) + # Var
    geom_line(dat = manual_dat, aes(x, y), col = "blue")  + #MEAN    
    theme_bw() + theme(plot.background = element_rect(fill = "transparent",colour = NA))
```

![plot of chunk gp-plot](http://carlboettiger.info/assets/figures/2012-12-10-d0529e905d-gp-plot.png) 



Look at the actual parameter estimates from the trace


```r
gp$trace$XX[[1]][1999,]
```

```
     index lambda    s2    tau2  beta0    nug      d b  ldetK
1999  1999  17.56 1.258 0.03925 -1.975 0.2159 0.1163 1 -1.601
```




Posteriors


```r
require(reshape2)
hyperparameters <- c("index", "s2", "tau2", "beta0", "nug", "d", "ldetK")
posteriors <- melt(gp$trace$XX[[1]][,hyperparameters], id="index")
ggplot(posteriors) + geom_histogram(aes(value)) + facet_wrap(~variable, scales="free")
```

![plot of chunk unnamed-chunk-11](http://carlboettiger.info/assets/figures/2012-12-10-d0529e905d-unnamed-chunk-11.png) 






## Compare to priors with `scale` misplaced as `rate` in `d`, `nug`


Define as random samples (avoids having to specify range)


```r
s2_rprior <- function(x) rinvgamma(x, s2.p[1], s2.p[2])
tau2_rprior <- function(x) rinvgamma(x, tau2.p[1], tau2.p[2])
d_rprior <- function(x) rgamma(x, d.p[1], d.p[2]) + rgamma(x, d.p[3], d.p[4])
nug_rprior <- function(x) rgamma(x, nug.p[1], nug.p[2]) + rgamma(x, nug.p[3], nug.p[4])
beta0_rprior <- function(x, tau) rnorm(x, 0, tau)
n = 1000
rpriors <- data.frame(s2 = s2_rprior(n), tau2 = tau2_rprior(n), beta0 = beta0_rprior(n, 1), nug = nug_rprior(n), d = d_rprior(n))
rpriors <- melt(rpriors)
ggplot(rpriors) + geom_histogram(aes(value)) + facet_wrap(~variable, scale="free")
```

![plot of chunk unnamed-chunk-12](http://carlboettiger.info/assets/figures/2012-12-10-d0529e905d-unnamed-chunk-12.png) 










```r
s2_prior <- function(x) dinvgamma(x, s2.p[1], s2.p[2])
tau2_prior <- function(x) dinvgamma(x, tau2.p[1], tau2.p[2])
d_prior <- function(x) dgamma(x, d.p[1], d.p[2]) + dgamma(x, d.p[3], d.p[4])
nug_prior <- function(x) dgamma(x, nug.p[1], nug.p[2]) + dgamma(x, nug.p[3], nug.p[4])
beta0_prior <- function(x, tau) dnorm(x, 0, tau)

xx <- seq(.0001, 2500, length.out=1000)
priors <- data.frame(x = xx, s2 = s2_prior(xx), tau2 = tau2_prior(xx), beta0 = beta0_prior(xx, 1), nug = nug_prior(xx), d = d_prior(xx))
priors <- melt(priors, id="x")
ggplot(priors) + geom_line(aes(x, value)) + facet_wrap(~variable, scale="free")
```

![plot of chunk unnamed-chunk-13](http://carlboettiger.info/assets/figures/2012-12-10-d0529e905d-unnamed-chunk-13.png) 



