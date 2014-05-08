---
layout: page
---## One-step ahead predictor plots and forecasting plots

Set our usual plotting options for the notebook


```r
opts_chunk$set(tidy = FALSE, warning = FALSE, message = FALSE, cache = FALSE, 
    comment = NA, fig.width = 7, fig.height = 4)
library(ggplot2)
opts_knit$set(upload.fun = socialR::flickr.url)
theme_set(theme_bw(base_size = 12))
theme_update(panel.background = element_rect(fill = "transparent", colour = NA), 
    plot.background = element_rect(fill = "transparent", colour = NA))
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
    "#D55E00", "#CC79A7")
```



Now we show how each of the fitted models performs on the training data (e.g. plot of the step-ahead predictors). For the GP, we need to predict on the training data first:



```r
gp_f_at_obs <- gp_predict(gp, x, burnin=1e4, thin=300)
```


For the parametric models, prediction is a matter of feeding in each of the observed data points to the fitted parameters, like so 



```r
step_ahead <- function(x, f, p){
  h = 0
  x_predict <- sapply(x, f, h, p)
  n <- length(x_predict) - 1
  y <- c(x[1], x_predict[1:n])
  y
}
```


Which we apply over each of the fitted models, including the GP, organizing the "expected" transition points (given the previous point) into a data frame.  



```r
df <- melt(data.frame(time = 1:length(x), stock = x, 
                GP = gp_f_at_obs$E_Ef,
                True = step_ahead(x,f,p),  
                MLE = step_ahead(x,f,est$p), 
                Parametric.Bayes = step_ahead(x, f, bayes_pars), 
                Ricker = step_ahead(x,alt$f, ricker_bayes_pars), 
                Myers = step_ahead(x, Myer_harvest, myers_bayes_pars)
                 ), id=c("time", "stock"))

ggplot(df) + geom_point(aes(time, stock)) + 
  geom_line(aes(time, value, col=variable)) +
    scale_colour_manual(values=colorkey) 
```

![plot of chunk unnamed-chunk-3](http://farm8.staticflickr.com/7294/8743766846_a5d6b20f15_o.png) 



## Posterior predictive curves

This shows only the mean predictions.  For the Bayesian cases, we can instead loop over the posteriors of the parameters (or samples from the GP posterior) to get the distribution of such curves in each case.  

We will need a vector version (`pmin` in place of `min`) of this function that can operate on the posteriors, others are vectorized already.


```r
ricker_f <- function(x,h,p){
  sapply(x, function(x){ 
    x <- pmax(0, x-h) 
    pmax(0, x * exp(p[2] * (1 - x / p[1] )) )
  })
}
```



Then we proceed as before, now looping over 100 random samples from the posterior for each Bayesian estimate.  We write this as a function for easy reuse.  


```r
require(MASS)

step_ahead_posteriors <- function(x){
gp_f_at_obs <- gp_predict(gp, x, burnin=1e4, thin=300)
df_post <- melt(lapply(sample(100), 
  function(i){
    data.frame(time = 1:length(x), stock = x, 
                GP = mvrnorm(1, gp_f_at_obs$Ef_posterior[,i], gp_f_at_obs$Cf_posterior[[i]]),
                True = step_ahead(x,f,p),  
                MLE = step_ahead(x,f,est$p), 
                Parametric.Bayes = step_ahead(x, allen_f, pardist[i,]), 
                Ricker = step_ahead(x, ricker_f, ricker_pardist[i,]), 
                Myers = step_ahead(x, myers_f, myers_pardist[i,]))
  }), id=c("time", "stock"))

}


df_post <- step_ahead_posteriors(x)
ggplot(df_post) + geom_point(aes(time, stock)) + 
  geom_line(aes(time, value, col=variable, group=interaction(L1,variable)), alpha=.1) + 
  facet_wrap(~variable) + 
  scale_colour_manual(values=colorkey, guide = guide_legend(override.aes = list(alpha = 1))) 
```

![plot of chunk unnamed-chunk-5](http://farm8.staticflickr.com/7290/8742650211_92d087672b_o.png) 



alternately, try the plot without facets


```r
ggplot(df_post) + geom_point(aes(time, stock)) + 
  geom_line(aes(time, value, col=variable, group=interaction(L1,variable)), alpha=.1) + 
  scale_colour_manual(values=colorkey, guide = guide_legend(override.aes = list(alpha = 1))) 
```

![plot of chunk unnamed-chunk-6](http://farm8.staticflickr.com/7284/8743767710_17a0bc6428_o.png) 



## Performance on data outside of observations

Of course it is hardly suprising that all models do reasonably well on the data on which they were trained.  A crux of the problem is the model performance on data outside the observed range.  (Though we might also wish to repeat the above plot on data in the observed range but a different sequence from the observed data).  

First we generate some data from the underlying model coming from below the tipping point: 



```r
Tobs <- 8
y <- numeric(Tobs)
y[1] = 4.5
for(t in 1:(Tobs-1))
    y[t+1] = z_g() * f(y[t], h=0, p=p)
```


Proceed as before on this data: 


```r
crash_data <- step_ahead_posteriors(y)

ggplot(crash_data) + geom_point(aes(time, stock)) + 
  geom_line(aes(time, value, col=variable, group=interaction(L1,variable)), alpha=.1) + 
  facet_wrap(~variable) + 
  scale_colour_manual(values=colorkey, guide = guide_legend(override.aes = list(alpha = 1))) 
```

![plot of chunk unnamed-chunk-8](http://farm8.staticflickr.com/7282/8742650625_efd126ef62_o.png) 






## Forecast Distributions

Another way to visualize this is to look directly at the distribution predicted under each model, one step and several steps into the future (say, at a fixed harvest level).  



```r
library(expm)
get_forecasts <- function(i, Tobs, h_i){

df <- data.frame(
  x = x_grid,  
  GP = matrices_gp[[h_i]][i,], 
  True = matrices_true[[h_i]][i,], 
  MLE = matrices_estimated[[h_i]][i,], 
  Parametric.Bayes = matrices_par_bayes[[h_i]][i,], 
  Ricker = matrices_alt[[h_i]][i,],
  Myers = matrices_myers[[h_i]][i,])

df2 <- data.frame(
  x = x_grid,
  GP = (matrices_gp[[h_i]] %^% Tobs)[i,], 
  True = (matrices_true[[h_i]] %^% Tobs)[i,], 
  MLE = (matrices_estimated[[h_i]] %^% Tobs)[i,],
  Parametric.Bayes = (matrices_par_bayes[[h_i]] %^% Tobs)[i,], 
  Ricker = (matrices_alt[[h_i]] %^% Tobs)[i,], 
  Myers = (matrices_myers[[h_i]] %^% Tobs)[i,])

T2 <- 2 * Tobs

df4 <- data.frame(
  x = x_grid,
  GP = (matrices_gp[[h_i]] %^% T2)[i,], 
  True = (matrices_true[[h_i]] %^% T2)[i,], 
  MLE = (matrices_estimated[[h_i]] %^% T2)[i,],
  Parametric.Bayes = (matrices_par_bayes[[h_i]] %^% T2)[i,], 
  Ricker = (matrices_alt[[h_i]] %^% T2)[i,], 
  Myers = (matrices_myers[[h_i]] %^% T2)[i,])

forecasts <- melt(list(T_start = df, T_mid = df2, T_end = df4), id="x")

}

i = 15
```


This takes `i`, an index to `x_grid` value (e.g. for `i`=15 corresponds to a starting postion `x` = 3.4286)



```r
forecasts <- get_forecasts(i = 15, Tobs = 5, h_i = 1)

ggplot(forecasts) +
  geom_line(aes(x, value, group=interaction(variable, L1), col=variable, lty=L1)) +
  facet_wrap(~ variable, scale="free_y", ncol=2) + 
  scale_colour_manual(values=colorkey) 
```

![plot of chunk unnamed-chunk-10](http://farm8.staticflickr.com/7283/8743768066_ae4cb1a63d_o.png) 



We can compare to a better starting stock,


```r
i<-30
```


where `i`=30 corresponds to a starting postion `x` = 7.102



```r
forecasts <- get_forecasts(i = i, Tobs = 5, h_i = 1)

ggplot(forecasts) +
  geom_line(aes(x, value, group=interaction(variable, L1), col=variable, lty=L1)) +
  facet_wrap(~ variable, scale="free_y", ncol=2) + 
  scale_colour_manual(values=colorkey) 
```

![plot of chunk unnamed-chunk-12](http://farm8.staticflickr.com/7289/8743768174_106af587ce_o.png) 



