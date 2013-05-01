






Basic regression in Gaussian processes  
---------------------------------------

Working out coded examples for basic Gaussian process regression using R.  I've just read through the first few chapters of [Rasmussen & Williams (2006)](http://www.GaussianProcess.org/gpml), this implements the examples discussed in Chapter 2.1-2.5.  

<!-- Also motivated by [James Keirstead](http://www.jameskeirstead.ca/r/gaussian-process-regression-with-r/), with cleaner code.  -->



```r
x_predict <- seq(-5,5,len=50)
obs <- data.frame(x = c(-4, -3, -1,  0,  2),
                  y = c(-2,  0,  1,  2, -1))
```





```r
gp <- gp_mcmc(obs$x, obs$y, d.p = c(2,.5), s2.p = c(2,.5))
```




```r
summary_gp_mcmc(gp)
```

![plot of chunk unnamed-chunk-3](http://farm9.staticflickr.com/8268/8699038171_f89f3d5ff4_o.png) ![plot of chunk unnamed-chunk-3](http://farm9.staticflickr.com/8277/8699038281_503ccd606b_o.png) 




```r
gp_dat <- gp_predict(gp, x_predict)
```




```r
tgp_dat <- 
    data.frame(  x = x_predict, 
                 y = gp_dat$E_Ef, 
                 ymin = gp_dat$E_Ef - 2 * sqrt(gp_dat$E_Vf), 
                 ymax = gp_dat$E_Ef + 2 * sqrt(gp_dat$E_Vf) )

ggplot(tgp_dat) + geom_ribbon(aes(x,y,ymin=ymin,ymax=ymax), fill="gray80") +
    geom_line(aes(x, y), lwd=2, alpha=0.8) + 
    geom_point(data=obs, aes(x,y), alpha=0.8) + 
    xlab(expression(X[t])) + ylab(expression(X[t+1])) 
```

![plot of chunk unnamed-chunk-5](http://farm9.staticflickr.com/8538/8699039223_863d2f447a_o.png) 

