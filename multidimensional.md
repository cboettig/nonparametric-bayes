---
layout: page
---multi-dimensional dynamics example
========================================================



```r
library(reshape2)
library(ggplot2)
library(data.table)
library(kernlab)
library(tgp)
```








$$\begin{align}
\frac{dM}{dt} = aMC - \frac{(g-h) M}{M+T} + \gamma M T \\
\frac{dC}{dt} = rTC - dC - a M C 
\end{align}$$


```r
mumby <- function(x, h, p= c(a = 0.1, g = .5, gamma = 0.8, r = 1, d = 0.44)){
  M <- x[1] # Macroalgae
  C <- x[2] # Corals
  a <- p[1] # algal growth rate
  g <- p[2] # grazing rate (reduced by harvesting)
  gamma <- p[3] # algal colinization of dead tufts
  r <- p[4] # Coral growth rate
  d <- p[5] # Coral death rate
  dt <- 0.025 # 
  M_t <- M * exp( dt * (a * M * C - (g - h) * M / (M + (1-M-C)) + gamma * M * (1-M-C)) )
  C_t <- C * exp( dt * (r * (1-M-C) * C - d * C - a * M * C) )
  c(M_t, C_t)
}
                  
```



```r
Tmax <- 2500
X <- matrix(NA, nrow=Tmax, ncol=2)
colnames(X) = c("M", "C")
```





```r
X[1,] = c(0.9,0.9)
for(i in 1:(Tmax-1))
  X[i+1,] = mumby(X[i,],0)
dat <- data.frame(time = 1:Tmax, X)
dat <- dat[seq(1,2500,by=25),]
```



```r
ggplot(dat) + geom_line(aes(time, M), col="red") + geom_line(aes(time, C), col="blue")
```

![plot of chunk unnamed-chunk-4](http://carlboettiger.info/assets/figures/2012-12-11-b74ccc4566-unnamed-chunk-4.png) 







```r
inits <- as.matrix(expand.grid(seq(0.05,0.6,length=11), seq(0.05,0.6, length=11)))
df <- lapply(1:round(length(inits[,1])), function(j){
  X[1,] = inits[j,]
  for(i in 1:(Tmax-1))
    X[i+1,] = mumby(X[i,],0.2)
  data.frame(time = 1:Tmax, X)
  })
df2 <- melt(df, id=c("time", "M", "C"))
ggplot(df2) + geom_line(aes(M, C, group=L1)) + 
  geom_point(data = subset(df2, time==1), aes(M, C, group=L1), col="red")
```

![plot of chunk phase-space](http://carlboettiger.info/assets/figures/2012-12-11-b74ccc4566-phase-space.png) 







Noisy data example


```r
sigma_g = 0.1
X[1,] = c(0.9,0.9)
for(i in 1:(Tmax-1))
  X[i+1,] = rlnorm(1,0, sigma_g) * mumby(X[i,],0)
dat <- data.frame(time = 1:Tmax, X)
dat <- dat[seq(1,2500,by=25),]
```



Multi-dimensional GP, `tgp` method, based on time-delay data:



```r
n <- length(dat[,1])
obs <- data.frame(x1=dat[2:(n-1),2], x2 = dat[1:(n-2),2], y=dat[3:n,2])
x_grid <- seq(0,1,length.out=100)
```



```r
require(rgl)
```

```
Loading required package: rgl
```

```r
plot3d(dat$time, dat$M, dat$C)
plot3d(obs$x1, obs$x2, obs$y)
```




```r
gp <- bgp(X=cbind(obs$x1, obs$x2), XX=cbind(x_grid,x_grid), Z=cbind(obs$y), verb=0,
          meanfn="linear", bprior="b0", BTE=c(2000,6000,2), m0r1=FALSE,
          corr="exp", trace=TRUE, beta = c(0,0,0),
          s2.p = c(50,50), d.p = c(10, 1/0.01, 10, 1/0.01), nug.p = c(10, 1/0.01, 10, 1/0.01),
          s2.lam = "fixed", d.lam = "fixed", nug.lam = "fixed", 
          tau2.lam = "fixed", tau2.p = c(50,1))
```



```r
plot(gp, layout="surf")
```

![plot of chunk unnamed-chunk-8](http://carlboettiger.info/assets/figures/2012-12-11-b74ccc4566-unnamed-chunk-8.png) 



Extract the posterior Gaussian process mean and the $\pm 2$ standard deviations over the predicted grid from the fit:


```r
V <- gp$ZZ.ks2
Ef = gp$ZZ.km
tgp_dat <- data.frame(x1   = gp$XX[[1]], 
                      x2 =  gp$XX[[2]],
                  y   = gp$ZZ.km, 
                 ymin = gp$ZZ.km - 1.96 * sqrt(gp$ZZ.ks2), 
                 ymax = gp$ZZ.km + 1.96 * sqrt(gp$ZZ.ks2))
```





```r
ggplot(tgp_dat)  + geom_ribbon(aes(x1,y,ymin=ymin,ymax=ymax), fill="gray80") +
  geom_line(aes(x1,y)) + geom_point(data=obs, aes(x1,y)) 
```

![plot of chunk x1](http://carlboettiger.info/assets/figures/2012-12-11-b74ccc4566-x1.png) 



```r
ggplot(tgp_dat)  + geom_ribbon(aes(x2,y,ymin=ymin,ymax=ymax), fill="gray80") +
  geom_line(aes(x2,y)) + geom_point(data=obs, aes(x2,y)) 
```

![plot of chunk x2](http://carlboettiger.info/assets/figures/2012-12-11-b74ccc4566-x2.png) 




`kernlab` method

