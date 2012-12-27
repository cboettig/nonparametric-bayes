










Fixed priors on hyperparameters, fixed model type.


```r
s2.p <- c(50,50)  # inverse gamma mean is 
tau2.p <- c(20,1)
d.p = c(10, 1/0.01, 10, 1/0.01)
nug.p = c(10, 1/0.01, 10, 1/0.01) # gamma mean
s2_prior <- function(x) dinvgamma(x, s2.p[1], s2.p[2])
tau2_prior <- function(x) dinvgamma(x, tau2.p[1], tau2.p[2])
d_prior <- function(x) dgamma(x, d.p[1], scale = d.p[2]) + dgamma(x, d.p[3], scale = d.p[4])
nug_prior <- function(x) dgamma(x, nug.p[1], scale = nug.p[2]) + dgamma(x, nug.p[3], scale = nug.p[4])
beta0_prior <- function(x, tau) dnorm(x, 0, tau)
beta = c(0)
priors <- list(s2 = s2_prior, tau2 = tau2_prior, beta0 = dnorm, nug = nug_prior, d = d_prior, ldetK = function(x) 0)
```



```r
profit = function(x,h) pmin(x, h)
delta <- 0.01
OptTime = 20  # stationarity with unstable models is tricky thing
reward = 0
xT <- 0
z_g = function() rlnorm(1, 0, sigma_g)
z_m = function() 1+(2*runif(1, 0,  1)-1) * sigma_m
```




```r
f <- RickerAllee
p <- c(1.1, 10, 5) 
K <- 10
allee <- 5
```




```r
sigma_g <- 0.05
sigma_m <- 0.02
x_grid <- seq(0, 1.5 * K, length=101)
h_grid <- x_grid
```


With parameters `1.1, 10, 5`. 



```r
require(snowfall)
sfInit(par=TRUE, cpu=8)
sfExportAll()
sfLibrary(pdgControl)
sfLibrary(nonparametricbayes)
sfLibrary(reshape2)
sfLibrary(ggplot2)
sfLibrary(data.table)
sfLibrary(tgp)
sfLibrary(kernlab)
sfLibrary(MCMCpack)
sfLibrary(plyr)
```




```r
seed <- round(runif(32) * 1e6)
seed[1] <- 1
seed
```

```
 [1]      1 327173 248867 598945 914776 498256 177673 814670 472764 464349
[11] 382117 378880 977121 467465 436839 959060  71648 425478  74179  55641
[21] 635656 425355 687443 104227 812015 252113 922965 396591 635429  18817
[31] 994700 821002
```

```r
yields <- lapply(seed, 
function(seed_i){
  Xo <- K # observations start from
  x0 <- Xo # simulation under policy starts from
  obs <- sim_obs(Xo, z_g, f, p, Tobs=40, seed = seed_i)
  alt <- par_est(obs)
  est <- par_est_allee(obs, f, p)
  gp <- bgp(X=obs$x, XX=x_grid, Z=obs$y, verb=0,
          meanfn="constant", bprior="b0", BTE=c(2000,16000,2),
          m0r1=FALSE, corr="exp", trace=FALSE, 
          beta = beta, s2.p = s2.p, d.p = d.p, nug.p = nug.p, tau2.p = tau2.p,
          s2.lam = "fixed", d.lam = "fixed", nug.lam = "fixed", tau2.lam = "fixed")      
  gp_plot(gp, f, p, est$f, est$p, alt$f, alt$p, x_grid, obs, seed_i)
#  posteriors_plot(gp, priors) # needs trace=TRUE!
  OPT <- optimal_policy(gp, f, est$f, alt$f,
                        p, est$p, alt$p,
                        x_grid, h_grid, sigma_g, 
                        est$sigma_g, alt$sigma_g, 
                        delta, xT, profit, reward, OptTime)
  plot_policies(x_grid, OPT$gp_D, OPT$est_D, OPT$true_D, OPT$alt_D)
  dt <- simulate_opt(OPT, f, p, x_grid, h_grid, x0, z_g, profit)
  sim_plots(dt, seed=seed_i)
  profits_stats(dt)
})
```

![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-11.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-12.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-13.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-14.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-15.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-16.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-17.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-18.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-19.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-110.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-111.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-112.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-113.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-114.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-115.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-116.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-117.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-118.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-119.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-120.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-121.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-122.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-123.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-124.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-125.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-126.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-127.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-128.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-129.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-130.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-131.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-132.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-133.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-134.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-135.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-136.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-137.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-138.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-139.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-140.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-141.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-142.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-143.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-144.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-145.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-146.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-147.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-148.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-149.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-150.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-151.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-152.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-153.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-154.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-155.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-156.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-157.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-158.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-159.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-160.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-161.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-162.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-163.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-164.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-165.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-166.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-167.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-168.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-169.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-170.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-171.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-172.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-173.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-174.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-175.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-176.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-177.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-178.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-179.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-180.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-181.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-182.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-183.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-184.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-185.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-186.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-187.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-188.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-189.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-190.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-191.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-192.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-193.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-194.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-195.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-196.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-197.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-198.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-199.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-1100.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-1101.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-1102.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-1103.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-1104.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-1105.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-1106.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-1107.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-1108.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-1109.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-1110.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-1111.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-1112.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-1113.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-1114.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-1115.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-1116.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-1117.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-1118.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-1119.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-1120.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-1121.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-1122.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-1123.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-1124.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-1125.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-1126.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-1127.png) ![plot of chunk unnamed-chunk-1](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-1128.png) 



```r
yields <- melt(yields, id=c("method", "V1", "sd"))
ggplot(yields) + geom_histogram(aes(V1)) + facet_wrap(~method, scale='free')
```

```
stat_bin: binwidth defaulted to range/30. Use 'binwidth = x' to adjust
this.
```

```
stat_bin: binwidth defaulted to range/30. Use 'binwidth = x' to adjust
this.
```

```
stat_bin: binwidth defaulted to range/30. Use 'binwidth = x' to adjust
this.
```

```
stat_bin: binwidth defaulted to range/30. Use 'binwidth = x' to adjust
this.
```

![plot of chunk unnamed-chunk-2](figure/2012-12-25-19-44-55-1c0ab7ec73-unnamed-chunk-2.png) 

```r
yields
```

```
        method      V1     sd L1
1           GP 11.1060 1.4953  1
2   Parametric  0.3405 0.2945  1
3         True 12.0060 1.7235  1
4   Structural  0.0825 0.1452  1
5           GP 10.3455 1.4586  2
6   Parametric  6.1170 1.6850  2
7         True 12.0060 1.7235  2
8   Structural 11.4330 1.5393  2
9           GP 10.7190 1.4716  3
10  Parametric  2.4645 1.1723  3
11        True 12.0060 1.7235  3
12  Structural  0.4335 0.3692  3
13          GP 11.1060 1.4953  4
14  Parametric  6.1170 1.6850  4
15        True 12.0060 1.7235  4
16  Structural  0.3405 0.2945  4
17          GP 11.1840 1.5199  5
18  Parametric  0.0255 0.1244  5
19        True 12.0060 1.7235  5
20  Structural  6.7196 0.1072  5
21          GP 11.1060 1.4953  6
22  Parametric  8.7844 2.0371  6
23        True 12.0060 1.7235  6
24  Structural  6.0755 0.5659  6
25          GP 11.1060 1.4953  7
26  Parametric 11.8005 1.7528  7
27        True 12.0060 1.7235  7
28  Structural  8.4319 1.6929  7
29          GP  9.3825 1.4459  8
30  Parametric 11.8005 1.7528  8
31        True 12.0060 1.7235  8
32  Structural  8.7855 1.6799  8
33          GP  9.2400 1.4307  9
34  Parametric 10.6365 1.8054  9
35        True 12.0060 1.7235  9
36  Structural  0.3405 0.2945  9
37          GP  9.3105 1.4167 10
38  Parametric 11.9850 1.6559 10
39        True 12.0060 1.7235 10
40  Structural  0.0000 0.0000 10
41          GP 10.3635 1.4759 11
42  Parametric  0.0000 0.0000 11
43        True 12.0060 1.7235 11
44  Structural  0.7815 0.4080 11
45          GP 11.1060 1.4953 12
46  Parametric  9.9210 1.8263 12
47        True 12.0060 1.7235 12
48  Structural  8.0265 1.4667 12
49          GP 10.4895 1.5019 13
50  Parametric  7.9110 1.7819 13
51        True 12.0060 1.7235 13
52  Structural  0.0000 0.0000 13
53          GP 11.7015 1.5626 14
54  Parametric  0.8850 0.6423 14
55        True 12.0060 1.7235 14
56  Structural  0.0000 0.0000 14
57          GP 11.1060 1.4953 15
58  Parametric  1.2330 0.8506 15
59        True 12.0060 1.7235 15
60  Structural  0.3450 0.2919 15
61          GP 10.5720 1.5127 16
62  Parametric 11.8005 1.7528 16
63        True 12.0060 1.7235 16
64  Structural  9.1725 1.7933 16
65          GP 10.8585 1.5051 17
66  Parametric  6.9570 1.7868 17
67        True 12.0060 1.7235 17
68  Structural  6.5330 0.9145 17
69          GP 11.7825 1.5720 18
70  Parametric  5.1855 1.6043 18
71        True 12.0060 1.7235 18
72  Structural 11.0790 1.5088 18
73          GP 12.0900 1.6192 19
74  Parametric 11.8005 1.7528 19
75        True 12.0060 1.7235 19
76  Structural  8.2890 1.4996 19
77          GP  9.7845 1.4153 20
78  Parametric  8.7105 1.7981 20
79        True 12.0060 1.7235 20
80  Structural  0.3405 0.2945 20
81          GP 11.4330 1.5286 21
82  Parametric  7.0455 1.7495 21
83        True 12.0060 1.7235 21
84  Structural  0.0000 0.0000 21
85          GP  9.7845 1.4153 22
86  Parametric  1.9155 0.8442 22
87        True 12.0060 1.7235 22
88  Structural  5.6055 0.1392 22
89          GP 10.7685 1.4533 23
90  Parametric 10.6365 1.8054 23
91        True 12.0060 1.7235 23
92  Structural  2.6460 1.1181 23
93          GP 10.7190 1.4716 24
94  Parametric  5.2155 1.5894 24
95        True 12.0060 1.7235 24
96  Structural  0.3405 0.2945 24
97          GP 11.1060 1.4953 25
98  Parametric  7.0455 1.7495 25
99        True 12.0060 1.7235 25
100 Structural 11.4330 1.5393 25
101         GP  9.9270 1.4467 26
102 Parametric  3.3555 1.3418 26
103       True 12.0060 1.7235 26
104 Structural  0.3405 0.2945 26
105         GP 10.2780 1.4369 27
106 Parametric  9.3900 1.7569 27
107       True 12.0060 1.7235 27
108 Structural  9.9579 2.0057 27
109         GP 11.1060 1.4953 28
110 Parametric  0.3405 0.2945 28
111       True 12.0060 1.7235 28
112 Structural  5.9051 0.2624 28
113         GP 11.1060 1.4953 29
114 Parametric  5.1855 1.6043 29
115       True 12.0060 1.7235 29
116 Structural  9.1766 1.7858 29
117         GP 11.5605 1.5659 30
118 Parametric  8.1135 1.8593 30
119       True 12.0060 1.7235 30
120 Structural  9.9579 2.0057 30
121         GP 11.1075 1.4950 31
122 Parametric  8.2425 1.7275 31
123       True 12.0060 1.7235 31
124 Structural  8.2890 1.4996 31
125         GP  7.5480 1.6134 32
126 Parametric  6.9570 1.7868 32
127       True 12.0060 1.7235 32
128 Structural  7.0838 0.1023 32
```




<p>Myers RA, Barrowman NJ, Hutchings JA and Rosenberg AA (1995).
&ldquo;Population Dynamics of Exploited Fish Stocks at Low Population Levels.&rdquo;
<EM>Science</EM>, <B>269</B>.
ISSN 0036-8075, <a href="http://dx.doi.org/10.1126/science.269.5227.1106">http://dx.doi.org/10.1126/science.269.5227.1106</a>.

