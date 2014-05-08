---
layout: page
---


```r
f_transition_matrix <- function(f, p, x_grid, h_grid, sigma_g, pardist) {
    lapply(h_grid, par_F, f, p, x_grid, sigma_g, pardist)
}


par_F <- function(h, f, p, x_grid, sigma_g, pardist, n_mc = 100) {
    
    # Set up monte carlo sampling
    d <- dim(pardist)
    indices <- round(runif(n_mc, 1, d[1]))
    
    # compute mean with monte carlo samples
    mu <- sapply(indices, function(i) {
        # Cols are parameter samples
        p <- unname(pardist[i, c(4, 1, 5)])  # Only, we need to add the measurement and process noise
        unname(sapply(x_grid, f, h, p))  # rows are possible values of x_t
    })
    mu <- data.frame(t(mu))  # Let's have cols as possible x_t, rows as parameters, and type as data.frame
    
    F_true <- sapply(mu, function(cur_state) {
        # For each x_t
        
        bypar <- sapply(cur_state, function(m) {
            # For each parameter value
            if (snap_to_grid(m, x_grid) < x_grid[2]) {
                #
                out <- numeric(length(x_grid))
                out[1] <- 1
                out
            } else {
                out <- dlnorm(x_grid/m, 0, sigma_g)  # should be using the estimate process noise!
            }
        })
        ave_over_pars <- apply(bypar, 1, sum)  # collapse by weighted average over possible parameters
        ave_over_pars/sum(ave_over_pars)
    })
    
    F_true <- t(F_true)
    
    
}


snap_to_grid <- function(x, grid) sapply(x, function(x) grid[which.min(abs(grid - 
    x))])


# internal helper function
rownorm <- function(M) t(apply(M, 1, function(x) {
    if (sum(x) > 0) {
        x/sum(x)
    } else {
        out <- numeric(length(x))
        out[1] <- 1
        out
    }
}))
```





```r
knit("allen.Rmd")
```

```
## 
## 
## processing file: allen.Rmd
```

```
##   |                                                                         |                                                                 |   0%  |                                                                         |>>                                                               |   3%
##   ordinary text without R code
## 
##   |                                                                         |>>>>                                                             |   7%
## label: libraries (with options) 
## List of 2
##  $ include: logi FALSE
##  $ cache  : logi FALSE
## 
##   |                                                                         |>>>>>>>                                                          |  10%
##   ordinary text without R code
## 
##   |                                                                         |>>>>>>>>>                                                        |  14%
## label: plotting-options (with options) 
## List of 3
##  $ echo   : logi FALSE
##  $ include: logi FALSE
##  $ cache  : logi TRUE
## 
##   |                                                                         |>>>>>>>>>>>                                                      |  17%
##   ordinary text without R code
## 
##   |                                                                         |>>>>>>>>>>>>>                                                    |  21%
## label: unnamed-chunk-3
##   |                                                                         |>>>>>>>>>>>>>>>>                                                 |  24%
##   ordinary text without R code
## 
##   |                                                                         |>>>>>>>>>>>>>>>>>>                                               |  28%
## label: sdp-pars (with options) 
## List of 1
##  $ dependson: chr "stateeq"
```

```
## Warning: code chunks must not depend on the uncached chunk "sdp-pars"
```

```
##   |                                                                         |>>>>>>>>>>>>>>>>>>>>                                             |  31%
##   ordinary text without R code
## 
##   |                                                                         |>>>>>>>>>>>>>>>>>>>>>>                                           |  34%
## label: obs (with options) 
## List of 1
##  $ dependson: chr "sdp-pars"
```

```
## Warning: code chunks must not depend on the uncached chunk "obs"
```

```
##   |                                                                         |>>>>>>>>>>>>>>>>>>>>>>>>>                                        |  38%
##   ordinary text without R code
## 
##   |                                                                         |>>>>>>>>>>>>>>>>>>>>>>>>>>>                                      |  41%
## label: mle (with options) 
## List of 1
##  $ dependson: chr "obs"
## 
##   |                                                                         |>>>>>>>>>>>>>>>>>>>>>>>>>>>>>                                    |  45%
##   ordinary text without R code
## 
##   |                                                                         |>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>                                  |  48%
## label: unnamed-chunk-4
##   |                                                                         |>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>                               |  52%
##   ordinary text without R code
## 
##   |                                                                         |>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>                             |  55%
## label: unnamed-chunk-5
##   |                                                                         |>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>                           |  59%
##   ordinary text without R code
## 
##   |                                                                         |>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>                         |  62%
## label: unnamed-chunk-6
##   |                                                                         |>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>                      |  66%
##   ordinary text without R code
## 
##   |                                                                         |>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>                    |  69%
## label: unnamed-chunk-7
```

```
## Warning: failed to tidy R code in chunk <unnamed-chunk-7> reason: Error in
## base::parse(text = text) : 1:65: unexpected SPECIAL 1: jags.params = c (
## "K" , "logr0" , "logtheta" , "iR" , "iQ" ) ; %InLiNe_IdEnTiFiEr% ^
```

```
##   |                                                                         |>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>                  |  72%
##   ordinary text without R code
## 
##   |                                                                         |>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>                |  76%
## label: unnamed-chunk-8
```

```
## Loading required package: Rflickr
```

```
## Loading required package: RCurl
```

```
## Loading required package: bitops
```

```
## Loading required package: XML
```

```
## Loading required package: digest
```

```
##   |                                                                         |>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>             |  79%
##   ordinary text without R code
## 
##   |                                                                         |>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>           |  83%
## label: unnamed-chunk-9
##   |                                                                         |>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>         |  86%
##   ordinary text without R code
## 
##   |                                                                         |>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>       |  90%
## label: unnamed-chunk-10
```

```
##   |                                                                         |>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    |  93%
##   ordinary text without R code
## 
##   |                                                                         |>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  |  97%
## label: unnamed-chunk-11
##   |                                                                         |>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>| 100%
##   ordinary text without R code
```

```
## output file:
## /home/cboettig/Documents/code/nonparametric-bayes/inst/examples/BUGS/allen.md
```

```
## [1] "allen.md"
```

```r

pardist <- mcmcall
# pardist[,1] = p[2] + rnorm(100, 0, 0.000001) pardist[,4] = p[1] +
# rnorm(100, 0, 0.000001) pardist[,5] = p[3] + rnorm(100, 0, 0.000001)

sdp = f_transition_matrix(f, p, x_grid, h_grid, sigma_g, pardist)
s_opt <- value_iteration(sdp, x_grid, h_grid, OptTime = 1000, xT, 
    profit, delta)
# opt <- find_dp_optim(sdp, x_grid, h_grid, OptTime, xT, profit, delta,
# reward=0)

SDP_Mat <- determine_SDP_matrix(f, p, x_grid, h_grid, sigma_g)
pars_fixed <- value_iteration(SDP_Mat, x_grid, h_grid, OptTime = 1000, 
    xT, profit, delta)


require(reshape2)
policies <- melt(data.frame(stock = x_grid, pars.uncert = x_grid[s_opt$D], 
    pars.fixed = x_grid[pars_fixed$D]), id = "stock")

ggplot(policies, aes(stock, value, color = variable)) + geom_point(alpha = 0.5) + 
    xlab("stock size") + ylab("harvest")
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-21.png) 

```r

ggplot(policies, aes(stock, stock - value, color = variable)) + geom_line(alpha = 1) + 
    xlab("stock size") + ylab("escapement")
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-22.png) 

