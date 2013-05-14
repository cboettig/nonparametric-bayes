One-step ahead predictor plots




```r
opts_chunk$set(tidy = FALSE, warning = FALSE, message = FALSE, cache = FALSE, 
    comment = NA, fig.width = 6, fig.height = 4)

library(knitcitations)

library(ggplot2)  # plotting
opts_knit$set(upload.fun = socialR::flickr.url)
theme_set(theme_bw(base_size = 12))
theme_update(panel.background = element_rect(fill = "transparent", colour = NA), 
    plot.background = element_rect(fill = "transparent", colour = NA))
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
    "#D55E00", "#CC79A7")
```




```r

gp_f_at_obs <- gp_predict(gp, x, burnin=1e4, thin=300)


step_ahead <- function(x, f, p){
  h = 0
  x_predict <- sapply(x, f, h, p)
  n <- length(x_predict) - 1
  y <- c(x[1], x_predict[1:n])
  y
}

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

![plot of chunk unnamed-chunk-1](http://farm8.staticflickr.com/7286/8739510164_c2edb036c6_o.png) 

