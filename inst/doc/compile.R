compile <- function(target, fmt = list("html", "pdf", "docx", "tex", "epub")){
  fmt <- match.arg(fmt)
  library(knitr)


## global chunk options: default uses high quality Cairo PNG device
  opts_chunk$set(cache=TRUE, fig.path='figure/', dev='Cairo_png',
    fig.width=5, fig.height=5, cache.path = 'cache-local/', par=TRUE)
#  opts_chunk$set(dev.args=list(bg="transparent")) # is this needed for Cario_png?
  
## xtable html format for github/html output.  
## otherwise this should be ignored (defaults to 'latex') for pdf output
  options(xtable.type = 'html')

## Global options 
  opts_chunk$set(warning=FALSE, message=FALSE, comment=NA, tidy=FALSE, echo=FALSE)

## verbose compile
  opts_knit$set(progress = TRUE, verbose = TRUE)


    if (fmt %in% c('github', 'html')) {
      # upload to flickr when output is for github or html
      opts_knit$set(upload.fun = socialR::notebook.url)
      opts_chunk$set(cache.path = 'cache-upload/')
    } else if (fmt == 'epub') {

    } else if (fmt %in% c('odt', 'docx', 'doc')) {
      options(xtable.type = 'html')
      # xtable shouldn't print comments 
      options(xtable.print.comment=FALSE)
      # Journal probably wants eps formatted graphics.  
      opts_chunk$set(dev = 'postscript', fig.width=5, fig.height=5, 
                     cache.path = 'cache-doc/')
      # try splitting the tables out with special options?
      options(knitr.include=TRUE)
      options(knitr.split=TRUE)
      # For Word docs, we'll want plots to be external only:
      knit_hooks$set(plot = function(x, options) "") 

    } else if (fmt %in% c('tex','pdf')) {
      options(xtable.print.comment=FALSE)
      options(xtable.type = 'latex', table.placement="H")
      opts_chunk$set(cache.path = 'cache-pdf/')
      # use pdf graphics for PDF output
      opts_chunk$set(dev = 'Cairo_pdf', fig.width=5, fig.height=5)
    }
      knit(target)
}
