Source for: Avoiding tipping points in fisheries management through Gaussian Process Dynamic Programming
========================================================================================================

[![Build Status](http://107.170.225.143:88/api/badge/github.com/cboettig/nonparametric-bayes/status.svg?branch=master)](http://107.170.225.143:88/github.com/cboettig/nonparametric-bayes)

<!--
Manuscript code is run automatically any time changes are made to
this repository.  The badge above summarizes the status and contains
links to further details. 
--> 

Quickstart
----------

[Install Docker](https://docs.docker.com/installation) on your laptop or cloud server. In Docker, run: 

```bash
docker run -d -p 8787:8787 cboettig/nonparametric-bayes
```

This downloads the computational environment necessary and launches RStudio-server to interact with it. Point your browser to:

```
http://localhost:8787
```

Mac and Windows users should replace `localhost` with the IP address returned by `boot2docker ip`.  Log in with user/pw rstudio/rstudio. 

Open the `manuscript.Rmd` file from the manuscripts directory and you're ready to explore.  See [rocker-org Wiki](https://github.com/rocker-org/rocker/wiki/Using-the-RStudio-image) for details like custom passwords or user names.  



Edit the corresponding `.Rmd` file (manuscript or supplement),
which contains both text and code.


Installation requirements
-------------------------

Compiling the manuscript requires additional software beyond the basic
requirements of the `nonparametric-bayes` package, such as `jags` (for
comparison to parametric models), but also sofware such as `LaTeX` and
`pandoc` for typesetting. If you run into trouble, use the Docker approach
described above instead, which provides all these components already installed.

- RStudio (includes `pandoc`)
- A LaTeX environment
- The `nonparametric-bayes` package [from github](http://github.com/cboettig/nonparametric-bayes), including all packages on the SUGGESTS list. 


Manuscripts in this package are written in [knitr]'s R Markdown format
(`.Rmd` files), and compiled into PDFs using [pandoc] with [LaTeX]. Thanks
to pandoc, this workflow can easily produce alternative formats, such
as Microsoft Office's `.docx` (or open standard `.odt`), HTML, various
flavors of `.tex` (for journal submission engines) or pure markdown
(e.g. to be rendered by Github). Because markdown is platform independent
plain text and easier to learn than LaTeX, the hope is that this workflow
is relatively portable across users.  Markdown also does a better job than
most alternatives (tex included) at separating content from formatting,
freeing the writer to just write.


[knitr]: http://yihui.name/knitr
[pandoc]: http://johnmacfarlane.net/pandoc/
[LaTeX]: http://www.latex-project.org/


Caching
-------

I've also enabled caching.  It can be annoying to have to have to rerun
all the R code just to make a textual change to the manuscript or readme.
The cache is ignored by git so the first time you run `make` all the
code will run, and thereafter you will have the cache. The caching is
modestly intelligent, in that if you edit a chunk it will be rerun by
default, (as will chunks with declared dependencies on it). See `knitr`'s
[caching documentation] for details.

You can clear the cache by deleting it, or just use `make clear-cache`.

**Restoring the default cache**

You can obtain my current cache by using `make restore-cache`.

**Picking up at a particular chunk**

For the sake of modular reproducibility, it can also be desirable
to pick up from somewhere in the middle of the manuscript and not
want to have to run all the previous code just to try out one line
(not really an issue in this paper, but more generally).  Caching is
therefore modular by chunk, allowing you to restore the results
of a particular chunk to investigate. Again see `knitr`'s [caching
documentation] for details.


[caching documentation]: http://yihui.name/knitr/demo/cache/
