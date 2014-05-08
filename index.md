---
layout: page
---

nonparametric-bayes
===================


### Collaborators

- Carl Boettiger
- Marc Mangel
- Steve Munch

### Funding & Support

- NSF [DBI-1306697](http://www.nsf.gov/awardsearch/showAward?AWD_ID=1306697) to CB
- Center for Stock Assessment Research, Southwest Fisheries Science Center, NOAA
- Department of Applied Mathematics and Statistics, Baskins School of Engineering, UC Santa Cruz


### Abstract

Decision-theoretic methods often rely on simple parametric models of ecological dynamics to compare the value of a potential sequence of actions. Unfortunately, such simple models rarely capture the complexity or uncertainty found in most real ecosystems. Non-parametric Bayesian methods offer a promising statistical approach for predictive modeling of ecological dynamics in regions of state space where the data is adequate, while at the same time offering more flexible patterns with greater uncertainty outside the observed data. This contrasts from simple parametric models which provide relatively constant level of uncertainty in regions with and without adequate data. The consequence of such misplaced confidence outside the data can lead to highly undesirable results that may be avoided with the more flexible non-parametric Bayesian approach.

### Publications

_in preparation, see [manuscripts](http://github.com/cboettig/nonparametric-bayes/tree/master/manuscripts) directory_

### Software

Source code for functions developed for this research is provided as an R package on the project's Github repository, [cboettig/nonparametric-bayes](http://github.com/cboettig/nonparametric-bayes). The package can be installed directly from `R` using:

```r
library(devtools)
install_github("nonparameteric-bayes", "cboettig")
```


### Data

Any necessary raw data ships with the package.  Results from the analyses are provided as binary R cache objects.

### Lab Notebook

* See my lab notebook entries under the [nonparametric-bayes](http://www.carlboettiger.info/tags.html#nonparametric-bayes) tag for ongoing discription of this research.
* The [issues tracker](https://github.com/cboettig/nonparametric-bayes/issues) lists both current and accomplished goals in this project, and steps towards their completion.
* See the repository [history](https://github.com/cboettig/nonparametric-bayes/commits/master) for a fine-grained view of progress and changes. The repository master branch is structured as an R package and can be installed directly from github using the `devtools` package.  Active research scripts are found on the `gh-pages` branch (previously in `inst/examples` directory of the master branch). Active scripts can be seen here:


Current scripts
---------------

({{ site.pages | size }} entries)

{% for post in site.pages %}
- [{{ post.path }}]({{ post.url }})
{% endfor %}

