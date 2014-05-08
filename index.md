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

### Code

### Data

### Notebook



A repository for my research in nonparametric Bayesian inference for improving ecosystem management under substantial structural uncertainty.  This work is being conducted as part of my postdoctoral research program under Marc Mangel and Steve Munch at UC Santa Cruz and the NOAA SW Marine Fisheries Service in Santa Cruz.

* See my lab notebook entries under the [nonparametric-bayes](http://www.carlboettiger.info/tags.html#nonparametric-bayes) tag for ongoing discription of this research.
* The [issues tracker](https://github.com/cboettig/nonparametric-bayes/issues) lists both current and accomplished goals in this project, and steps towards their completion.
* See the repository [history](https://github.com/cboettig/nonparametric-bayes/commits/master) for a fine-grained view of progress and changes

This repository is structured as an R package and can be installed directly from github using the `devtools` package.  Active research scripts are found on the `gh-pages` branch.


Current scripts
---------------

({{ site.pages | size }} entries)

{% for post in site.pages %}
- [{{ post.path }}]({{site.repo}}{{ post.url }})
{% endfor %}

