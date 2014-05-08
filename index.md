---
layout: page
---

# nonparametric-bayes

* Author: Carl Boettiger

A repository for my research in nonparametric Bayesian inference for improving ecosystem management under substantial structural uncertainty.  This work is being conducted as part of my postdoctoral research program under Marc Mangel and Steve Munch at UC Santa Cruz and the NOAA SW Marine Fisheries Service in Santa Cruz.

* See my lab notebook entries under the [nonparametric-bayes](http://www.carlboettiger.info/tags.html#nonparametric-bayes) tag for ongoing discription of this research.
* The [issues tracker](https://github.com/cboettig/nonparametric-bayes/issues) lists both current and accomplished goals in this project, and steps towards their completion.
* See the repository [history](https://github.com/cboettig/nonparametric-bayes/commits/master) for a fine-grained view of progress and changes

This repository is structured as an R package and can be installed directly from github using the `devtools` package.  Active research scripts are found on the `gh-pages` branch and can be browsed here:

### {{ site.pages | size }} entries

{% for post in site.pages %}
- [{{ post.path }}]({{ post.url }})
{% endfor %}

