---
output: github_document
---

<!-- badges: start -->
  [![Travis build status](https://travis-ci.com/anacarolinaleote/ADImpute.svg?branch=master)](https://travis-ci.com/anacarolinaleote/ADImpute)
  <!-- badges: end -->

# ADImpute
ADImpute predicts unmeasured gene expression values from single cell RNA-sequencing data (dropout imputation). This R-package combines multiple dropout imputation methods, including a novel gene regulatory network-based method.

## Getting started
ADImpute requires R version 3.4.4 and the R packages scImpute, SCRABBLE and DrImpute to be installed.
The following R commands should allow you to install ADImpute from github (< 10 seconds):

```
# install.packages("devtools")
devtools::install_github("anacarolinaleote/ADImpute")
library(ADImpute)
```

For an example of the usage of the ADImpute package, follow the 20-minute [tutorial](https://github.com/anacarolinaleote/ADImpute/blob/master/vignettes/ADImpute_tutorial.Rmd).

ADImpute was developed and tested in R 3.4.4 on Fedora release 27.

For further questions, please contact: ana.carolina.leote@uni-koeln.de
