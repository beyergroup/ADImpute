
<!-- badges: start -->
  [![Travis build status](https://travis-ci.com/anacarolinaleote/ADImpute.svg?branch=master)](https://travis-ci.com/anacarolinaleote/ADImpute)
  [![Build status](https://ci.appveyor.com/api/projects/status/qsslj60tuvcg75vr?svg=true)](https://ci.appveyor.com/project/anacarolinaleote/adimpute)


  <!-- badges: end -->

# ADImpute
ADImpute predicts unmeasured gene expression values from single cell RNA-sequencing data (dropout imputation). This R-package combines multiple dropout imputation methods, including a novel gene regulatory network-based method.

## Getting started
ADImpute requires R version 4.0.
The following R commands should allow you to install ADImpute from github (< 10 seconds):

```
# install.packages("devtools")
devtools::install_github("anacarolinaleote/ADImpute")
library(ADImpute)
```

For an example of the usage of the ADImpute package, follow the < 10min [tutorial](https://github.com/anacarolinaleote/ADImpute/blob/master/vignettes/ADImpute_tutorial.Rmd).

ADImpute was developed in R 4.0.2, under Linux Mint 20, and tested in Linux, OS X and Windows.

For further questions, please contact: ana.carolina.leote@uni-koeln.de
