
<!-- badges: start -->
  [![Travis build status](https://travis-ci.com/anacarolinaleote/ADImpute.svg?branch=master)](https://travis-ci.com/anacarolinaleote/ADImpute)
  [![Build status](https://ci.appveyor.com/api/projects/status/qsslj60tuvcg75vr?svg=true)](https://ci.appveyor.com/project/anacarolinaleote/adimpute)
<!-- badges: end -->

# ADImpute
ADImpute predicts unmeasured gene expression values from single cell RNA-sequencing data (dropout imputation). This R-package combines multiple dropout imputation methods, including a novel gene regulatory network-based method.

ADImpute was developed in R 4.0.2, under Linux Mint 20, and tested in Linux, OS X and Windows.
For further questions, please contact: ana.carolina.leote@uni-koeln.de

## Installation
ADImpute requires R version 4.0.
The following R commands should allow you to install ADImpute from github (< 10 seconds):
```
# install.packages("devtools")
devtools::install_github("anacarolinaleote/ADImpute")
library(ADImpute)
```

## Quick start

### Imputation with method(s) of choice
ADImpute currently supports, by default, DrImpute and two novel imputation
methods: a baseline method, "Baseline", where genes are imputed with their
average quantified expression across the dataset, and a regulatory-network-
-based method, which uses previously learnt regulatory models of gene expression
to infer the expression of dropout genes from the expression of other relevant
genes in the cell.
```
RPM <- ADImpute::NormalizeRPM(ADImpute::demo_data)
imputed <- Impute(data = RPM,
                  do = c("Baseline","Network","DrImpute"),
                  cores = 2,
                  net.coef = ADImpute::demo_net)
```

### Imputation with ensemble
In addition to running different methods on the data, ADImpute can also
determine which of these performs best for each gene and perform an "Ensemble"
imputation, which combines the best performing methods for different genes.
First, evaluate methods to determine the best performing imputation method for
each gene. This step sets a fraction of the quantified entries in the input data
to zero, applies different normalization methods to the data and compares the
imputation results to the original values. This allows ADImpute to determine
which method imputes values with the lowest errors for each gene.
```
RPM <- ADImpute::NormalizeRPM(ADImpute::demo_data)
methods_pergene <- EvaluateMethods(data = RPM,
                                   do = c("Baseline", "DrImpute",
                                          "Network"), # these are the default
                                                      # methods to test. Exclude
                                                      # any of them by removing
                                                      # them from the vector
                                   cores = 2,
                                   train.ratio = .7, mask.ratio = .2,
                                   net.coef = ADImpute::demo_net)
```
After determining which method performs best for each gene, ADImpute redoes the
imputation on the original data and combines the results of different methods
into an ensemble.
```
imputed <- Impute(do = "Ensemble",
                  method.choice = methods_pergene,
                  data = RPM,
                  cores = 2,
                  net.coef = ADImpute::demo_net)
```

### Determination of biological zeros
Some zeros in the data correspond to genes expressed in the cell, but not
captured upon sequencing - the technical dropouts - while others correspond to
genes truly not expressed in the cell - the biological zeros. In order to avoid
imputation of biological zeros, ADImpute adapts the well-established approach of
scImpute for the computation of the probability of each entry to be a technical
dropout. A matrix of such probabilities, of the same size as the original data,
can be provided by the user, or computed by ADImpute using scImpute's approach,
as below. To activate this option, provide a value for true.zero.thr in the call
to Impute(), as exemplified below:
```
imputed <- Impute(do = "Baseline",
                  data = RPM,
                  cores = 2,
                  true.zero.thr = .3)
```

### Imputation with SCRABBLE and scImpute
In order to use SCRABBLE and/or scImpute imputation, please follow these steps:
1) install scImpute and/or SCRABBLE from their github repositories
2) clone the ADImpute repository
2) open the file Impute_extra.R in the source R/ folder of ADImpute
3) copy the commented lines at the end of the script to the file Wrap.R in the
source R/ folder of ADImpute, line #276.
4) re-load ADImpute using devtools::load_all() on ADImpute's folder

