
<!-- badges: start -->
[![Build Status](https://travis-ci.com/anacarolinaleote/ADImpute.svg?branch=master)](https://travis-ci.com/anacarolinaleote/ADImpute)
  [![Build status](https://ci.appveyor.com/api/projects/status/qsslj60tuvcg75vr?svg=true)](https://ci.appveyor.com/project/anacarolinaleote/adimpute)
[![Codecov test coverage](https://codecov.io/gh/anacarolinaleote/ADImpute/branch/master/graph/badge.svg)](https://codecov.io/gh/anacarolinaleote/ADImpute?branch=master)
<!-- badges: end -->

# ADImpute
ADImpute predicts unmeasured gene expression values from single cell
RNA-sequencing data (dropout imputation). This R-package combines multiple
dropout imputation methods. ADImpute currently supports, by default,
```DrImpute``` and ```SAVER``` and two novel imputation methods:
```Baseline``` and ```Network```. ```Baseline``` imputes dropouts with the
average quantified expression level of the respective gene across the dataset.
```Network``` uses previously learnt regulatory models of gene expression to
infer the expression of dropout genes from the expression of other relevant
(predictive) genes in the cell.
For more details refer to the [preprint](https://www.biorxiv.org/content/10.1101/611517v2).

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
The ```Impute()``` function allows for dropout imputation with one or more of
the supported methods. In order to specify which imputation method(s) to run,
pass them via the ```do``` argument:
```
RPM <- ADImpute::NormalizeRPM(ADImpute::demo_data)
imputed <- Impute(data = RPM, do = c("Network"), cores = 2,
    net.coef = ADImpute::demo_net)
```


### Imputation with ensemble
In addition to running different methods on the data, ADImpute can also
determine which of these performs best for each gene and perform an "Ensemble"
imputation, which combines the best performing methods for different genes.
First, evaluate methods using ```EvaluateMethods``` to determine the best
performing imputation method for each gene. This step sets a fraction of the
quantified entries in the input data to zero, applies different imputation
methods to the data and compares the imputation results to the original values.
This allows ADImpute to determine which method imputes values with the lowest
errors for each gene.
```
RPM <- ADImpute::NormalizeRPM(ADImpute::demo_data)
methods_pergene <- EvaluateMethods(data = RPM,
    do = c("Baseline", "DrImpute", "Network"),
    cores = 2, net.coef = ADImpute::demo_net)
```
After determining which method performs best for each gene, the imputation can
be re-done on the original data and the results of different methods combined
into an ensemble:
```
imputed <- Impute(do = "Ensemble", method.choice = methods_pergene,
    data = RPM, cores = 2, net.coef = ADImpute::demo_net)
```
Both the method-specific imputations and the final ensemble results are
available for further examination.


### Determination of biological zeros
Some zeros in the data correspond to genes expressed in the cell, but not
captured upon sequencing - the technical dropouts - while others correspond to
genes truly not expressed in the cell - the biological zeros. In order to avoid
imputation of biological zeros, ```ADImpute``` adapts the well-established
approach of ```scImpute``` for the computation of the probability of each entry
to be a technical dropout. A matrix of such probabilities, of the same size as
the original data, can be provided by the user, or computed by ```ADImpute```
using ```scImpute```'s approach, as below. To activate this option, provide a
value for ```true.zero.thr``` in the call to ```Impute()```, as exemplified
below:
```
imputed <- Impute(do = "Baseline",
    data = RPM,
    cores = 2,
    true.zero.thr = .3)
```


### Imputation of a SingleCellExperiment
```ADImpute``` can also take a ```SingleCellExperiment``` object as input.
In this case, ```EvaluateMethods()``` will result in new internal row metadata
being added to the ```SingleCellExperiment``` object, with the best performing
methods per gene. ```Impute()``` results in new assays being added to the
object. If ```true.zero.thr``` is specified, only the results after setting
biological zeros back to zero will be added to the ```SingleCellExperiment```
object.
```
sce <- NormalizeRPM(sce = ADImpute::demo_sce)
sce <- EvaluateMethods(sce = sce)
sce <- Impute(sce = sce)
```

### Additional imputation methods
```ADImpute``` is built in a modular way, which facilitates the addition of
custom functions supporting other imputation methods. Two such methods are
```scImpute``` and ```SCRABBLE```, with wrapper functions already contained
within ```ADImpute```. To call these methods, please follow these steps:
1) install scImpute and/or SCRABBLE from their github repositories
2) clone the ADImpute repository
3) copy the lines below to the file Wrap.R in the source R/ folder of ADImpute,
line #309.
4) re-load ADImpute using devtools::load_all() on ADImpute's folder
```
# # call to scImpute
if('scimpute' %in% tolower(do)){
    message('Make sure you have previously installed scImpute via GitHub.\n')
    res <- tryCatch(ImputeScImpute(count_path, labeled = is.null(labels),
            Kcluster = cell.clusters, labels = labels, drop_thre = drop_thre,
            cores = cores, type = type, tr.length = tr.length),
        error = function(e){ stop(paste('Error:', e$message,
            '\nTry sourcing the Impute_extra.R file.'))})
    imputed$scImpute <- log2( (res / scale) + pseudo.count)
}

# call to SCRABBLE
if('scrabble' %in% tolower(do)){
    message('Make sure you have previously installed SCRABBLE via GitHub.\n')
    res <- tryCatch(ImputeSCRABBLE(data, bulk),
                    error = function(e) { stop(paste('Error:', e$message,
                        '\nTry sourcing the Impute_extra.R file.'))})
    imputed$SCRABBLE <- log2( (res / scale) + pseudo.count)
    rm(res);gc()
}
```
After these steps, ```scImpute``` and ```SCRABBLE``` can be run with
```EvaluateMethods``` or ```Impute()``` with
```do = c("scImpute","SCRABBLE")```.

