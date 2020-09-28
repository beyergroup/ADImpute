## ----setup, include = FALSE---------------------------------------------------
is_check <- ("CheckExEnv" %in% search()) || any(c("_R_CHECK_TIMINGS_",
             "_R_CHECK_LICENSE_") %in% names(Sys.getenv()))
knitr::opts_chunk$set(eval = !is_check, collapse = TRUE, comment = "#>")

## -----------------------------------------------------------------------------
library(ADImpute)

## -----------------------------------------------------------------------------
RPM <- NormalizeRPM(ADImpute::demo_data) # normalize for library size

## ---- warning=FALSE-----------------------------------------------------------
set.seed(1)
methods_pergene <- EvaluateMethods(data = RPM,
                                   do = c("Baseline", "DrImpute",
                                          "Network"), # these are the default
                                                      # methods to test. Exclude
                                                      # any of them by removing
                                                      # them from the vector
                                   cores = 2,
                                   train.ratio = .7, mask.ratio = .2,
                                   network.coefficients = ADImpute::demo_net)
head(methods_pergene) # show the best performing method for the first 6 genes

## ---- warning=FALSE-----------------------------------------------------------
imputed <- Impute(do = "Ensemble",
                  method.choice = methods_pergene,
                  data = RPM,
                  cores = 2,
                  network.coefficients = ADImpute::demo_net)
str(imputed)

