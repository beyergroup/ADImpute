# ADImpute predicts unmeasured gene expression values from single cell
# RNA-sequencing data (dropout imputation). This R-package combines multiple
# dropout imputation methods, including a novel gene regulatory network-based
# method.
# Copyright (C) 2020  Ana Carolina Leote
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


#' @title Computation of MSE per gene
#'
#' @usage ComputeMSEGenewise(real, masked, imputed, baseline)
#'
#' @description \code{ComputeMSEGenewise} computes the MSE of dropout
#' imputation for a given gene.
#'
#' @param real numeric; vector of original expression of a given gene (before
#' masking)
#' @param masked logical; vector indicating which entries were masked for a
#' given gene
#' @param imputed matrix; imputation results for a given imputation method
#' @param baseline logical; is this baseline imputation?
#'
#' @return MSE of all imputations indicated by \code{masked}
#'
ComputeMSEGenewise <- function(real,
                               masked,
                               imputed,
                               baseline){

  if(baseline){
    index <- masked # Baseline always imputes
    if(sum(index) == 0){
      mse <- NA
    } else{
      mse <- sum((imputed[index]-real[index])^2,
                 na.rm = TRUE)/sum(index, na.rm = TRUE)
    }
  } else{
    index <- masked & (imputed != 0) # do not assess if value is not imputed
    if(sum(index) == 0){
      mse <- NA
    } else{
      mse <- sum((imputed[index]-real[index])^2,
                 na.rm = TRUE)/sum(index, na.rm = TRUE)
    }
  }

  return(mse)
}


#' @title Method choice per gene
#'
#' @usage ChooseMethod(real, masked, imputed, write.to.file = TRUE)
#'
#' @description \code{ChooseMethod} determines the method for dropout
#' imputation based on performance on each gene in training data
#'
#' @param real matrix; original gene expression data, i.e. before masking
#' (genes as rows and samples as columns)
#' @param masked matrix, logical indicating which entries were masked
#' (genes as rows and samples as columns)
#' @param imputed list; list of matrices with imputation results for all
#' considered methods
#' @param write.to.file logical; should the output be written to a file?
#'
#' @details The imputed values are compared to the real ones for every masked
#' entry in \code{real}. The Mean Squared Error
#' is computed for all masked entries per gene and the method with the best
#' performance is chosen for each gene.
#'
#' @return character; best performing method in the training set for each gene
#'
#' @seealso \code{\link{ComputeMSEGenewise}}
#'
ChooseMethod <- function(real,
                         masked,
                         imputed,
                         write.to.file = TRUE){

  if (!all(colnames(real) == colnames(masked)) &&
      all(rownames(real) == rownames(masked)))
    stop("Error! Colnames / rownames before and after masking don't match.\n")

  which_masked <- (real != 0) & (masked == 0) # distinguishes masked values from
                                              # dropouts in the original data

  MSE <- lapply(imputed,
                function(x) sapply(rownames(real),
                                   function(gene_name) {if(gene_name %in% rownames(x)){
                                     ComputeMSEGenewise(real = real[gene_name,],
                                                        masked = which_masked[gene_name,],
                                                        imputed = x[gene_name,],
                                                        baseline = identical(x, imputed$Baseline))
                                            } else{ NA }}))

  MSE <- do.call(cbind, MSE)

  # keep cases where at least 2 methods are available for comparison
  MSE <- MSE[which(rowSums(!is.na(MSE)) >= 2), ]

  cat("Imputation errors computed for", nrow(MSE), "genes\n")

  best_method <- sapply(apply(MSE, 1, which.min),
                        function(x) colnames(MSE)[x])

  if (write.to.file)
    WriteTXT(cbind(MSE, best_method), "method_choices.txt")

  return(best_method)
}


#' @title Masking of entries for performance evaluation
#'
#' @usage MaskData(data, write.to.file = TRUE, filename, mask = .1,
#' seed = NULL)
#'
#' @description \code{MaskData} sets a portion (\code{mask}) of the non-zero
#' entries of each row of \code{data} to zero
#'
#' @param data matrix; raw counts (genes as rows and samples as columns)
#' @param write.to.file logical; should the output be written to a file?
#' @param filename character; name of output file to write masked data to
#' @param mask numeric; ratio of total non-zero samples to be masked per gene
#' (defaults to .1)
#' @param seed integer; optional seed (defaults to NULL, random selection of
#' samples to mask)
#'
#' @details Sets a portion (\code{mask}) of the non-zero entries of each row of
#' \code{data} to zero. For reproducible results, provide seed.
#' \code{seed} is summed to row number so that every gene has a different seed
#' (different samples are masked for each gene).
#' Result is written to \code{filename}.
#'
#' @return matrix containing masked raw counts (genes as rows and samples as
#' columns)
#'
MaskData <- function(data,
                     write.to.file = TRUE,
                     filename,
                     mask = .1,
                     seed = NULL){

  # Convert zeros to NAs
  data[data == 0] <- NA

  data <- as.matrix(data)

  # Original matrix sparcity
  sparcity <- sprintf("%.2f", sum(is.na(data)) / length(data))
  cat("original sparcity", sparcity, "\n")

  rowmask <- round(mask * ncol(data)) # samples to be masked per gene
  maskable <- !is.na(data) # maskable samples (originally not dropouts)

  MaskerPerGene <- function(x, rowmask, seed){

    if (sum(x) > rowmask){

      # For reproducibility
      if (!is.null(seed))
        set.seed(seed)

      # Randomly (or not) pick samples to mask
      tomask <- sample(which(x), rowmask)
      out <- rep(FALSE, length(x))
      out[tomask] <- TRUE
      names(out) <- names(x)

      return(out)

    } else{

      # Not enough non-dropout samples to mask up to defined ratio, mask all
      return(x)
    }

  }

  maskidx <- t(sapply(seq_len(nrow(maskable)),
                      function(x)
                        if(!is.null(seed)){
                          MaskerPerGene(maskable[x, ],
                                        rowmask = rowmask,
                                        seed = seed + x)}else{
                                          MaskerPerGene(maskable[x,],
                                                        rowmask = rowmask,
                                                        seed = seed)}))
  data[maskidx] <- NA
  sparcity <- sprintf("%.2f", sum(is.na(data)) / length(data))
  cat("final sparcity", sparcity, "\n")

  data[is.na(data)] <- 0

  # Write to file
  if (write.to.file)
    WriteTXT(data, filename)

  return(data)
}


#' @title Selection of samples for training
#'
#' @description \code{SplitData} selects a portion (\code{ratio}) of samples
#' (columns in \code{data}) to be used as training set
#'
#' @usage SplitData(data, ratio = .7, training.only = TRUE, seed = NULL)
#'
#' @param data matrix; raw counts (genes as rows and samples as columns)
#' @param ratio numeric; ratio of the samples to be used for training
#' @param training.only logical; if TRUE define only a training dataset, if
#' FALSE writes both training and validation sets (defaults to TRUE)
#' @param seed integer; optional seed (defaults to NULL, random selection of
#' samples to use for training)
#'
#' @details Selects a portion (\code{ratio}) of samples (columns in \code{data})
#' to be used as training set and writes to file "training_raw.txt".
#' For reproducible results, provide seed.
#'
#' @return matrix containing raw counts (genes as rows and samples as columns)
#'
SplitData <- function(data,
                      ratio = .7,
                      training.only = TRUE,
                      seed = NULL){
  # For reproducibility
  if (!is.null(seed))
    set.seed(seed)
  # Randomly select samples according to given ratio
  training_samples <- sample(colnames(data))[seq_len(round(ncol(data) * ratio))]
  training_data    <- data[, training_samples]

  # Write to file
  WriteTXT(training_data, "training.txt")

  # Optionally create and write validation dataset
  if (!training.only){
    validation_data <- data[, !(colnames(data) %in% training_samples)]
    WriteTXT(validation_data, "validation.txt")
  }

  return(training_data)
}
