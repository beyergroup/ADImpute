#' @title Selection of samples for training
#'
#' @description \code{SplitData} selects a portion (\code{ratio}) of samples
#' (columns in \code{data}) to be used as training set
#'
#' @usage \code{SplitData(data, ratio = .7, training.only = T, seed = NULL)}
#'
#' @param data matrix; raw counts (genes as rows and samples as columns)
#' @param ratio numeric; ratio of the samples to be used for training
#' @param training.only logical; if TRUE define only a training dataset, if
#' FALSE writes both training and validation sets (defaults to T)
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
                      training.only = T,
                      seed = NULL){
  # For reproducibility
  if (!is.null(seed))
      set.seed(seed)
  # Randomly select samples according to given ratio
  training_samples <- sample(colnames(data))[1:round(ncol(data) * ratio)]
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



#' @title Masking of entries for performance evaluation
#'
#' @usage \code{MaskData(data, write.to.file = T, filename, mask = .1,
#' seed = NULL)}
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
                     write.to.file = T,
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
      out <- rep(F, length(x))
      out[tomask] <- T
      names(out) <- names(x)

      return(out)

    } else{

      # Not enough non-dropout samples to mask up to defined ratio, mask all
      return(x)
    }

  }

  maskidx <- t(sapply(1:nrow(maskable),
                      function(x) MaskerPerGene(maskable[x, ],
                                                rowmask = rowmask,
                                                seed = seed + x)))
  data[maskidx] <- NA
  sparcity <- sprintf("%.2f", sum(is.na(data)) / length(data))
  cat("final sparcity", sparcity, "\n")

  data[is.na(data)] <- 0

  # Write to file
  if (write.to.file)
    WriteTXT(data, filename)

  return(data)
}



#' @title Computation of MSE per gene
#'
#' @usage \code{ComputeMSEGenewise(real, masked, scimputed, baseline, network)}
#'
#' @description \code{ComputeMSEGenewise} computes the MSE of dropout
#' imputation for a given gene.
#'
#' @param real numeric; vector of original expression of a given gene (before
#' masking)
#' @param masked logical; vector indicating which entries were masked for a
#' given gene
#' @param scimputed numeric; vector of results of scImpute for a given gene
#' @param baseline numeric; vector of results of baseline imputation for a
#' given gene
#' @param network numeric; vector of results of network imputation for a given
#' gene
#'
#' @return MSE of all imputations for a given gene, for each of the three
#' methods: scImpute, average expression (baseline) and network.
#'
ComputeMSEGenewise <- function(real,
                               masked,
                               scimputed,
                               baseline,
                               network){

  index <- masked & (scimputed != 0) # do not assess performance if value is not imputed
  if (sum(index) == 0){
    scimputed_mse <- NA
  } else{
    scimputed_mse <- sum( (scimputed[index] - real[index]) ^ 2, na.rm = T) /
      sum(index, na.rm = T)
  }

  # baseline always imputes, even if with 0 - check only if there are masked values
  if (sum(masked) == 0){
    baseline_mse <- NA
  } else{
    baseline_mse <- sum( (baseline[masked] - real[masked]) ^ 2, na.rm = T) /
      sum(masked, na.rm = T)
  }

  index <- masked & (network != 0) # do not assess performance if value is not imputed
  if (sum(index) == 0){
    net_mse <- NA
  } else{
    net_mse <- sum( (network[index] - real[index]) ^ 2, na.rm = T) /
      sum(index, na.rm = T)
  }

  return(c("scimpute" = scimputed_mse,
           "baseline" = baseline_mse,
           "network"  = net_mse))
}


#' @title Method choice per gene
#'
#' @usage \code{ChooseMethod(real, masked, scimpute, baseline, network,
#' write.to.file = T)}
#'
#' @description \code{ChooseMethod} determines the method for dropout
#' imputation based on performance on each gene in training data
#'
#' @param real matrix; original gene expression data, i.e. before masking
#' (genes as rows and samples as columns)
#' @param masked matrix, logical indicating which entries were masked
#' (genes as rows and samples as columns)
#' @param scimputed matrix; results of scImpute
#' @param baseline matrix; results of baseline imputation
#' @param network matrix; results of network imputation
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
                         scimpute,
                         baseline,
                         network,
                         write.to.file = T){

  if (!all(colnames(real) == colnames(masked)) &&
       all(rownames(real) == rownames(masked)))
    break

  which_masked <- (real != 0) & (masked == 0) # distinguishes masked values from dropouts in the original data

  MSE <- t(sapply(rownames(real),
                  function(gene_name) ComputeMSEGenewise(real[gene_name, ],
                                                         which_masked[gene_name, ],
                                                         scimpute[gene_name, ],
                                                         baseline[gene_name, ],
                                                         network[gene_name, ])))
  # entries in MSE that are NAs:
  # 1) genes that are not expressed in any cell (nothing to mask) - these are
  #    NA for all
  # 2) genes not imputable by the network and with very low expr - these are
  #    imputed only by baseline, but low confidence
  # 3) genes that are not imputable either by scImpute or by the Network - NA
  #    for only that method, since baseline always imputes

  MSE <- MSE[-which(rowSums(is.na(MSE)) >= 2), ] # remove cases 1) and 2)

  cat("Imputation errors computed for", nrow(MSE), "genes\n")

  best_method <- sapply(apply(MSE, 1, which.min),
                        function(x) colnames(MSE)[x])

  if (write.to.file)
    WriteTXT(cbind(MSE, best_method), "method_choices.txt")

  return(best_method)

}
