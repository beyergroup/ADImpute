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


#' @title Combine imputation methods
#'
#' @usage Combine(data, imputed, method.choice, write.to.file = FALSE)
#'
#' @param data matrix with entries equal to zero to be imputed, already
#' normalized (genes as rows and samples as columns)
#' @param method.choice named character; vector with the best performing method
#' per gene
#' @param imputed list; list of matrices with imputation results for all
#' considered methods
#' @param write.to.file logical; should a file with the imputation results be
#' written?
#'
#' @details Combines imputation results from all methods according to training
#' results provided in \code{method.choice}
#'
#' @return matrix; imputation results combining the best performing method
#' per gene
#'
#'
Combine <- function(data,
                    imputed,
                    method.choice,
                    write.to.file = FALSE){

  # all zeros are imputed
  dropouts <- data == 0

  best <- lapply(as.list(unique(method.choice)),
                 function(m) names(method.choice)[method.choice == m])
  # some genes could not have a method assigned to them, because they were too
  # lowly expressed for any masking to be done. These are assigned to Network:
  best$Network <- c(best$Network, setdiff(rownames(data), unlist(best)))

  # combine imputations for best performing methods
  combined <- matrix(nrow = nrow(data), ncol = ncol(data),
                     dimnames = dimnames(data))
  for(method in names(best))
    combined[best[[method]],] <- imputed[[method]][best[[method]],]
  combined <- DataCheck_Matrix(combined)

  # replace dropouts in the data with combined imputations
  imputed_combined <- data
  imputed_combined[dropouts] <- combined[dropouts]

  if (write.to.file)
    WriteTXT(imputed_combined, "final_imputation.txt")

  return(imputed_combined)
}


#' @title Impute using average expression across all cells
#'
#' @usage ImputeBaseline(data, write.to.file = TRUE, drop.exclude = TRUE, ...)
#'
#' @description \code{ImputeBaseline} imputes dropouts using gene averages
#' across cells
#'
#' @param data matrix with entries equal to zero to be imputed, normalized
#' and log2-transformed (genes as rows and samples as columns)
#' @param write.to.file logical; should a file with the imputation results be
#' written?
#' @param drop.exclude logical; should zeros be discarded for the calculation
#' of genewise average expression levels? (defaults to TRUE)
#' @param ... additional arguments to \code{saveRDS}
#'
#' @return matrix; imputation results considering the average expression
#' values of genes
#'
ImputeBaseline <- function(data,
                           write.to.file = TRUE,
                           drop.exclude = TRUE,
                           ...){

  cat("Imputing data using average expression\n")

  dropouts <- data == 0

  # Compute mean expression levels across cells
  if (drop.exclude){
    gene_avg <- apply(data, 1, function(x) mean(x[x != 0], na.rm = TRUE))
    gene_avg[is.na(gene_avg)] <- 0

  } else{
    gene_avg <- apply(data, 1, function(x) mean(x, na.rm = TRUE))
  }

  # Impute when dropout
  gene_avg_mat <- matrix(data = gene_avg,
                         nrow = nrow(data),
                         ncol = ncol(data),
                         byrow = FALSE,
                         dimnames = dimnames(data))
  res <- data
  res[dropouts] <- gene_avg_mat[dropouts]

  if (write.to.file){
    dir.create("Baseline")
    saveRDS(res, "Baseline/baseline_imputed.rds", ...)
  }

  return(res)
}



#' @title Use DrImpute
#'
#' @usage ImputeDrImpute(data, write.to.file = TRUE)
#'
#' @description \code{ImputeDrImpute} uses the DrImpute package for dropout
#' imputation
#'
#' @param data matrix with entries equal to zero to be imputed, normalized
#' and log2-transformed (genes as rows and samples as columns)
#' @param write.to.file logical; should a file with the imputation results be
#' written?
#'
#' @return matrix; imputation results from DrImpute
#'
#' @seealso \code{\link[DrImpute]{DrImpute}}
#'
ImputeDrImpute <- function(data, write.to.file = TRUE){

  cat("Imputing data using DrImpute\n")

  res <- DrImpute::DrImpute(as.matrix(data))
  colnames(res) <- colnames(data)

  if(write.to.file){
    dir.create("DrImpute")
    saveRDS(res, "DrImpute/drimputed.rds")
  }

  return(res)
}


#' @title Network-based imputation
#'
#' @usage ImputeNetwork(data, network.coefficients = NULL, network.path = NULL,
#' cores = 4, type = "iteration", write.to.file = TRUE,
#' drop.exclude = TRUE, ...)
#'
#' @param data matrix with entries equal to zero to be imputed, normalized
#' and log2-transformed (genes as rows and samples as columns)
#' @param network.coefficients matrix; network coefficients. Please provide
#' either \code{network.coefficients} or \code{network.path}.
#' @param network.path character; path to .txt or .rds file with network
#' coefficients
#' @param cores integer; number of cores to use
#' @param type character; either "iteration", for an iterative solution, or
#' "pseudoinv", to use Moore-Penrose pseudo-inversion as a solution.
#' @param write.to.file logical; should a file with the imputation results be
#' written?
#' @param drop.exclude logical; should zeros be discarded for the calculation
#' of genewise average expression levels? (defaults to TRUE)
#' @param ... additional arguments to \code{ImputeNetParallel}
#'
#' @details Imputes dropouts using a gene regulatory network trained on external
#' data, as provided in \code{network.path}. Dropout expression values are
#' estimated from the expression of their predictor genes and the network
#' coefficients.
#'
#' @return matrix; imputation results incorporating network information
#'
#' @seealso \code{\link{ImputeNetParallel}}
#'
ImputeNetwork <- function(data,
                          network.coefficients = NULL,
                          network.path = NULL,
                          cores = 4,
                          type = "iteration",
                          write.to.file = TRUE,
                          drop.exclude = TRUE,
                          ...){

  cat("Imputing data using network information\n")

  # Check arguments
  Check <- CreateArgCheck(match = list("type" = type),
                          acceptable = list("type"=c("iteration","pseudoinv")))
  ArgumentCheck::finishArgCheck(Check)

  # Limit data and network to genes common to both
  arranged <- ArrangeData(data, network.path, network.coefficients)

  cat("Data dim:", paste(dim(arranged$data), collapse = " x "), "\n")
  cat("Network dim:", paste(dim(arranged$network), collapse = " x "), "\n")

  # Center expression of each gene
  centered <- CenterData(arranged$data, drop.exclude)
  arranged$centered <- as.matrix(centered$data)

  dropout_mat <- arranged$data == 0 # dropout indexes in the data matrix

  cat("Starting network-based imputation\n")
  new_imp <- ImputeNetParallel(dropout_mat, arranged, cores, type, ...)

  res <- data
  # add back gene means
  res[rownames(new_imp), ][dropout_mat[rownames(new_imp), ]] <-
    (centered$center[rownames(new_imp)] + new_imp)[
      dropout_mat[rownames(new_imp), ]]

  res[res < 0] <- 0

  if (write.to.file){
    dir.create("Network")
    saveRDS(res, "Network/network_imputed.rds")
  }

  return(res)
}



#' @title Use SAVER
#'
#' @usage ImputeSAVER(data, cores, try.mean = FALSE, write.to.file = TRUE)
#'
#' @description \code{ImputeSAVER} uses the SAVER package for dropout
#' imputation
#'
#' @param data matrix with entries equal to zero to be imputed, normalized
#' (genes as rows and samples as columns)
#' @param cores integer; number of cores to use
#' @param try.mean logical; whether to additionally use mean gene expression
#' as prediction
#' @param write.to.file logical; should a file with the imputation results be
#' written?
#'
#' @return matrix; imputation results from SAVER
#'
#' @seealso \code{\link[SAVER]{saver}}
#'
ImputeSAVER <- function(data, cores, try.mean = FALSE, write.to.file = TRUE){

  cat("Imputing data using SAVER\n")

  if(write.to.file)
    dir.create("SAVER")

  if(try.mean){
    imp_mean <- SAVER::saver(data, size.factor = 1, ncores = cores,
                             null.model = TRUE)
    if(write.to.file)
      saveRDS(object = imp_mean, file = "SAVER/SAVER_nullmodel.rds")
  }

  res <- SAVER::saver(data, ncores = cores, size.factor = 1)

  if (write.to.file){
    dir.create("SAVER")
    saveRDS(object = res, file = "SAVER/SAVER.rds")
  }

  return(res$estimate)
}
