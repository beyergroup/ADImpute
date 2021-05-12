# ADImpute predicts unmeasured gene expression values from single cell
# RNA-sequencing data (dropout imputation). This R-package combines multiple
# dropout imputation methods, including a novel gene regulatory
# network-based method.  Copyright (C) 2020 Ana Carolina Leote This program
# is free software: you can redistribute it and/or modify it under the terms
# of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.  You should have received a copy of the
# GNU General Public License along with this program.  If not, see
# <https://www.gnu.org/licenses/>.


#' @title MSE aggregation across folds
#'
#' @usage AggregateMSEArray(MSE_array, aggr.method = "median")
#'
#' @description \code{AggregateMSEArray} aggregates the results of imputation
#' across different cross-validation folds into one value of imputation error
#' per gene and imputation method
#'
#' @param MSE_array array; array of matrices with genes in rows and methods in
#' columns, where entries correspond to the the method's MSE for the specific
#' gene.
#' @param aggr.method character; method used to aggregate MSE results from
#' different folds. One of 'mean', 'median' (defaults to 'median')
#'
#' @details The MSE results from different folds are aggregated into a single
#' value per gene, according to \code{aggr.method}.
#'
#' @return matrix; matrix with genes in rows and methods in columns, where
#' entries correspond to the MSE of the method when imputing the given gene
#'
#' @seealso \code{\link{ComputeMSE}}
#'
AggregateMSEArray <- function(MSE_array, aggr.method = "median"){

    MSE <- switch(aggr.method,
        "mean" = apply(MSE_array, 2, rowMeans),
        "median" = apply(MSE_array, 1:2, stats::median))

    return(MSE)
}


#' @title Method choice per gene
#'
#' @usage ChooseMethod(MSE, write.to.file = TRUE)
#'
#' @description \code{ChooseMethod} determines the method for dropout
#' imputation based on performance on each gene in training data
#'
#' @param MSE matrix; matrix with genes in rows and methods in columns, where
#' entries correspond to the the method's MSE for the specific gene.
#' @param write.to.file logical; should the output be written to a file?
#'
#' @details The method with the best performance (lowest MSE) is chosen for each
#' gene.
#'
#' @return character; best performing method in the training set for each gene
#'
#' @seealso \code{\link{ComputeMSE}}
#'
ChooseMethod <- function(MSE, write.to.file = TRUE) {

    # keep cases where at least 2 methods are available for comparison
    MSE <- MSE[which(rowSums(!is.na(MSE)) >= 2), ]

    message(paste("Imputation errors compared for", nrow(MSE), "genes\n"))

    best_method <- vapply(apply(MSE, 1, which.min),
        function(x) colnames(MSE)[x], FUN.VALUE = "Method_name")

    if (write.to.file)
        WriteTXT(cbind(MSE, best_method), "method_choices.txt")

    return(best_method)
}


#' @title MSE computation
#'
#' @usage ComputeMSE(real, masked, imputed)
#'
#' @description \code{ComputeMSE} determines the method for dropout
#' imputation based on performance on each gene in training data
#'
#' @param real matrix; original gene expression data, i.e. before masking
#' (genes as rows and samples as columns)
#' @param masked matrix, logical indicating which entries were masked
#' (genes as rows and samples as columns)
#' @param imputed list; list of matrices with imputation results for all
#' considered methods
#'
#' @details The imputed values are compared to the real ones for every masked
#' entry in \code{real}. The Mean Squared Error is computed for all masked
#' entries per gene.
#'
#' @return matrix; matrix with genes in rows and methods in columns, where
#' entries correspond to the the method's MSE for the specific gene.
#'
#' @seealso \code{\link{ComputeMSEGenewise}}
#'
ComputeMSE <- function(real, masked, imputed){

    if (!(all(colnames(real) == colnames(masked)) &&
        all(rownames(real) == rownames(masked))))
        stop("Error! Dimnames before and after masking don't match.\n")

    if(!is.list(imputed))
        stop("Error! Imputation result for MSE computation is not a list.\n")

    which_masked <- (real != 0) & (masked == 0)  # distinguishes masked values
    # from dropouts in the original data

    MSE <- lapply(imputed, function(x) vapply(rownames(real), function(g) {
        if (g %in% rownames(x)) {
            ComputeMSEGenewise(real = real[g, ], masked = which_masked[g, ],
                imputed = x[g, ], baseline = identical(x, imputed$Baseline))
        } else {
            NA
        }
    }, FUN.VALUE = 1))
    MSE <- do.call(cbind, MSE)

    return(MSE)
}


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
ComputeMSEGenewise <- function(real, masked, imputed, baseline) {
    if (baseline) {
        index <- masked  # Baseline always imputes
        if (sum(index) == 0) {
            mse <- NA
        } else {
            mse <- sum((imputed[index] - real[index])^2, na.rm = TRUE)/
                sum(index, na.rm = TRUE)
        }
    } else {
        index <- masked & (imputed != 0)  # not assess if value is not imputed
        if (sum(index) == 0) {
            mse <- NA
        } else {
            mse <- sum((imputed[index] - real[index])^2, na.rm = TRUE)/
                sum(index, na.rm = TRUE)
        }
    }

    return(mse)
}


#' @title Preparation of training data for method evaluation
#'
#' @usage CreateTrainData(data, train.ratio = .7, train.only = TRUE, mask = .1,
#' write = FALSE)
#'
#' @description \code{CreateTrainingData} selects a subset of cells to use as
#' training set and sets a portion (\code{mask}) of the non-zero entries in each
#' row of the subset to zero
#'
#' @param data matrix; raw counts (genes as rows and samples as columns)
#' @param train.ratio numeric; ratio of the samples to be used for training
#' @param train.only logical; if TRUE define only a training dataset, if
#' FALSE writes both training and validation sets (defaults to TRUE)
#' @param mask numeric; ratio of total non-zero samples to be masked per gene
#' (defaults to .1)
#' @param write logical; should the output be written to a file?
#'
#' @return list with resulting matrix after subsetting and after masking
#'
CreateTrainData <- function(data, train.ratio = 0.7, train.only = TRUE,
                            mask = 0.1, write = FALSE) {

    # Split data for training
    train_norm <- SplitData(data, ratio = train.ratio, write.to.file = write,
        train.only = train.only)

    # Mask selected training data
    masked_train_norm <- MaskData(train_norm, write.to.file = write,
        mask = mask)

    return(list(train = train_norm, mask = masked_train_norm))
}


#' @title Preparation of training data for method evaluation
#'
#' @usage CrossValidateImputation(data, train.ratio = 0.7, train.only = TRUE,
#' mask.ratio = 0.1, do = c("Baseline", "DrImpute", "Network"), write = FALSE,
#' scale = 1, pseudo.count = 1, labels = NULL, cell.clusters = 2,
#' drop_thre = NULL, type = "count", cores = BiocParallel::bpworkers(BPPARAM),
#' BPPARAM = BiocParallel::SnowParam(type = "SOCK"),
#' net.coef = ADImpute::network.coefficients, net.implementation = "iteration",
#' tr.length = ADImpute::transcript_length, bulk = NULL, ...)
#'
#' @description \code{CreateTrainingData} selects a subset of cells to use as
#' training set and sets a portion (\code{mask}) of the non-zero entries in each
#' row of the subset to zero
#'
#' @param data matrix; raw counts (genes as rows and samples as columns)
#' @param train.ratio numeric; ratio of samples to be used for training
#' @param train.only logical; if TRUE define only a training dataset, if
#' FALSE writes and returns both training and validation sets (defaults to TRUE)
#' @param mask.ratio numeric; ratio of samples to be masked per gene
#' @param do character; choice of methods to be used for imputation. Currently
#' supported methods are \code{'Baseline'}, \code{'DrImpute'} and
#' \code{'Network'}. Not case-sensitive. Can include one or more methods. Non-
#' supported methods will be ignored.
#' @param write logical; write intermediary and imputed objects to files?
#' @param scale integer; scaling factor to divide all expression levels by
#' (defaults to 1)
#' @param pseudo.count integer; pseudo-count to be added to expression levels
#' to avoid log(0) (defaults to 1)
#' @param labels character; vector specifying the cell type of each column of
#' \code{data}
#' @param cell.clusters integer; number of cell subpopulations
#' @param drop_thre numeric; between 0 and 1 specifying the threshold to
#' determine dropout values
#' @param type A character specifying the type of values in the expression
#' matrix. Can be 'count' or 'TPM'
#' @param cores integer; number of cores used for paralell computation
#' @param BPPARAM parallel back-end to be used during parallel computation.
#' See \code{\link[BiocParallel]{BiocParallelParam-class}}.
#' @param net.coef matrix; network coefficients. Please provide if you don't
#' want to use ADImpute's network model. Must contain one first column 'O'
#' acconting for the intercept of the model and otherwise be an adjacency matrix
#' with hgnc_symbols in rows and columns. Doesn't have to be squared. See
#' \code{ADImpute::demo_net} for a small example.
#' @param net.implementation character; either 'iteration', for an iterative
#' solution, or 'pseudoinv', to use Moore-Penrose pseudo-inversion as a
#' solution. 'pseudoinv' is not advised for big data.
#' @param tr.length matrix with at least 2 columns: 'hgnc_symbol' and
#' 'transcript_length'
#' @param bulk vector of reference bulk RNA-seq, if available (average across
#' samples)
#' @param ... additional parameters from \code{EvaluateMethods}
#'
#' @return matrix; matrix with genes in rows and methods in columns, where
#' entries correspond to the the method's MSE for the specific gene.
#'
CrossValidateImputation <- function(data, train.ratio = 0.7, train.only = TRUE,
    mask.ratio = 0.1, do = c("Baseline", "DrImpute", "Network"), write = FALSE,
    scale = 1, pseudo.count = 1, labels = NULL, cell.clusters = 2,
    drop_thre = NULL, type = "count", cores = BiocParallel::bpworkers(BPPARAM),
    BPPARAM = BiocParallel::SnowParam(type = "SOCK"),
    net.coef = ADImpute::network.coefficients, net.implementation = "iteration",
    tr.length = ADImpute::transcript_length, bulk = NULL, ...){

    # create training data
    train_data <- CreateTrainData(data, train.ratio = train.ratio,
        train.only = train.only, mask = mask.ratio)

    # impute training data
    train_imputed <- Impute(data = train_data$mask, sce = NULL,
        do = do[do != "Ensemble"], write = write, outdir = getwd(),
        scale = scale, pseudo.count = pseudo.count, labels = labels,
        cell.clusters = cell.clusters, drop_thre = drop_thre, type = type,
        tr.length = tr.length, bulk = bulk, cores = cores, BPPARAM = BPPARAM,
        net.coef = net.coef, net.implementation = net.implementation, ...)

    # compute MSE on training data
    MSE <- ComputeMSE(real = round(train_data$train, 2),
        masked = round(train_data$mask, 2), imputed = train_imputed)

    return(MSE)
}


#' @title Masking of entries for performance evaluation
#'
#' @usage MaskData(data, write.to.file = FALSE, mask = .1)
#'
#' @description \code{MaskData} sets a portion (\code{mask}) of the non-zero
#' entries of each row of \code{data} to zero
#'
#' @param data matrix; raw counts (genes as rows and samples as columns)
#' @param write.to.file logical; should the output be written to a file?
#' @param mask numeric; ratio of total non-zero samples to be masked per gene
#' (defaults to .1)
#'
#' @details Sets a portion (\code{mask}) of the non-zero entries of each row of
#' \code{data} to zero. Result is written to \code{filename}.
#'
#' @return matrix containing masked raw counts (genes as rows and samples as
#' columns)
#'
MaskData <- function(data, write.to.file = FALSE, mask = 0.1) {

    message("Masking training data\n")

    data <- DataCheck_Matrix(data)

    # Original matrix sparcity
    sparcity <- sprintf("%.2f", sum(data == 0)/length(data))
    message("original sparcity ", sparcity, "\n")

    rowmask <- round(mask * ncol(data))  # samples to be masked per gene
    maskable <- data != 0  # maskable samples (originally not dropouts)

    maskidx <- t(vapply(seq_len(nrow(maskable)),
        function(x) MaskerPerGene(maskable[x, ], rowmask = rowmask),
        FUN.VALUE = rep(FALSE, ncol(data))))

    data[maskidx] <- 0
    sparcity <- sprintf("%.2f", sum(data == 0)/length(data))
    message("final sparcity ", sparcity, "\n")

    # Write to file
    if (write.to.file)
        WriteTXT(data, "masked_data.txt")

    return(data)
}


#' @title Helper mask function
#'
#' @usage MaskerPerGene(x, rowmask)
#'
#' @description Helper mask function, per feature.
#'
#' @param x logical; data to mask
#' @param rowmask numeric; number of samples to be masked per gene
#'
#' @return logical containing positions to mask
#'
MaskerPerGene <- function(x, rowmask) {
    if (sum(x) > rowmask) {
        # Randomly pick samples to mask
        tomask <- sample(which(x), rowmask)
        out <- rep(FALSE, length(x))
        out[tomask] <- TRUE
        names(out) <- names(x)

        return(out)

    } else {
        # Not enough non-dropout samples to mask up to defined ratio, mask all
        return(x)
    }
}


#' @title Selection of samples for training
#'
#' @description \code{SplitData} selects a portion (\code{ratio}) of samples
#' (columns in \code{data}) to be used as training set
#'
#' @usage SplitData(data, ratio = .7, write.to.file = FALSE, train.only = TRUE)
#'
#' @param data matrix; raw counts (genes as rows and samples as columns)
#' @param ratio numeric; ratio of the samples to be used for training
#' @param write.to.file logical; should the output be written to a file?
#' @param train.only logical; if TRUE define only a training dataset, if
#' FALSE writes both training and validation sets (defaults to TRUE)
#'
#' @details Selects a portion (\code{ratio}) of samples (columns in \code{data})
#' to be used as training set and writes to file 'training_raw.txt'.
#'
#' @return matrix containing raw counts (genes as rows and samples as columns)
#'
SplitData <- function(data, ratio = 0.7, write.to.file = FALSE,
    train.only = TRUE) {

    message("Selecting training data\n")

    # Randomly select samples according to given ratio
    train_samples <- sample(colnames(data))[seq_len(round(ncol(data) * ratio))]
    if(length(train_samples) == 0)
        stop("No samples in training data - increase train.ratio\n")
    train_data <- data[, train_samples, drop = FALSE]

    # Write to file
    if (write.to.file)
        WriteTXT(train_data, "training.txt")

    # Optionally create and write validation dataset
    if (!train.only) {
        validation_data <- data[, !(colnames(data) %in% train_samples),
            drop = FALSE]
        if (write.to.file)
            WriteTXT(validation_data, "validation.txt")
    }

    return(train_data)
}
