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


#' @title Data trimming
#'
#' @usage ArrangeData(data, net.coef = NULL)
#'
#' @description \code{ArrangeData} finds common genes to the network and
#' provided data and limits both datasets to these
#'
#' @param data matrix with entries equal to zero to be imputed (genes as rows
#' and samples as columns)
#' @param net.coef matrix; object containing network coefficients
#'
#' @return list; data matrix, network coefficients matrix and intercept for
#' genes common between the data matrix and the network
#'
ArrangeData <- function(data, net.coef = NULL) {

    if (is.null(data))
        stop("Please provide an input data matrix.\n")

    if (is.null(net.coef))
        stop("Please provide valid network coefficients.\n")

    data <- DataCheck_Matrix(data)
    net.coef <- DataCheck_Network(net.coef)

    O <- net.coef[, 1]  # network intercept
    network_matrix <- net.coef[, -1]  # network coefficients

    comm_targ <- intersect(rownames(network_matrix), rownames(data))
    comm_pred <- intersect(colnames(network_matrix), rownames(data))

    comm <- union(comm_targ, comm_pred)

    network_matrix <- round(network_matrix[comm_targ, comm_pred], 2)
    O <- round(O[comm_targ], 2)
    data <- data[comm, ]

    return(list(data = data, network = network_matrix, O = O))
}


#' @title Data centering
#'
#' @usage CenterData(data)
#'
#' @description \code{CenterData} centers expression of each gene at 0
#'
#' @param data matrix of gene expression to be centered row-wise (genes as rows
#' and samples as columns)
#'
#' @return list; row-wise centers and centered data
#'
CenterData <- function(data) {

    message("Centering expression of each gene at 0\n")


    center <- apply(data, 1, function(x) mean(x[x != 0], na.rm = TRUE))
    center[is.na(center)] <- 0
    center <- round(center, 2)

    names(center) <- rownames(data)

    data <- data - center

    return(list(center = center, data = data))
}


#' @title Network-based parallel imputation
#'
#' @usage ImputeNetParallel(drop.mat, arranged, cores =
#' BiocParallel::bpworkers(BPPARAM), type = 'iteration', max.iter = 50,
#' BPPARAM = BiocParallel::SnowParam(type = "SOCK"))
#' #'
#' @description \code{ImputeNetParallel} implements network-based imputation
#' in parallel
#'
#' @param drop.mat matrix, logical; dropout entries in the data matrix
#' (genes as rows and samples as columns)
#' @param arranged list; output of \code{\link{ArrangeData}}
#' @param cores integer; number of cores used for paralell computation
#' @param type character; either 'iteration', for an iterative solution, or
#' 'pseudoinv', to use Moore-Penrose pseudo-inversion as a solution.
#' @param max.iter numeric; maximum number of iterations for network
#' imputation. Set to -1 to remove limit (not recommended)
#' @param BPPARAM parallel back-end to be used during parallel computation.
#' See \code{\link[BiocParallel]{BiocParallelParam-class}}.
#'
#' @return matrix; imputation results incorporating network information
#'
ImputeNetParallel <- function(drop.mat, arranged,
    cores = BiocParallel::bpworkers(BPPARAM), type = "iteration", max.iter = 50,
    BPPARAM = BiocParallel::SnowParam(type = "SOCK")) {

    d <- intersect(rownames(drop.mat)[rowSums(drop.mat) != 0],
        rownames(arranged$network))  # dropouts in >= 1 cell
    p <- colnames(arranged$network)[Matrix::colSums(
        arranged$network[d, ]) != 0]  # all available predictors

    arranged$network <- arranged$network[d, p]

    if (type == "iteration") {
        message("Starting network iterative imputation\n")
        i <- 1
        repeat {
            if ((max.iter != -1) & (i > max.iter)) {
                break
            }
            if ((i%%5) == 0) {
                message(paste("Iteration", i, "/", max.iter, "\n"))
            }

            new <- round(arranged$network %*% arranged$centered[p, ], 2)
            # expression = network coefficients * predictor expr.

            # Check convergence
            if (any(new[drop.mat[d, ]] !=
                arranged$centered[d, ][drop.mat[d, ]])) {
                arranged$centered[d, ][drop.mat[d, ]] <- new[drop.mat[d, ]]
            } else {
                break
            }
            i <- i + 1
        }
        imp <- arranged$centered
    } else {
        arranged$network <- as.matrix(arranged$network)
        imp <- BiocParallel::bplapply(seq_len(ncol(arranged$centered)),
            function(i) { PseudoInverseSolution_percell(arranged$centered[, i],
                arranged$network, drop.mat[, i]) }, BPPARAM = BPPARAM)
        imp <- do.call(cbind, imp)
        colnames(imp) <- colnames(arranged$centered)
    }
    message("Network imputation complete\n")
    return(imp)
}


#' @title Helper function to PseudoInverseSolution_percell
#'
#' @usage ImputeNPDropouts(net, expr)
#'
#' @description \code{ImputeNPDropouts} computes the non-dropout-
#' dependent solution of network imputation for each cell
#'
#' @param net matrix, logical; network coefficients for all dropout (to be
#' imputed) genes that are predictive of the expression of other dropout genes
#' @param expr numeric; vector of gene expression for all genes in the cell at
#' hand
#'
#' @return vector; imputation results for the non-dropout-dependent genes
#'
ImputeNPDropouts <- function(net, expr) {

    if (!all(colnames(net) %in% names(expr)))
        stop("Not all predictors are included in the expression vector.\n")

    net <- net[, which(Matrix::colSums(net != 0) != 0)]
    solution <- net %*% expr[colnames(net)]

    return(solution)
}


#' @title Helper function to PseudoInverseSolution_percell
#'
#' @usage ImputePredictiveDropouts(net, thr = 0.01, expr)
#'
#' @description \code{ImputePredictiveDropouts} applies Moore-Penrose
#' pseudo-inversion to compute the dropout-dependent solution of network
#' imputation for each cell
#'
#' @param net matrix, logical; network coefficients for all dropout (to be
#' imputed) genes that are predictive of the expression of other dropout genes
#' @param thr numeric; tolerance threshold to detect zero singular values
#' @param expr numeric; vector of gene expression for all genes in the cell at
#' hand
#'
#' @return vector; imputation results for the dropout-dependent genes
#'
ImputePredictiveDropouts <- function(net, thr = 0.01, expr) {

    if (!all(colnames(net) %in% names(expr)))
        stop("Not all predictors are quantified in the expression vector.\n")

    # Y = inv(I - squared_A).C

    # dropouts that are predictors and targets
    squared_A <- net[rownames(net), rownames(net), drop = FALSE]

    # I - squared_A can be inverted
    message("Computing pseudoinverse\n")
    pinv <- MASS::ginv(diag(nrow(squared_A)) - as.matrix(squared_A), tol = thr)
    dimnames(pinv) <- dimnames(squared_A)
    rm(squared_A)
    gc()

    # find C (quantified predictors)
    message("Computing constant contribution C\n")
    net <- net[, !(colnames(net) %in% rownames(net)), drop = FALSE]
    net <- net[, Matrix::colSums(net != 0) != 0, drop = FALSE]
    C <- net %*% expr[colnames(net)]
    rm(net)
    gc()

    # Y = pinv.C
    solution <- pinv %*% C

    return(solution)
}


#' @title Network-based parallel imputation - Moore-Penrose pseudoinversion
#'
#' @usage PseudoInverseSolution_percell(expr, net, drop_ind, thr = 0.01)
#'
#' @description \code{PseudoInverseSolution_percell} applies Moore-Penrose
#' pseudo-inversion to compute the solution of network imputation for each
#' cell
#'
#' @param expr numeric; expression vector for cell at hand
#' @param net matrix; network coefficients
#' @param drop_ind logical; dropout entries in the cell at hand
#' @param thr numeric; tolerance threshold to detect zero singular values
#'
#' @return matrix; imputation results incorporating network information
#'
PseudoInverseSolution_percell <- function(expr, net, drop_ind, thr = 0.01) {

    message("Starting cell pseudo-inversion\n")

    # restrict network rows to the dropouts in the data
    net <- net[intersect(rownames(net), names(expr)[drop_ind]),
        intersect(colnames(net), names(expr))]

    # restrict to predictors that are predictive of the dropouts
    net <- net[, which(Matrix::colSums(net != 0) != 0)]

    # dropouts can themselves be predictors of other genes: Y_drop =
    # A_drop*Y_drop + A_quant*Y_quant the first term depends only on quantified
    # gene expression and known network coefficients, thus is constant; the
    # second term depends on itself as the dropout genes can themselves be
    # predictive of other genes

    # Dropout-based contribution (variable)
    pd <- intersect(rownames(net), colnames(net))
    pd_imputed <- ImputePredictiveDropouts(net[pd, , drop = FALSE], thr, expr)

    # Non-dropout-based contribution
    message("Adding remaining genes\n")
    npd <- names(expr)[drop_ind][!(names(expr)[drop_ind] %in%
        rownames(pd_imputed))]
    if (length(npd) > 0) {
        npd_imputed <- ImputeNPDropouts(net[intersect(npd, rownames(net)),
            , drop = FALSE], expr)
        full_solution <- rbind(pd_imputed, npd_imputed)  # join solutions
    } else {
        full_solution <- pd_imputed
    }

    final <- expr
    final[rownames(full_solution)] <- full_solution

    return(final)
}
