# Helper functions to ADImpute

#' @title Argument check
#'
#' @usage CreateArgCheck(missing = NULL, match = NULL, acceptable = NULL,
#' null = NULL)
#'
#' @description \code{CreateArgCheck} creates tests for argument correctness.
#'
#' @param missing named list; logical. Name corresponds to variable name, and
#' corresponding entry to whether it was missing from the function call.
#' @param match named list. Name corresponds to variable name, and corresponding
#' entry to its value.
#' @param acceptable named list. Name corresponds to variable name, and
#' corresponding entry to its acceptable values.
#' @param null named list; logical. Name corresponds to variable name, and
#' corresponding entry to whether it was NULL in the function call.
#'
#' @return argument check object.
#'
CreateArgCheck <- function(missing = NULL, match = NULL, acceptable = NULL,
    null = NULL) {

    Check <- ArgumentCheck::newArgCheck()

    # errors for missing arguments
    if (!is.null(missing)) {
        for (varname in names(missing)) {
            if (missing[[varname]]) {
                ArgumentCheck::addError(paste("A value for ", varname,
                    " was not provided", sep = "'"), Check)
            }
        }
    }

    # errors for arguments outside of predefined options
    if (!(is.null(match)) & !(is.null(acceptable))) {
        for (varname in names(match)) {
            if (!(match[[varname]] %in% acceptable[[varname]]))
                ArgumentCheck::addError(paste(NULL, varname, " must be one of ",
                    paste(acceptable[[varname]], collapse = "', '"), NULL,
                    sep = "'"), Check)
        }
    }

    # errors for NULL arguments
    if (!is.null(null)) {
        for (varname in names(null)) {
            if (null[[varname]]) {
                ArgumentCheck::addError(paste(NULL, varname,
                    " must have non-NULL value", sep = "'"), Check)
            }
        }
    }

    return(Check)
}


CheckArguments_Impute <- function(data, method.choice, do, tr.length, labels,
    cell.clusters, true.zero.thr, drop_thre) {

    if (is.null(tr.length))
        tr.length <- ADImpute::transcript_length

    if (is.null(method.choice) & ("ensemble" %in% tolower(do)))
        stop(paste0("Please provide method.choice for Ensemble imputation. ",
            "Consider running EvaluateMethods()\n"))

    if (is.null(labels) & is.null(cell.clusters))
        stop(paste0("Please provide cell type labels ('labels') or number of",
            " cell clusters ('cell.clusters')\n"))

    if (is.null(do))
        stop("Please provide appropriate imputation methods\n")
    l <- tolower(do) %in% c("baseline", "drimpute", "ensemble", "network",
                            "saver", "scimpute", "scrabble")
    if (any(!l))
        warning(paste0("The following methods were detected as input but are",
            " not supported and will be ignored: ", paste(NULL, paste(do[!l],
                collapse = "', '"), NULL, sep = "'")))
    if ("scimpute" %in% tolower(do))
        warning(paste0("You are trying to run scImpute, which is not supported",
            " by default. Make sure you have installed scImpute from GitHub ",
            "and added the call to ImputescImpute() to the ADImpute code. ",
            "To add this call, move the code at Impute_extra.R#145 to ",
            "Wrap.R#270 and uncomment it.\n"))
    if ("scrabble" %in% tolower(do))
        warning(paste0("You are trying to run SCRABBLE, which is not supported",
            " by default. Make sure you have installed SCRABBLE from GitHub ",
            "and added the call to ImputeSCRABBLE() to the ADImpute code. ",
            "To add this call, move the code at Impute_extra.R#157 to ",
            "Wrap.R#270 and uncomment it.\n"))
    if (sum(l) == 0)
        stop("Please provide at least one supported method\n")

    if (!is.null(true.zero.thr)) {
        if ((true.zero.thr < 0) | (true.zero.thr > 1))
            stop("true.zero.thr must be a numeric between 0 and 1")
    }
    if (!is.null(drop_thre)) {
        if ((drop_thre < 0) | (drop_thre > 1))
            stop("drop_thre must be a numeric between 0 and 1")
    }

    return(NULL)
}


#' @title Data check (matrix)
#'
#' @usage DataCheck_Matrix(data)
#'
#' @description \code{DataCheck_Matrix} tests for potential format and storage
#' issues with matrices. Helper function to ADImpute.
#'
#' @param data data object to check
#'
#' @return data object with needed adjustments
#'
DataCheck_Matrix <- function(data) {

    # check format
    if (!is.matrix(data)) {
        cat("Converting input to matrix.\n")
        data <- as.matrix(data)
    }

    # check if entries are numeric
    if (!is.numeric(data)) {
        stop("Input must be numeric.\n")
    }

    # check dimnames
    if (any(is.null(dimnames(data)))) {
        stop("Input has NULL dimnames.\n")
    }

    # look for NAs
    if (any(is.na(data))) {
        cat("Converting NAs to zero.\n")
        if (all(is.na(data))) {
            stop("Input has only NAs.\n")
        }
        data[is.na(data)] <- 0
    }

    return(data)
}


#' @title Data check (network)
#'
#' @usage DataCheck_Network(network)
#'
#' @description \code{DataCheck_Network} tests for potential format and storage
#' issues with the network coefficient matrix. Helper function to ADImpute.
#'
#' @param network data object containing matrix coefficients
#'
#' @return network data object with needed adjustments
#'
DataCheck_Network <- function(network) {

    # 1st column is O
    if (colnames(network)[1] != "O")
        stop("First column of network matrix must be the intercept (O).\n")

    # is a matrix
    if (!methods::is(network, "sparseMatrix")) {
        cat("Converting network to dgCMatrix.\n")
        network <- methods::as(network, "dgCMatrix")
    }

    # contents are not all 0
    if (!(any(network != 0)))
        stop("No non-zero coefficients in network matrix.\n")

    return(network)
}


#' @title Data check (transcript length)
#'
#' @usage DataCheck_TrLength(trlength)
#'
#' @description \code{DataCheck_TrLength} tests for potential format and
#' storage issues with the object encoding transcript length, for e.g. TPM
#' normalization. Helper function to ADImpute.
#'
#' @param trlength data object containing transcript length information
#'
#' @return transcript length object with needed adjustments
#'
DataCheck_TrLength <- function(trlength) {

    # coerce to data.frame
    if (!is.data.frame(trlength)) {
        cat("Converting input to data.frame.\n")
        trlength <- as.data.frame(trlength)
    }

    if (nrow(trlength) < 1)
        stop("Not enough rows in transcript length data.\n")

    # check if required colnames are present
    if (!all(c("hgnc_symbol", "transcript_length") %in% colnames(trlength)))
        stop(cat("Transcript length data must contain the following colnames:",
            "hgnc_symbol, transcript_length\n"))

    # check if columns have the right format
    if (!(is.factor(trlength$hgnc_symbol) | is.character(trlength$hgnc_symbol)))
        stop("hgnc_symbol column must be character/factor.\n")

    storage.mode(trlength$transcript_length) <- "numeric"

    # remove lines with '' hgnc_symbol
    hgnc_symbol <- NULL  # get rid of check note
    trlength <- subset(trlength, subset = hgnc_symbol != "")
    if (nrow(trlength) < 1)
        stop("Not enough non-empty gene symbols in transcript length data.\n")

    return(trlength)
}
