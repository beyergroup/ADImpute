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

    coll <- checkmate::makeAssertCollection()

    # errors for missing arguments
    if (!is.null(missing)){
        for(varname in names(missing)){
            if (missing[[varname]]){
                coll$push(paste0("A value for ", varname,
                    " was not provided", sep = "'"))
            }
        }
    }

    # errors for arguments outside of predefined options
    if (!is.null(match) & !is.null(acceptable)){
        for (varname in names(match)) {
            checkmate::assert_subset(x = match[[varname]],
                choices = acceptable[[varname]],
                .var.name = varname,
                add = coll)
        }
    }

    # error for NULL arguments
    if (!is.null(null)){
        for (varname in names(null)) {
            if (null[[varname]]){
                coll$push(paste(NULL, varname, " must have a non-NULL value",
                                sep = "'"))
            }
        }
    }

    return(coll)
}


#' @title Argument check to Impute()
#'
#' @usage CheckArguments_Impute(data, method.choice, do, tr.length, labels,
#' cell.clusters, true.zero.thr, drop_thre)
#'
#' @param data matrix; raw counts (genes as rows and samples as columns)
#' @param method.choice character; best performing method in training data for
#' each gene
#' @param do character; choice of methods to be used for imputation. Currently
#' supported methods are \code{'Baseline'}, \code{'DrImpute'},
#' \code{'Network'}, and \code{'Ensemble'}. Defaults to \code{'Ensemble'}.
#' Not case-sensitive. Can include one or more methods. Non-supported methods
#' will be ignored.
#' @param tr.length matrix with at least 2 columns: 'hgnc_symbol' and
#' 'transcript_length'
#' @param labels character; vector specifying the cell type of each column of
#' \code{data}
#' @param cell.clusters integer; number of cell subpopulations
#' @param true.zero.thr if set to NULL (default), no true zero estimation is
#' performed. Set to numeric value between 0 and 1 for estimation. Value
#' corresponds to the threshold used to determine true zeros: if the probability
#' of dropout is lower than \code{true.zero.thr}, the imputed entries are set
#' to zero.
#' @param drop_thre numeric; between 0 and 1 specifying the threshold to
#' determine dropout values
#'
#' @description \code{CheckArguments_Impute} checks whether the arguments passed
#' to \code{Impute} are correct.
#'
#' @return NULL object
#'
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

    # wrongly passed a SingleCellExperiment object
    if(class(data)[1] == "SingleCellExperiment"){
        stop(paste("SingleCellExperiment object passed as matrix. Please use",
            "argument 'sce' instead of 'data'.\n"))
    }

    # check format
    if (!is.matrix(data)) {
        message("Converting input to matrix.\n")
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
        message("Converting NAs to zero.\n")
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
        message("Converting network to dgCMatrix.\n")
        network <- methods::as(network, "dgCMatrix")
    }

    # contents are not all 0
    if (!(any(network != 0)))
        stop("No non-zero coefficients in network matrix.\n")

    return(network)
}


#' @title Data check (SingleCellExperiment)
#'
#' @usage DataCheck_SingleCellExperiment(sce, normalized = TRUE)
#'
#' @description \code{DataCheck_SingleCellExperiment} tests for existence of the
#' appropriate assays in \code{sce}. Helper function to ADImpute.
#'
#' @param sce SingleCellExperiment; data for normalization or imputation
#' @param normalized logical; is the data expected to be normalized?
#'
#' @return NULL object.
#'
DataCheck_SingleCellExperiment <- function(sce, normalized = TRUE){

    if(normalized & (!("normcounts" %in%
        names(SummarizedExperiment::assays(sce))))){
        stop(paste("Normalized counts not found in 'normcounts' assay of",
            "SingleCellExperiment object.\n"))
    } else if((!normalized) & (!("counts" %in%
        names(SummarizedExperiment::assays(sce))))){
        stop(paste("Raw counts not found in 'counts' assay of",
            "SingleCellExperiment object.\n"))
    }

    return(NULL)
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
        message("Converting input to data.frame.\n")
        trlength <- as.data.frame(trlength)
    }

    if (nrow(trlength) < 1)
        stop("Not enough rows in transcript length data.\n")

    # check if required colnames are present
    if (!all(c("hgnc_symbol", "transcript_length") %in% colnames(trlength)))
        stop(paste("Transcript length data must contain the following",
            "colnames: hgnc_symbol, transcript_length\n"))

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


#' @title Wrapper for return of EvaluateMethods()
#'
#' @usage ReturnChoice(sce, choice)
#'
#' @description \code{ReturnChoice} Adjusts the output of \code{EvaluateMethods}
#' to a character vector or a SingleCellExperiment object. Helper function to
#' ADImpute.
#'
#' @param sce SingleCellExperiment; a SingleCellExperiment object if available;
#' NULL otherwise
#' @param choice character; best performing method in the training set for each
#' gene
#'
#' @return
#'  \itemize{
#'     \item if \code{sce} is provided: returns a SingleCellExperiment with the
#'     best performing method per gene stored as row-features. Access via
#'     \code{SingleCellExperiment::int_elementMetadata(sce)$ADImpute$methods}.
#'     \item if \code{sce} is not provided: returns a character with the best
#'     performing method in the training set for each gene
#' }
#'
ReturnChoice <- function(sce, choice){
    if(is.null(sce)){
        message("Returning method choice vector.\n")
        return(choice)
    } else{
        message("Returning SingleCellExperiment object.\n")
        v <- rep(NA, nrow(sce))
        names(v) <- rownames(sce)
        v[names(choice)] <- choice
        SingleCellExperiment::int_elementMetadata(sce)$ADImpute <-
            S4Vectors::DataFrame(method = v)
        return(sce)
    }
}


#' @title Wrapper for return of Impute()
#'
#' @usage ReturnOut(result, sce)
#'
#' @description \code{ReturnOut} Adjusts the output of \code{Impute} to a list
#' of matrices or a SingleCellExperiment object. Helper function to ADImpute.
#'
#' @param result list; imputation result
#' @param sce SingleCellExperiment; a SingleCellExperiment object if available;
#' NULL otherwise
#'
#' @return imputation results. A SingleCellExperiment if \code{!is.null(sce)},
#' or a list with imputed results in matrix format otherwise.
#'
ReturnOut <- function(result, sce){

    if(is.null(sce)){
        message("Returning matrix lists.\n")
        return(result)
    } else{
        message("Returning SingleCellExperiment object.\n")
        if(length(result) == 0){
            return(sce)
        }
        if(!is.list(result[[1]])){
            SummarizedExperiment::assays(sce) <-
                c(SummarizedExperiment::assays(sce), result)
        } else{
            SummarizedExperiment::assays(sce) <-
                c(SummarizedExperiment::assays(sce), result$zerofiltered)
        }
        return(sce)
    }
}
