# Helper functions to ADImpute


DataCheck_Arranged <- function(arranged){


  return()
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
DataCheck_Matrix <- function(data){

  # check format
  if(!is.matrix(data)){
    cat("Converting input to matrix.\n")
    data <- as.matrix(data)
  }

  # check if entries are numeric
  if(!is.numeric(data)){stop("Input must be numeric.\n")}

  # check dimnames
  if(any(is.null(dimnames(data)))){stop("Input has NULL dimnames.\n")}

  # look for NAs
  if(any(is.na(data))){
    cat("Converting NAs to zero.\n")
    if(all(is.na(data))){stop("Input has only NAs.\n")}
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
DataCheck_Network <- function(network){

  # 1st column is O
  if(colnames(network)[1] != "O")
    stop("First column of network matrix must be the intercept (O).\n")

  # is a matrix
  if(!is.matrix(network)){
    cat("Converting network to matrix.\n")
    network <- as.matrix(network)
  }

  # contents are numeric
  storage.mode(network) <- "numeric"

  # contents are not all 0
  if(!(any(network != 0)))
    stop("No non-zero coefficients in network matrix.\n")

  return(network)
}


#' @title Data check (transcript length)
#'
#' @usage DataCheck_TranscriptLength(trlength)
#'
#' @description \code{DataCheck_TranscriptLength} tests for potential format and
#' storage issues with the object encoding transcript length, for e.g. TPM
#' normalization. Helper function to ADImpute.
#'
#' @param trlength data object containing transcript length information
#'
#' @return transcript length object with needed adjustments
#'
DataCheck_TranscriptLength <- function(trlength){

  # coerce to data.frame
  if(!is.data.frame(trlength)){
    cat("Converting input to data.frame.\n")
    trlength <- as.data.frame(trlength)
  }

  if(nrow(trlength) < 1)
    stop("Not enough rows in transcript length data.\n")

  # check if required colnames are present
  if(!all(c("hgnc_symbol","transcript_length") %in% colnames(trlength)))
    stop(cat("Transcript length data must contain the following colnames:",
             "hgnc_symbol, transcript_length\n"))

  # check if columns have the right format
  if(!(is.factor(trlength$hgnc_symbol) | is.character(trlength$hgnc_symbol)))
    stop("hgnc_symbol column must be character/factor.\n")

  storage.mode(trlength$transcript_length) <- "numeric"

  # remove lines with "" hgnc_symbol
  trlength <- subset(trlength, hgnc_symbol != "")
  if(nrow(trlength) < 1)
    stop("Not enough non-empty gene symbols in transcript length data.\n")

  return(trlength)
}
