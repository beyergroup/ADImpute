# Helper functions to ADImpute

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


DataCheck_TranscriptLength <- function(trlength){

  # coerce to data.frame
  trlength <- as.data.frame(trlength)

  # check if required colnames are present
  if(!(c("hgnc_symbol","transcript_length") %in% colnames(trlength)))
    stop(cat("Transcript length data must contain the following colnames:",
             "hgnc_symbol, transcript_length\n"))

  # check if columns have the right format
  if(!(is.factor(trlength$hgnc_symbol) | is.character(tr_length$hgnc_symbol)))
    stop("hgnc_symbol column must be character/factor.\n")

  storage.mode(trlength$transcript_length) <- "numeric"

  return(NULL)
}
