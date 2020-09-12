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
