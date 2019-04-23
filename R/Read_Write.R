#' @title Data read
#'
#' @description \code{ReadData} reads data from raw input file (.txt or .csv)
#'
#' @usage \code{ReadData(path, ...)}
#'
#' @param path character; path to input file
#'
#' @return data.frame; raw counts (genes as rows and samples as columns)
#'
ReadData <- function(path, ...){

  # Read data
  if (grepl(path, pattern = ".txt")){

    data <- data.frame(data.table::fread(file = path,
                                         header = T,
                                         sep = "\t",
                                         showProgress = T,
                                         ...),
                       row.names = 1)

    cat("Data loaded\n")

  } else if (grepl(path, pattern = ".csv")){

    data <- data.frame(data.table::fread(file = path,
                                         header = T,
                                         sep = ",",
                                         showProgress = T,
                                         ...),
                       row.names = 1)

    cat("Data loaded\n")

  } else {

    stop("Please input txt or csv file")
  }

  # Check for NAs - if yes convert to 0
  if (sum(is.na(data)) > 0){

    data[is.na(data)] <- 0
    cat("NAs converted to 0\n")
  }

  return(as.matrix(data))
}




WriteTXT <- function(object, file){

  write.table(object, file, quote = F, sep = "\t")

  return(NULL)
}



WriteCSV <- function(object, file){

  write.csv(object, file, sep = ",")

  return(NULL)
}
