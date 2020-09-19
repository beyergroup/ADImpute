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


#' @title Data read
#'
#' @description \code{ReadData} reads data from raw input file (.txt or .csv)
#'
#' @usage ReadData(path, ...)
#'
#' @param path character; path to input file
#' @param ... additional arguments to \code{data.table::fread()}
#'
#' @return matrix; raw counts (genes as rows and samples as columns)
#'
ReadData <- function(path, ...){

  # Read data
  if (grepl(path, pattern = ".txt")){

    data <- data.frame(data.table::fread(file = path,
                                         header = TRUE,
                                         sep = "\t",
                                         showProgress = TRUE,
                                         ...),
                       row.names = 1)

    cat("Data loaded\n")

  } else if (grepl(path, pattern = ".csv")){

    data <- data.frame(data.table::fread(file = path,
                                         header = TRUE,
                                         sep = ",",
                                         showProgress = TRUE,
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



#' @title Write txt file
#'
#' @description \code{WriteTXT} writes data to a tab-delimited output file
#'
#' @usage WriteTXT(object, file)
#'
#' @param object R object to write
#' @param file character; path to output file
#'
#' @return Returns NULL
#'
#' @examples
#' file <- tempfile()
#' WriteTXT(iris, file = file)
#'
#' @export
#'
WriteTXT <- function(object, file){

  utils::write.table(object, file, quote = FALSE, sep = "\t")

  return(NULL)
}



#' @title Write csv file
#'
#' @description \code{WriteCSV} writes data to a comma-delimited output file
#'
#' @usage WriteCSV(object, file)
#'
#' @param object R object to write
#' @param file character; path to output file
#'
#' @return Returns NULL
#'
#' @examples
#' file <- tempfile()
#' WriteCSV(iris, file = file)
#'
#' @export
#'
WriteCSV <- function(object, file){

  utils::write.csv(object, file)

  return(NULL)
}

