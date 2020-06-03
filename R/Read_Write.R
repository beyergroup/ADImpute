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
#' @usage \code{ReadData(path, ...)}
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

