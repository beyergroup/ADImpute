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


#' @title RPM normalization
#'
#' @description \code{NormalizeRPM} performs RPM normalization, with possibility
#' to log the result
#'
#' @usage NormalizeRPM(data, log = FALSE, scale = 1, pseudo.count = 1)
#'
#' @param data matrix; raw data (genes as rows and samples as columns)
#' @param log logical; log RPMs?
#' @param scale integer; scale factor to divide RPMs by
#' @param pseudo.count numeric; if \code{log = TRUE}, value to add to RPMs in
#' order to avoid taking \code{log(0)}
#'
#' @return matrix; library size normalized data
#'
#' @examples
#' demo <- NormalizeRPM(demo_data_50cells)
#'
#' @export
#'
NormalizeRPM <- function(data,
                         log = FALSE,
                         scale = 1,
                         pseudo.count = 1){

  # divide by million RPK in sample - transcripts per million (TPM)
  data <- sweep(data, 2, STATS = colSums(data) / ( 10 ^ 6), FUN = "/")

  data <- data / scale

  if (log)
    data <- log2(data + pseudo.count)

  return(data)
}


#' @title TPM normalization
#'
#' @description \code{NormalizeTPM} performs TPM normalization, with possibility
#' to log the result
#'
#' @usage NormalizeTPM(data, transcript.length = NULL,
#' log = FALSE, scale = 1, pseudo.count = 1)
#'
#' @param data matrix; raw data (genes as rows and samples as columns)
#' @param transcript.length matrix with at least 2 columns: "hgnc_symbol" and
#' "transcript_length"
#' @param log logical; log TPMs?
#' @param scale integer; scale factor to divide TPMs by
#' @param pseudo.count numeric; if \code{log = T}, value to add to TPMs in order
#' to avoid taking \code{log(0)}
#'
#' @details Gene length is estimated as the median of the lengths of all
#' transcripts for each gene, as obtained from biomaRt.
#' Genes for which length information cannot be found in biomaRt are dropped.
#'
#' @return matrix; normalized data (for transcript length and library size)
#'
#' @examples
#' demo <- NormalizeTPM(demo_data_50cells)
#'
#' @export
#'
NormalizeTPM <- function(data,
                         transcript.length = NULL,
                         log = FALSE,
                         scale = 1,
                         pseudo.count = 1){

  if (is.null(transcript.length)){

    tr_length <- ADImpute::transcript_length

  } else {
    tr_length <- transcript.length
  }

  # Median length of all transcripts for a given gene
  med_length <- stats::aggregate(x = tr_length$transcript_length,
                                 by = list("hgnc_symbol" = tr_length$hgnc_symbol),
                                 FUN = stats::median)
  common <- intersect(rownames(data), med_length$hgnc_symbol)
  data   <- data[common, ]

  med_length <- med_length[match(common, med_length$hgnc_symbol), ]

  # divide by length of transcript in kb - reads per kilobase (RPK)
  data <- sweep(data, 1, STATS = med_length$x / 1000, FUN = "/")
  # divide by million RPK in sample - transcripts per million (TPM)
  data <- sweep(data, 2, STATS = colSums(data) / ( 10 ^ 6), FUN = "/")

  data <- data / scale

  if (log)
    data <- log2(data + pseudo.count)

  return(data)
}
