#' @title TPM normalization
#'
#' @description \code{NormalizeTPM} performs TPM normalization, with possibility
#' to log the result
#'
#' @usage \code{NormalizeTPM(data, transcript.length = NULL,
#' log = F, scale = 1, pseudo.count = 1)}
#'
#' @param data matrix; raw data (genes as rows and samples as columns)
#' @param transcript.length matrix with at least 2 columns: "hgnc_symbol" and
#' "transcript_length"
#' @param log logical; log TPMs?
#' @param scale integer; scale factor to divide TPMs by
#' @param pseudocount numeric; if \code{log = T}, value to add to TPMs in order
#' to avoid taking \code{log(0)}
#'
#' @details Gene length is estimated as the median of the lengths of all
#' transcripts for each gene, as obtained from biomaRt.
#' Genes for which length information cannot be found in biomaRt are dropped.
#'
#' @return matrix; normalized data (for transcript length and library size)
#'
#'
NormalizeTPM <- function(data,
                         transcript.length = NULL,
                         log = F,
                         scale = 1,
                         pseudo.count = 1){

  if (is.null(transcript.length)){

    data("transcript_length", package = "ADImpute")
    tr_length <- transcript_length
    rm(transcript_length)

  } else {
    tr_length <- transcript.length
  }

  # Median length of all transcripts for a given gene
  med_length <- aggregate(x = tr_length$transcript_length,
                          by = list("hgnc_symbol" = tr_length$hgnc_symbol),
                          FUN = median)
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
