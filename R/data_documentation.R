#' @title Small dataset for example purposes
#'
#' @description A small dataset to use on vignettes and examples (50 cells).
#'
#' @format matrix; a subset of the Grun pancreas dataset, obtained with the
#' \code{scRNAseq} R package, to use in the vignette and examples.
#'
#' @references Grun D et al. (2016). De novo prediction of stem cell identity
#' using single-cell transcriptome data. Cell Stem Cell 19(2), 266-277.
#'
"demo_data"

#' @title Small regulatory network for example purposes
#'
#' @description Subset of the Gene Regulatory Network used by ADImpute's Network
#' imputation method.
#'
#' @format matrix; subset of the Gene Regulatory Network installed along with
#' ADImpute.
#'
"demo_net"

#' @title Transcriptome wide gene regulatory network
#'
#' @description Gene Regulatory Network used by ADImpute's Network imputation
#' method. First column, \code{O}, corresponds to the intercept of a gene-
#' -specific predicion model. The remaining rows and columns correspond to the
#' adjacency matrix of the inferred network, where rows are target genes and
#' columns are predictors. Genes are identified by their hgnc_symbol.
#'
#' @format dgCMatrix
#'
"network.coefficients"

#' @title Table for transcript length calculations
#'
#' @description A data.frame to be used for transcript length computations.
#' May be necessary upon TPM normalization, or as input to \code{scImpute}.
#' All data was retrieved from \code{biomaRt}.
#'
#' @format A data.frame with 2 columns:
#' \describe{
#' \item{hgnc_symbol}{Gene symbol identifier}
#' \item{transcript length}{Length of transcript}
#' }
#'
"transcript_length"
