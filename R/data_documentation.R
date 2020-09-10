#' @title Dataset for example purposes
#'
#' @description A dataset to use on vignettes and examples (100 cells).
#'
#' @format A matrix of 100 random cells sampled from a publicly
#' available dataset on hESC differention by Chu et al.
#'
#' @references Chu, L., Leng, N., Zhang, J. et al. Single-cell RNA-seq reveals
#' novel regulators of human embryonic stem cell differentiation to definitive
#' endoderm. Genome Biol 17, 173 (2016).
#' \href{https://doi.org/10.1186/s13059-016-1033-x}{DOI}
#'
#' @source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75748}{GEO accession GSE75748}
#'
"demo_data"

#' @title Small dataset for example purposes
#'
#' @description A small dataset to use on vignettes and examples (50 cells).
#'
#' @format A matrix of 50 random cells sampled from a publicly
#' available dataset on hESC differention by Chu et al.
#'
#' @references Chu, L., Leng, N., Zhang, J. et al. Single-cell RNA-seq reveals
#' novel regulators of human embryonic stem cell differentiation to definitive
#' endoderm. Genome Biol 17, 173 (2016).
#' \href{https://doi.org/10.1186/s13059-016-1033-x}{DOI}
#'
#' @source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75748}{GEO accession GSE75748}
#'
"demo_data_50cells"

#' @title Table for transcript length calculations
#'
#' @description A data.frame to be used for transcript length computations.
#' May be necessary upon TPM normalization, or as input to \code{scImpute}.
#' All data was retrieved from \code{biomaRt}.
#'
#' @format A data.frame with 3 columns:
#' \describe{
#' \item{hgnc_symbol}{Gene symbol identifier}
#' \item{ensembl_transcript_id}{ENSEMBL transcript ID}
#' \item{transcript length}{Length of transcript}
#' }
#'
"transcript_length"
