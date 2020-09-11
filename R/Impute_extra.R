#' @title Use scImpute
#'
#' @usage ImputeScImpute(count_path, infile, outfile = "rds", out_dir,
#' labeled, drop_thre, Kcluster, labels = NULL, ncores = 4, type = "TPM",
#' transcript.length = NULL, genelen = NULL)
#'
#' @description \code{ImputeScImpute} uses the scImpute package for dropout
#' imputation
#'
#' @param count_path A character specifying the full path of the raw count
#' matrix
#' @param infile A character specifying the type of file storing the raw count
#' matrix; can be "csv", "txt", or "rds". The input file shoule have rows
#' representing genes and columns representing cells, with its first row as cell
#' names and first column as gene names
#' @param outfile A character specifying the type of file storing the imputed
#' count matrix; can be "csv", "txt", or "rds"
#' @param out_dir A character specifying the full path of the output directory,
#' which is used to store all intermdediate and final outputs
#' @param labeled A logical value indicating whether cell type information is
#' available. \code{labels} must be specified if \code{labeled = TRUE}
#' @param drop_thre A number between 0 and 1, specifying the threshold to
#' determine dropout values
#' @param Kcluster An integer specifying the number of cell subpopulations.
#' This parameter can be determined based on prior knowledge or clustering of
#' raw data. Kcluster is used to determine the candidate neighbors of each cell
#' @param labels A character vector specifying the cell type of each column in
#' the raw count matrix. Only needed when \code{labeled = TRUE}. Each cell type
#' should have at least two cells for imputation
#' @param ncores A integer specifying the number of cores used for parallel
#' computation
#' @param type A character specifying the type of values in the expression
#' matrix. Can be "count" or "TPM"
#' @param transcript.length matrix with at least 2 columns: "hgnc_symbol" and
#' "transcript_length"
#' @param genelen An integer vector giving the length of each gene. Order must
#' match the gene orders in the expression matrix. genelen must be specified if
#' \code{type = "count"}
#'
#' @return matrix; imputation results from scImpute
#'
#' @seealso \code{\link[scImpute]{scimpute}}
#'
ImputeScImpute <- function(count_path,
                           infile,
                           outfile = "rds",
                           out_dir,
                           labeled,
                           drop_thre,
                           Kcluster,
                           labels = NULL,
                           ncores = 4,
                           type = "TPM",
                           transcript.length = NULL,
                           genelen = NULL){



  # Get genlen if needed
  if (type == "TPM"){

    if (is.null(transcript.length)){

      tr_length <- ADImpute::transcript_length

    } else{
      tr_length <- transcript.length
    }

    # Median length of all transcripts for a given gene
    med_length <- stats::aggregate(x = tr_length$transcript_length,
                  by = list("hgnc_symbol" = tr_length$hgnc_symbol),
                                   FUN = stats::median)
    if (infile == "txt"){
      data <- utils::read.table(count_path)
    } else if (infile == "rds"){
      data <- readRDS(count_path)
    }else{
      data <- utils::read.csv(count_path)
    }
    common <- intersect(rownames(data), med_length$hgnc_symbol)
    data   <- data[common, ]
    count_path <- paste0(strsplit(count_path,
                                  split = paste0("\\.", infile))[[1]][1],
                         "_red.csv")
    infile <- "csv"
    WriteCSV(data, count_path)
    genelen <- as.integer(med_length[match(common, med_length$hgnc_symbol), 2])
    saveRDS(genelen, paste0(out_dir,"genelength.rds"))
  }

  # Call scImpute
  scImpute::scimpute(count_path = count_path,
                     infile = infile,
                     outfile = outfile,
                     out_dir = out_dir,
                     labeled = labeled,
                     drop_thre = drop_thre,
                     Kcluster = Kcluster,
                     labels = labels,
                     ncores = ncores,
                     type = type,
                     genelen = genelen)

  # Read scImpute output
  res <- as.matrix(readRDS(paste0(out_dir, "scimpute_count.rds")))

  return(res)
}



#' @title Use SCRABBLE
#'
#' @usage ImputeSCRABBLE(data, bulk = NULL, write.to.file = T)
#'
#' @description \code{ImputeSCRABBLE} uses the SCRABBLE package for dropout
#' imputation
#'
#' @param data matrix with entries equal to zero to be imputed, normalized
#' (genes as rows and samples as columns)
#' @param bulk vector of reference bulk RNA-seq, if available (average across
#' samples)
#' @param write.to.file logical; should a file with the imputation results be
#' written?
#'
#' @return matrix; imputation results from SCRABBLE
#'
#' @seealso \code{\link[SCRABBLE]{scrabble}}
#'
ImputeSCRABBLE <- function(data, bulk = NULL, write.to.file = TRUE){

  if(is.null(bulk)){
    cat("Taking average of single cell data as reference bulk for
        SCRABBLE imputation.\n")
    bulk <- rowMeans(data)

    res <- SCRABBLE::scrabble(list(data, bulk), parameter = c(1,1e-6,1e-4))
    rownames(res) <- rownames(data)
    colnames(res) <- colnames(data)

  } else{
    # Match the rownames and order
    common <- intersect(rownames(data),names(bulk))
    if(length(common) == 0)
      stop("No common genes between single cell and bulk data.", call.=FALSE)
    data <- data[common,]
    bulk <- bulk[common]

    res <- SCRABBLE::scrabble(list(data, bulk), parameter = c(1,1e-6,1e-4))
    rownames(res) <- common
    colnames(res) <- colnames(data)
  }

  if (write.to.file){
    dir.create("SCRABBLE")
    saveRDS(res, "SCRABBLE/SCRABBLE_imputed.rds")
  }

  return(res)
}
