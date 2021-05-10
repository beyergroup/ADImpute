#' @title Use scImpute
#'
#' @usage ImputeScImpute(data, labeled, drop_thre, Kcluster, labels = NULL,
#' cores = 4, type = 'count', outdir = tempdir(), tr.length =
#' ADImpute::transcript_length, genelen = NULL)
#'
#' @description \code{ImputeScImpute} uses the scImpute package for dropout
#' imputation
#'
#' @param data matrix; data to be imputed
#' @param outdir A character specifying the full path of the output directory,
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
#' @param cores A integer specifying the number of cores used for parallel
#' computation
#' @param type A character specifying the type of values in the expression
#' matrix. Can be 'count' or 'TPM'
#' @param tr.length matrix with at least 2 columns: 'hgnc_symbol' and
#' 'transcript_length'
#' @param genelen An integer vector giving the length of each gene. Order must
#' match the gene orders in the expression matrix. genelen must be specified if
#' \code{type = 'count'}
#'
#' @return matrix; imputation results from scImpute
#'
#' @seealso \code{\link[scImpute]{scimpute}}
#'
ImputeScImpute <- function(data, labeled, drop_thre, Kcluster, labels = NULL,
    cores = 4, type = "count", outdir = tempdir(),
    tr.length = ADImpute::transcript_length, genelen = NULL) {

    # Get genlen if needed
    if (type == "TPM") {
        if (is.null(tr.length)) {
            tr_length <- ADImpute::transcript_length
        } else {
            tr_length <- tr.length
        }
        common <- intersect(rownames(data), tr_length$hgnc_symbol)
        data <- data[common, ]
        genelen <- tr_length[match(common, tr_length$hgnc_symbol),
                             "transcript_length"]
    }

    if (labeled == TRUE & is.null(labels))
        stop("'labels' must be specified when 'labeled = TRUE'!")
    if (labeled == FALSE & is.null(Kcluster))
        stop("'Kcluster' must be specified when 'labeled = FALSE'!")
    if (!(type %in% c("count", "TPM")))
        stop("expression values can be either 'count' or 'TPM'!")
    if (type == "TPM" & is.null(genelen))
        stop("'genelen' must be specified when type = 'TPM'!")

    count_lnorm = read_count(raw_count = data, type = type, genelen = genelen)
    saveRDS(count_lnorm, paste0(outdir, "/count_lnorm.rds"))

    if (labeled == TRUE) {
        if (length(labels) != ncol(count_lnorm))
            stop("number of cells does not match number of labels !")
    }

    # Call scImpute's functions
    scImpute::scimpute(count_path = paste0(outdir, "/count_lnorm.rds"),
        infile = "rds", outfile = "rds", out_dir = outdir, labeled = labeled,
        drop_thre = drop_thre, Kcluster = Kcluster, labels = labels,
        cores = cores, type = type, genelen = genelen)

    # Read scImpute output
    res <- as.matrix(readRDS(paste0(outdir, "scimpute_count.rds")))

    tryCatch(unlink(outdir, recursive = T))

    return(res)
}



#' @title Use SCRABBLE
#'
#' @usage ImputeSCRABBLE(data, bulk = NULL, write = FALSE)
#'
#' @description \code{ImputeSCRABBLE} uses the SCRABBLE package for dropout
#' imputation
#'
#' @param data matrix with entries equal to zero to be imputed, normalized
#' (genes as rows and samples as columns)
#' @param bulk vector of reference bulk RNA-seq, if available (average across
#' samples)
#' @param write logical; should a file with the imputation results be written?
#'
#' @return matrix; imputation results from SCRABBLE
#'
#' @seealso \code{\link[SCRABBLE]{scrabble}}
#'
ImputeSCRABBLE <- function(data, bulk = NULL, write = FALSE) {

    if (is.null(bulk)) {
        message("Taking average of single cell data as reference bulk for
        SCRABBLE imputation.\n")
        bulk <- rowMeans(data)

        res <- tryCatch(SCRABBLE::scrabble(list(data, bulk),
            parameter = c(1, 1e-06, 1e-04)), error = function(e) {
                stop(paste("Error:", e$message, "\nIs SCRABBLE installed?"))
            })
        rownames(res) <- rownames(data)
        colnames(res) <- colnames(data)

    } else {
        # Match the rownames and order
        common <- intersect(rownames(data), names(bulk))
        if (length(common) == 0)
            stop("No common genes between single cell and bulk data.",
                call. = FALSE)
        data <- data[common, ]
        bulk <- bulk[common]

        res <- SCRABBLE::scrabble(list(data, bulk), parameter = c(1, 1e-06,
            1e-04))
        rownames(res) <- common
        colnames(res) <- colnames(data)
    }

    if (write) {
        dir.create("SCRABBLE")
        saveRDS(res, "SCRABBLE/SCRABBLE_imputed.rds")
    }

    return(res)
}

