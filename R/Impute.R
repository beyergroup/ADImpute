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


#' @title Combine imputation methods
#'
#' @usage \code{Combine(data, method.choice, imputed, write.to.file = T)}
#'
#' @param data matrix with entries equal to zero to be imputed, already
#' normalized (genes as rows and samples as columns)
#' @param method.choice named character; vector with the best performing method
#' per gene
#' @param imputed list; list of matrices with imputation results for all
#' considered methods
#' @param write.to.file logical; should a file with the imputation results be
#' written?
#'
#' @details Combines imputation results from all methods according to training
#' results provided in \code{method.choice}
#'
#' @return matrix; imputation results combining the best performing method
#' per gene
#'
#'
Combine <- function(data,
                    imputed,
                    method.choice,
                    write.to.file = T){

  # all zeros are imputed
  dropouts <- data == 0

  best <- sapply(unique(method.choice),
                 function(m) names(method.choice)[method.choice == m])
  # some genes could not have a method assigned to them, because they were too
  # lowly expressed for any masking to be done. These are assigned to Network:
  best$Network <- c(best$Network, setdiff(rownames(data), unlist(best)))

  # combine imputations for best performing methods
  combined <- sapply(names(best), function(m) imputed[[m]][best[[m]],])
  combined <- do.call(rbind, combined)
  combined <- combined[rownames(data),colnames(data)]

  # replace dropouts in the data with combined imputations
  imputed_combined <- data
  imputed_combined[dropouts] <- combined[dropouts]

  if (write.to.file)
    WriteTXT(imputed_combined, "final_imputation.txt")

  return(imputed_combined)
}


#' @title Impute using average expression across all cells
#'
#' @usage ImputeBaseline(data, write.to.file = T, drop.exclude = T, ...)
#'
#' @description \code{ImputeBaseline} imputes dropouts using gene averages
#' across cells
#'
#' @param data matrix with entries equal to zero to be imputed, normalized
#' and log2-transformed (genes as rows and samples as columns)
#' @param write.to.file logical; should a file with the imputation results be
#' written?
#' @param drop.exclude logical; should zeros be discarded for the calculation
#' of genewise average expression levels? (defaults to T)
#' @param ... additional arguments to \code{saveRDS}
#'
#' @return matrix; imputation results considering the average expression
#' values of genes
#'
ImputeBaseline <- function(data,
                           write.to.file = T,
                           drop.exclude = T,
                           ...){

  dropouts <- data == 0

  # Compute mean expression levels across cells
  if (drop.exclude){
    gene_avg <- apply(data, 1, function(x) mean(x[x != 0], na.rm = T))
    gene_avg[is.na(gene_avg)] <- 0

  } else{
    gene_avg <- apply(data, 1, function(x) mean(x, na.rm = T))
  }

  # Impute when dropout
  gene_avg_mat <- matrix(data = gene_avg,
                         nrow = nrow(data),
                         ncol = ncol(data),
                         byrow = F,
                         dimnames = dimnames(data))
  res <- data
  res[dropouts] <- gene_avg_mat[dropouts]

  if (write.to.file){
    dir.create("Baseline")
    saveRDS(res, "Baseline/baseline_imputed.rds", ...)
  }

  return(res)
}


#' @title Use DrImpute
#'
#' @usage \code{ImputeDrImpute(data, write.to.file = T)}
#'
#' @description \code{ImputeDrImpute} uses the DrImpute package for dropout
#' imputation
#'
#' @param data matrix with entries equal to zero to be imputed, normalized
#' and log2-transformed (genes as rows and samples as columns)
#' @param write.to.file logical; should a file with the imputation results be
#' written?
#'
#' @return matrix; imputation results from DrImpute
#'
#' @seealso \code{\link{DrImpute::DrImpute}}
#'
ImputeDrImpute <- function(data, write.to.file = T){

  res <- DrImpute(as.matrix(data))
  colnames(res) <- colnames(data)

  if(write.to.file){
    dir.create("DrImpute")
    saveRDS(res, "DrImpute/drimputed.rds")
  }

  return(res)
}


#' @title Network-based imputation
#'
#' @usage \code{ImputeNetwork(data, network.path = NULL, cores = 4,
#' cluster.type = "SOCK", write.to.file = T, drop.exclude = T, ...)}
#'
#' @param data matrix with entries equal to zero to be imputed, normalized
#' and log2-transformed (genes as rows and samples as columns)
#' @param network.path character; path to .txt or .rds file with network
#' coefficients
#' @param cores integer; number of cores to use
#' @param cluster.type character; either "SOCK" or "MPI"
#' @param write.to.file logical; should a file with the imputation results be
#' written?
#' @param drop.exclude logical; should zeros be discarded for the calculation
#' of genewise average expression levels? (defaults to T)
#' @param ... additional arguments to \code{ImputeNetParallel}
#'
#' @details Imputes dropouts using a gene regulatory network trained on external
#' data, as provided in \code{network.path}. Dropout expression values are
#' estimated from the expression of their predictor genes and the network
#' coefficients.
#'
#' @return matrix; imputation results incorporating network information
#'
#' @seealso \code{\link{ImputeNetParallel}}
#'
ImputeNetwork <- function(data,
                          network.path = NULL,
                          cores = 4,
                          cluster.type = "SOCK",
                          write.to.file = T,
                          drop.exclude = T,
                          ...){

  # Limit data and network to genes common to both
  arranged <- ArrangeData(data, network.path)

  cat("Dimensions of data matrix:",
      dim(arranged$data)[1], "x", dim(arranged$data)[2],
      "\n")
  cat("Dimensions of network matrix:",
      dim(arranged$network)[1], "x", dim(arranged$network)[2],
      "\n")

  # Center expression of each gene
  centered <- CenterData(arranged$data, drop.exclude)
  arranged$centered <- as.matrix(centered$data)

  dropout_mat <- arranged$data == 0 # dropout indexes in the data matrix

  cat("Starting network-based imputation\n")

  new_imp <- ImputeNetParallel(dropout_mat,
                               arranged,
                               cores,
                               cluster.type,
                               ...)

  res <- data
  # add back gene means
  res[rownames(new_imp), ][dropout_mat[rownames(new_imp), ]] <-
    (centered$center[rownames(new_imp)] + new_imp)[
      dropout_mat[rownames(new_imp), ]]

  res[res < 0] <- 0

  if (write.to.file){
    dir.create("Network")
    saveRDS(res, "Network/network_imputed.rds")
  }

  return(res)
}


#' @title Use SAVER
#'
#' @usage \code{ImputeSAVER(data, cores = 4, try.mean = F,
#' write.to.file = T)}
#'
#' @description \code{ImputeSAVER} uses the SAVER package for dropout
#' imputation
#'
#' @param data matrix with entries equal to zero to be imputed, normalized
#' (genes as rows and samples as columns)
#' @param cores integer; number of cores to use
#' @param try.mean logical; whether to additionally use mean gene expression
#' as prediction
#' @param write.to.file logical; should a file with the imputation results be
#' written?
#'
#' @return matrix; imputation results from SAVER
#'
#' @seealso \code{\link{SAVER::saver}}
#'
ImputeSAVER <- function(data, cores, try.mean = F, write.to.file = T){

  dir.create("SAVER")

  if(try.mean){
    imp_mean <- SAVER::saver(data, size.factor = 1, ncores = cores, null.model = T)
    saveRDS(object = imp_mean, file = "SAVER/SAVER_nullmodel.rds")
  }

  res <- SAVER::saver(data, ncores = cores, size.factor = 1)

  if (write.to.file){
    dir.create("SAVER")
    saveRDS(object = res, file = "SAVER/SAVER.rds")
  }

  return(res$estimate)
}


#' @title Use scImpute
#'
#' @usage \code{ImputeScImpute(count_path, infile, outfile = "rds", out_dir,
#' labeled, drop_thre, Kcluster, labels = NULL, ncores = 4, type = "TPM",
#' transcript.length = NULL, genelen = NULL)}
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
#' @seealso \code{\link{scImpute::scimpute}}
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

      data("transcript_length", package = "ADImpute")
      tr_length <- transcript_length
      rm(transcript_length)

    } else{
      tr_length <- transcript.length
    }

    # Median length of all transcripts for a given gene
    med_length <- aggregate(x = tr_length$transcript_length,
                            by = list("hgnc_symbol" = tr_length$hgnc_symbol),
                            FUN = median)
    if (infile == "txt"){
      data <- read.table(count_path)
    } else if (infile == "rds"){
      data <- readRDS(count_path)
    }else{
      data <- read.csv(count_path)
    }
    common <- intersect(rownames(data), med_length$hgnc_symbol)
    data   <- data[common, ]
    count_path <- paste0(strsplit(count_path,
                                  split = paste0("\\.", infile))[[1]][1],
                         "_red.csv")
    infile <- "csv"
    WriteCSV(data, count_path)
    genelen <- as.integer(med_length[match(common, med_length$hgnc_symbol), 2])
    saveRDS(genelength, paste0(out_dir,"genelength.rds"))
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
#' @usage \code{ImputeSCRABBLE(data, bulk = NULL, write.to.file = T)}
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
#' @seealso \code{\link{SCRABBLE::scrabble}}
#'
ImputeSCRABBLE <- function(data, bulk = NULL, write.to.file = T){

  if(is.null(bulk)){
    cat("Taking average of single cell data as reference bulk for SCRABBLE imputation.\n")
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
