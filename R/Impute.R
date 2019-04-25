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

    } else{

      data <- read.csv(count_path)

    }
    common <- intersect(rownames(data), med_length$hgnc_symbol)
    data   <- data[common, ]
    count_path <- paste0(strsplit(count_path,
                                  split = paste0("\\.", infile))[[1]][1],
                         "_red",
                         paste0(".", infile))
    infile <- "csv"
    WriteCSV(data, count_path)
    genelen <- as.integer(med_length[match(common, med_length$hgnc_symbol), 2])

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
  out <- as.matrix(readRDS(paste0(out_dir, "scimpute_count.rds")))

  return(out)
}


#' @title Impute using average expression across all cells
#'
#' @usage ImputeBaseline(data, write.to.file = T, drop.exclude = T, ...)
#'
#' @description \code{ImputeBaseline} imputes dropouts using gene averages
#' across cells
#'
#' @param data matrix with entries equal to zero to be imputed (genes as rows
#' and samples as columns)
#' @param write.to.file logical; should a file with the imputation results be
#' written?
#' @param drop.exclude logical; should zeros be discarded for the calculation
#' of genewise average expression levels? (defaults to T)
#' @return matrix; imputation results considering the average expression
#' values of genes
#'
ImputeBaseline <- function(data, write.to.file = T, drop.exclude = T, ...){

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
  baseline_imputed <- data
  baseline_imputed[dropouts] <- gene_avg_mat[dropouts]

  if (write.to.file){

    dir.create("Baseline")
    WriteTXT(baseline_imputed, "Baseline/baseline_imputed.txt", ...)
  }

  return(baseline_imputed)
}



#' @title Network-based imputation
#'
#' @usage \code{ImputeNetwork(data, network.path = NULL, cores = 4,
#' cluster.type = "SOCK", write.to.file = T, drop.exclude = T, ...)}
#'
#' @param data matrix with entries equal to zero to be imputed, already
#' normalized (genes as rows and samples as columns)
#' @param network.path character; path to .txt or .rds file with network
#' coefficients
#' @param cores integer; number of cores to use (parallel computation if cores
#' > 1)
#' @param write.to.file logical; should a file with the imputation results be
#' written?
#' @param drop.exclude logical; should zeros be discarded for the calculation
#' of genewise average expression levels? (defaults to T)
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

  if (cores >= 1){

    new_imp <- ImputeNetParallel(dropout_mat,
                                 arranged,
                                 cores,
                                 cluster.type,
                                 ...)
  } else{
    # non parallel - Not implemented yet
  }

  imputed <- data
  imputed[rownames(new_imp), ][dropout_mat[rownames(new_imp), ]] <-
    (centered$center[rownames(new_imp)] + new_imp)[
      dropout_mat[rownames(new_imp), ]]

  imputed[imputed < 0] <- 0

  if (write.to.file){

    dir.create("Network")
    WriteTXT(imputed, "Network/network_imputed.txt")
  }

  return(imputed)
}



#' @title Combine imputation methods
#'
#' @usage \code{Combine(data, method.choice, scimputed, baseline, net,
#' write.to.file = T)}
#'
#' @param data matrix with entries equal to zero to be imputed, already
#' normalized (genes as rows and samples as columns)
#' @param method.choice named character; vector with the best performing method
#' per gene
#' @param scimputed matrix; data imputed by scimputed
#' @param baseline matrix; data imputed by baseline method
#' @param net matrix; data imputed by network method
#' @param write.to.file logical; should a file with the imputation results be
#' written?
#'
#' @details Combines imputation results from all methods according to training
#' results provided in \code{method.choice}
#'
#' @return matrix; imputation results combining the best performing method
#' per gene
#'
#' @seealso \code{\link{ImputeBaseline}},
#' \code{\link{ImputeNetwork}},
#' \code{\link{ImputeScImpute}},
#' \code{\link{scImpute::scimpute}}
#'
Combine <- function(data,
                    method.choice,
                    scimputed,
                    baseline,
                    net,
                    write.to.file = T){

  dropouts <- data == 0

  scimpute_genes <- names(method.choice)[method.choice == "scimpute"]
  baseline_genes <- names(method.choice)[method.choice == "baseline"]
  network_genes <- rownames(data)[!(rownames(data) %in% c(scimpute_genes, baseline_genes))]

  scimpute_index <- dropouts & (matrix(rownames(data),
                                       byrow = F,
                                       nrow = nrow(data),
                                       ncol = ncol(data))
                                %in% scimpute_genes)

  baseline_index <- dropouts & (matrix(rownames(data),
                                       byrow = F,
                                       nrow = nrow(data),
                                       ncol = ncol(data))
                                %in% baseline_genes)

  network_index  <- dropouts & (matrix(rownames(data),
                                       byrow = F,
                                       nrow = nrow(data),
                                       ncol = ncol(data))
                                %in% network_genes)

  imputed <- data
  imputed[scimpute_index] <- scimputed[scimpute_index]
  imputed[baseline_index] <- baseline[baseline_index]
  imputed[network_index]  <- net[network_index]

  if (write.to.file)
    WriteTXT(imputed, "final_imputation.txt")

  return(imputed)
}



ImputeSAVER <- function(data, cores, try.mean = T){

  library(SAVER)

  dir.create("SAVER")

  if(try.mean){
    imp_mean <- saver(data, size.factor = 1, ncores = cores, null.model = T)
    write.table(imp_mean$estimate, "SAVER/SAVER_nullmodel_imputed.txt", quote = F, sep = "\t")
    saveRDS(object = imp_mean, file = "SAVER/SAVER_nullmodel.rds")
  }

  imp <- saver(data, ncores = cores, size.factor = 1)
  write.table(imp$estimate, "SAVER/SAVER_imputed.txt", quote = F, sep = "\t")
  saveRDS(object = imp, file = "SAVER/SAVER.rds")

  return(imp$estimate)
}

