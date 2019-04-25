#' @title Imputation method evaluation on training set
#'
#' @usage \code{EvaluateMethods(data, training.ratio = .7, training.only = T,
#' mask.ratio = .1, split.seed = NULL, mask.seed = NULL, scale = 1,
#' pseudo.count = 1, labels = NULL, cell.clusters = NULL, drop_thre = NULL,
#' type = "TPM", cores = 4, cluster.type = "SOCK", network.path = NULL,
#' transcript.length = NULL, drop.exclude = T, ...)}
#'
#' @description \code{EvaluateMethods} returns the best-performing imputation
#' method for each gene in the dataset
#'
#' @param data matrix; normalized counts, not logged (genes as rows and samples
#' as columns)
#' @param training.ratio numeric; ratio of samples to be used for training
#' @param training.only logical; if TRUE define only a training dataset, if
#' FALSE writes and returns both training and validation sets (defaults to T)
#' @param mask.ratio numeric; ratio of samples to be masked per gene
#' @param split.seed integer; optional seed (defaults to NULL, random selection
#' of samples to use for training)
#' @param mask.seed integer; optional seed (defaults to NULL, random selection
#' of samples to mask)
#' @param scale integer; scaling factor to divide all expression levels by
#' (defaults to 1)
#' @param pseudo.count integer; pseudo-count to be added to expression levels
#' to avoid log(0) (defaults to 1)
#' @param labels character; vector specifying the cell type of each column of
#' \code{data}
#' @param cell.clusters integer; number of cell subpopulations
#' @param drop_thre numeric; between 0 and 1 specifying the threshold to
#' determine dropout values
#' @param type A character specifying the type of values in the expression
#' matrix. Can be "count" or "TPM"
#' @param cores integer; number of cores used for paralell computation
#' @param cluster.type character; either "SOCK" or "MPI"
#' @param network.path character; path to .txt or .rds file with network
#' coefficients
#' @param transcript.length matrix with at least 2 columns: "hgnc_symbol" and
#' "transcript_length"
#' @param drop.exclude logical; should zeros be discarded for the calculation
#' of genewise average expression levels? (defaults to T)
#'
#' @details For each gene, a fraction (\code{mask.ratio}) of the quantified
#' expression values are set to zero and imputed according to 3 different
#' methods: scImpute, baseline (average gene expression across all cells) or a
#' network-based method. The imputation error is computed for each of the
#' values in the original dataset that was set to 0, for each method. The
#' method resulting in a lowest imputation error for each gene is chosen.
#'
#' @return character; best performing method in the training set for each gene
#'
#' @seealso \code{\link{ImputeBaseline}},
#' \code{\link{ImputeNetwork}},
#' \code{\link{ImputeScImpute}},
#' \code{\link{scimpute}}
#'

EvaluateMethods <- function(data,
                            training.ratio = .7,
                            training.only = T,
                            mask.ratio = .1,
                            split.seed = NULL,
                            mask.seed = NULL,
                            scale = 1,
                            pseudo.count = 1,
                            labels = NULL,
                            cell.clusters = NULL,
                            drop_thre = NULL,
                            type = "TPM",
                            cores = 4,
                            cluster.type = "SOCK",
                            network.path = NULL,
                            transcript.length = NULL,
                            drop.exclude = T,
                            ...){

  # Check arguments
  if (is.null(transcript.length)){
    data("transcript_length", package = "ADImpute")
    transcript.length <- transcript_length
    rm(transcript_length)
  }

  dir.create("training")
  setwd("training")
  on.exit(setwd("../"))

  # Select training data
  cat("Selecting training data\n")
  training_norm <- SplitData(data,
                             ratio = training.ratio,
                             training.only = training.only,
                             seed = split.seed)

  # Mask selected training data
  cat("Masking training data\n")
  masked_training_norm <- MaskData(training_norm,
                                   write.to.file = T,
                                   filename = "masked_training_norm.txt",
                                   mask = mask.ratio,
                                   seed = mask.seed)

  # Run scImpute
  cat("Running scImpute on training data\n")
  dir.create("scImpute")

  labeled <- F # by default consider the cell types are unknown
  if (!is.null(labels))
    labeled <- T

  scimputed <- ImputeScImpute(count_path = "masked_training_norm.txt",
                              infile = "txt",
                              out_dir = "scImpute/",
                              labeled = labeled,
                              Kcluster = cell.clusters,
                              labels = labels,
                              drop_thre = drop_thre,
                              ncores = cores,
                              type = type,
                              transcript.length = transcript.length)

  # Log masked data
  log_masked_training_norm <- log2( (masked_training_norm / scale) +
                                      pseudo.count)
  WriteTXT(log_masked_training_norm, "log_masked_training_norm.txt")

  # Run Baseline
  cat("Imputing training data using average expression\n")
  baseline <- ImputeBaseline(log_masked_training_norm,
                             drop.exclude = drop.exclude,
                             write.to.file = T)
  baseline_norm <- round(scale * ( (2 ^ baseline) - pseudo.count), 2)
  WriteTXT(baseline_norm, "Baseline/baseline_imputed_norm.txt")

  # Run Net
  cat("Imputing training data using network information\n")
  net <- ImputeNetwork(log_masked_training_norm,
                       network.path,
                       cores,
                       cluster.type,
                       drop.exclude = drop.exclude,
                       write.to.file = T,
                       ...)
  net_norm <- round(scale * ( (2 ^ net) - pseudo.count), 2)
  WriteTXT(net_norm, "Network/network_imputed_norm.txt")

  # Run optimum choice
  choice <- ChooseMethod(real = round(training_norm, 2),
                         masked = round(masked_training_norm, 2),
                         scimpute = scimputed,
                         baseline = baseline_norm,
                         network = net_norm)

  return(choice)
}


#' @title Dropout imputation using gene-specific best-performing methods
#'
#' @usage \code{Impute(data, method.choice, scale = 1, pseudo.count = 1,
#' count_path, labels = NULL, cell.clusters = NULL, drop_thre = NULL,
#' type = "TPM", cores = 4, cluster.type = "SOCK", network.path = NULL,
#' transcript.length = NULL, drop.exclude = T, ...)}
#'
#' @description \code{Impute} performs dropout imputation based on the
#' performance results obtained in the training data, coupled to normalization
#' using \code{normalization.function}
#'
#' @param data matrix; raw counts (genes as rows and samples as columns)
#' @param method.choice character; best performing method in training data for
#' each gene
#' @param scale integer; scaling factor to divide all expression levels by
#' (defaults to 1)
#' @param pseudo.count integer; pseudo-count to be added to expression levels
#' to avoid log(0) (defaults to 1)
#' @param count_path character; path to data file
#' @param labels character; vector specifying the cell type of each column of
#' \code{data}
#' @param cell.clusters integer; number of cell subpopulations
#' @param drop_thre numeric; between 0 and 1 specifying the threshold to
#' determine dropout values
#' @param type A character specifying the type of values in the expression
#' matrix. Can be "count" or "TPM"
#' @param cores integer; number of cores used for paralell computation
#' @param cluster.type character; either "SOCK" or "MPI"
#' @param network.path character; path to .txt or .rds file with network
#' coefficients
#' @param transcript.length matrix with at least 2 columns: "hgnc_symbol" and
#' "transcript_length"
#' @param drop.exclude logical; should zeros be discarded for the calculation
#' of genewise average expression levels? (defaults to T)
#'
#' @return matrix; imputed and normalized expression values
#'
#' @details Values that are 0 in \code{data} are imputed according to the
#' best-performing methods indicated in \code{methods}.
#'
#' @seealso \code{\link{EvaluateMethods}},
#' \code{\link{ImputeBaseline}},
#' \code{\link{ImputeNetwork}},
#' \code{\link{ImputeScImpute}},
#' \code{\link{scimpute}}
#'
Impute <- function(data,
                   method.choice,
                   scale = 1,
                   pseudo.count = 1,
                   count_path,
                   labels = NULL,
                   cell.clusters = NULL,
                   drop_thre = NULL,
                   type = "TPM",
                   cores = 4,
                   cluster.type = "SOCK",
                   network.path = NULL,
                   transcript.length = NULL,
                   drop.exclude = T,
                   ...){

  # Check arguments
  if (is.null(transcript.length)){
    data("transcript_length", package = "ADImpute")
    transcript.length <- transcript_length
    rm(transcript_length)
  }

  dir.create("imputation")
  setwd("imputation")
  on.exit(setwd("../"))

  # Run scImpute
  cat("Running scImpute\n")
  dir.create("scImpute")
  labeled = F # by default consider the cell types are unknown
  if (!is.null(labels)){
    labeled = T
  }

  scimputed <- ImputeScImpute(count_path,
                              infile = "txt",
                              outfile = "rds",
                              out_dir = "scImpute/",
                              labeled = labeled,
                              Kcluster = cell.clusters,
                              labels = labels,
                              drop_thre = drop_thre,
                              ncores = 1,
                              type = type,
                              transcript.length = transcript.length)

  # Log masked data
  log_masked_norm <- log2( (data / scale) + pseudo.count)
  WriteTXT(log_masked_norm, "log_masked_norm.txt")

  # Run Baseline
  cat("Imputing data using average expression\n")
  baseline <- ImputeBaseline(log_masked_norm,
                             drop.exclude = drop.exclude,
                             write.to.file = T)
  baseline_norm <- round(scale * ( (2 ^ baseline) - pseudo.count), 2)
  WriteTXT(baseline_norm, "Baseline/baseline_imputed_norm.txt")

  # Run Net
  cat("Imputing data using network information\n")
  net <- ImputeNetwork(log_masked_norm,
                       network.path,
                       cores,
                       cluster.type,
                       drop.exclude = drop.exclude,
                       write.to.file = T,
                       ...)
  net_norm <- round(scale * ( (2 ^ net) - pseudo.count), 2)
  WriteTXT(net_norm, "Network/network_imputed_norm.txt")

  # Combine results
  imputed <- Combine(round(log2( (data / scale) + pseudo.count), 2),
                     method.choice,
                     round(log2( (scimputed / scale) + pseudo.count), 2),
                     round(baseline, 2),
                     round(net, 2))

  return(imputed)
}
