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


#' @title Imputation method evaluation on training set
#'
#' @usage EvaluateMethods(data, do = c("Baseline", "DrImpute", "Network"),
#' training.ratio = .7, training.only = T, mask.ratio = .1, split.seed = NULL,
#' mask.seed = NULL, scale = 1, pseudo.count = 1, labels = NULL,
#' cell.clusters = 2, drop_thre = NULL, type = "TPM", cores = 4,
#' cluster.type = "SOCK", network.path = NULL, transcript.length = NULL,
#' drop.exclude = T, bulk = NULL, ...)
#'
#' @description \code{EvaluateMethods} returns the best-performing imputation
#' method for each gene in the dataset
#'
#' @param data matrix; normalized counts, not logged (genes as rows and samples
#' as columns)
#' @param do character; choice of methods to be used for imputation. Currently
#' supported methods are \code{"Baseline"}, \code{"DrImpute"} and
#' \code{"Network"}. Not case-sensitive. Can include one or more methods.
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
#' @param bulk vector of reference bulk RNA-seq, if available (average across
#' samples)
#' @param ... additional parameters to pass to network-based imputation
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
#' \code{\link{ImputeDrImpute}},
#' \code{\link{ImputeNetwork}}
#'
#' @export
#'
EvaluateMethods <- function(data,
                            do = c("Baseline", "DrImpute", "Network"),
                            training.ratio = .7,
                            training.only = T,
                            mask.ratio = .1,
                            split.seed = NULL,
                            mask.seed = NULL,
                            scale = 1,
                            pseudo.count = 1,
                            labels = NULL,
                            cell.clusters = 2,
                            drop_thre = NULL,
                            type = "TPM",
                            cores = 4,
                            cluster.type = "SOCK",
                            network.path = NULL,
                            transcript.length = NULL,
                            drop.exclude = T,
                            bulk = NULL,
                            ...){

  # Check arguments
  if (is.null(transcript.length)){
    transcript.length <- ADImpute::transcript_length
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

  train_imputed <- list()

  # Run scImpute
  if("scimpute" %in% tolower(do)){
    cat("scImpute is not supported by default. Please check XXXX for details.\n")
    # cat("Imputing training data using scImpute\n")
    # dir.create("scImpute")
    # labeled = F # by default consider the cell types are unknown
    # if (!is.null(labels))
    #   labeled <- T
    # train_imputed$scImpute <- ImputeScImpute(count_path = "masked_training_norm.txt",
    #                                          infile = "txt",
    #                                          outfile = "rds",
    #                                          out_dir = "scImpute/",
    #                                          labeled = labeled,
    #                                          Kcluster = cell.clusters,
    #                                          labels = labels,
    #                                          drop_thre = drop_thre,
    #                                          ncores = cores,
    #                                          type = type,
    #                                          transcript.length = transcript.length)
    # train_imputed$scImpute <- log2( (train_imputed$scImpute / scale) + pseudo.count)
  }

  # Run SCRABBLE
  if("scrabble" %in% tolower(do)){
    cat("SCRABBLE is not supported by default. Please check XXXX for details.\n")
    # cat("Imputing training data using SCRABBLE\n")
    # dir.create("SCRABBLE")
    # train_imputed$SCRABBLE <- ImputeSCRABBLE(masked_training_norm, bulk)
    # train_imputed$SCRABBLE <- log2( (train_imputed$SCRABBLE / scale) + pseudo.count)
  }

  # Log masked data
  log_masked_training_norm <- log2( (masked_training_norm / scale) +
                                      pseudo.count)
  WriteTXT(log_masked_training_norm, "log_masked_training_norm.txt")

  # Run Baseline
  if("baseline" %in% tolower(do)){
    cat("Imputing training data using average expression\n")
    train_imputed$Baseline <- ImputeBaseline(log_masked_training_norm,
                                             drop.exclude = drop.exclude,
                                             write.to.file = T)
    baseline_norm <- round(scale * ( (2 ^ train_imputed$Baseline) - pseudo.count), 2)
    WriteTXT(baseline_norm, "Baseline/baseline_imputed_norm.txt")
  }

  # Run Net
  if("network" %in% tolower(do)){
    cat("Imputing training data using network information\n")
    train_imputed$Network <- ImputeNetwork(log_masked_training_norm,
                                           network.path,
                                           cores,
                                           cluster.type,
                                           drop.exclude = drop.exclude,
                                           write.to.file = T,
                                           ...)
    net_norm <- round(scale * ( (2 ^ train_imputed$Network) - pseudo.count), 2)
    WriteTXT(net_norm, "Network/network_imputed_norm.txt")
  }

  # Run DrImpute
  if("drimpute" %in% tolower(do)){
    cat("Imputing training data using DrImpute\n")
    train_imputed$DrImpute <- ImputeDrImpute(log_masked_training_norm)
  }

  # Run optimum choice
  choice <- ChooseMethod(real = round(training_norm, 2),
                         masked = round(masked_training_norm, 2),
                         imputed = train_imputed)

  return(choice)
}


#' @title Dropout imputation using gene-specific best-performing methods
#'
#' @usage Impute(data, do = "Ensemble", method.choice = NULL, scale = 1,
#' pseudo.count = 1, count_path, labels = NULL, cell.clusters = 2,
#' drop_thre = NULL, type = "TPM", cores = 4, cluster.type = "SOCK",
#' network.path = NULL, transcript.length = NULL, drop.exclude = T,
#' bulk = NULL, true.zero.thr = NULL, ...)
#'
#' @description \code{Impute} performs dropout imputation based on the
#' performance results obtained in the training data, coupled to normalization
#' using \code{normalization.function}
#'
#' @param data matrix; raw counts (genes as rows and samples as columns)
#' @param do character; choice of methods to be used for imputation. Currently
#' supported methods are \code{"Baseline"}, \code{"DrImpute"},
#' \code{"Network"}, and \code{"Ensemble"}.
#' Not case-sensitive. Can include one or more methods.
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
#' @param bulk vector of reference bulk RNA-seq, if available (average across
#' samples)
#' @param true.zero.thr if set to NULL (default), no true zero estimation is
#' performed. Set to numeric value between 0 and 1 for estimation. Value
#' corresponds to the threshold used to determine true zeros: if the probability
#' of dropout is lower than \code{true.zero.thr}, the imputed entries are set
#' to zero.
#' @param ... additional parameters to pass to network-based imputation
#'
#' @return list of imputation results (normalized, log-transformed) for all
#' selected methods in \code{do}
#'
#' @details Values that are 0 in \code{data} are imputed according to the
#' best-performing methods indicated in \code{method.choice}. Currently
#' supported methods are:
#' \itemize{
#'  \item \code{Baseline}: imputation with average expression across all cells
#'  in the dataset. See \code{\link{ImputeBaseline}}.
#'  \item Previously published approaches: \code{DrImpute}, \code{scImpute} and
#'  \code{SCRABBLE}.
#'  \item \code{Network}: leverages information from a gene regulatory network
#'  to predicted expression of genes that are not quantified based on
#'  quantified interacting genes, in the same cell. See
#'  \code{\link{ImputeNetwork}}.
#'  \item \code{Ensemble}: is based on results on a training subset of the data at
#'  hand, indicating which method best predicts the expression of each gene.
#'  These results are supplied via \code{method.choice}. Applies the imputation
#'  results of the best performing method to the zero entries of each gene.
#' }
#' If \code{"Ensemble"} is included in \code{do}, \code{method.choice} has to
#' be provided (use output from \code{EvaluateMethods()}).
#' \code{Impute} creates a directory \code{imputation} containing the
#' imputation results of all methods in \code{do}.
#'
#' @examples
#' # Normalize demo data
#' norm_data <- NormalizeRPM(demo_data_50cells[1:100,])
#' # Impute with particular methods
#' imputed_data <- Impute(do = c("Baseline","DrImpute"), data = norm_data,
#' cores = 1)
#'
#' @seealso \code{\link{EvaluateMethods}},
#' \code{\link{ImputeBaseline}},
#' \code{\link{ImputeDrImpute}},
#' \code{\link{ImputeNetwork}},
#' \code{\link{ImputeSAVER}}
#'
#' @export
#'
Impute <- function(data,
                   do = "Ensemble",
                   method.choice = NULL,
                   scale = 1,
                   pseudo.count = 1,
                   count_path = NULL,
                   labels = NULL,
                   cell.clusters = 2,
                   drop_thre = NULL,
                   type = "TPM",
                   cores = 4,
                   cluster.type = "SOCK",
                   network.path = NULL,
                   transcript.length = NULL,
                   drop.exclude = T,
                   bulk = NULL,
                   true.zero.thr = NULL,
                   ...){

  imputed <- list()

  # Check arguments
  if (is.null(transcript.length)){
    transcript.length <- ADImpute::transcript_length
  }
  if(is.null(method.choice) & ("ensemble" %in% tolower(do))){
    stop("Please provide method.choice for Ensemble imputation. Consider running EvaluateMethods()\n")
  }

  dir.create("imputation")
  setwd("imputation")
  on.exit(setwd("../"))

  if("ensemble" %in% tolower(do)){
    do <- union(do, unique(method.choice))
  }

  # Run scImpute
  if("scimpute" %in% tolower(do)){
    cat("scImpute is not supported by default. Please check XXXX for details.\n")
    # cat("Imputing data using scImpute\n")
    # dir.create("scImpute")
    # labeled = F # by default consider the cell types are unknown
    # if (!is.null(labels)){
    #   labeled = T
    # }
    #
    # # Get type of input masked file
    # if(any(grep(".txt", count_path))){
    #   infile <- "txt"
    # } else if(any(grep(".rds", count_path))){
    #   infile <- "rds"
    # } else if (any(grep(".csv", count_path))){
    #   infile <- "csv"
    # } else{
    #   stop("Cannot recognize input count path to scImpute. Please provide .txt, .rds or .csv file\n")
    # }
    #
    # imputed$scImpute <- ImputeScImpute(count_path,
    #                                    infile = infile,
    #                                    outfile = "rds",
    #                                    out_dir = "scImpute/",
    #                                    labeled = labeled,
    #                                    Kcluster = cell.clusters,
    #                                    labels = labels,
    #                                    drop_thre = drop_thre,
    #                                    ncores = cores,
    #                                    type = type,
    #                                    transcript.length = transcript.length)
    # imputed$scImpute <- log2( (imputed$scImpute / scale) + pseudo.count)
  }

  # Run SAVER
  if("saver" %in% tolower(do)){
    cat("Imputing data using SAVER\n")
    imputed$SAVER <- ImputeSAVER(data, cores, try.mean = F)
    imputed$SAVER <- log2( (imputed$SAVER / scale) + pseudo.count)
    saver_norm <- round(scale * ( (2 ^ imputed$SAVER) - pseudo.count), 2)
    WriteTXT(saver_norm, "SAVER/SAVER_imputed_norm.txt")
  }

  # Run SCRABBLE
  if("scrabble" %in% tolower(do)){
    cat("scImpute is not supported by default. Please check XXXX for details.\n")
    # cat("Imputing data using SCRABBLE\n")
    # imputed$SCRABBLE <- ImputeSCRABBLE(data, bulk)
    # imputed$SCRABBLE <- log2( (imputed$SCRABBLE / scale) + pseudo.count)
  }

  # Log masked data
  log_masked_norm <- log2( (data / scale) + pseudo.count)
  WriteTXT(log_masked_norm, "log_masked_norm.txt")

  # Run Baseline
  if("baseline" %in% tolower(do)){
    cat("Imputing data using average expression\n")
    imputed$Baseline <- ImputeBaseline(log_masked_norm,
                                       drop.exclude = drop.exclude,
                                       write.to.file = T)
    baseline_norm <- round(scale * ( (2 ^ imputed$Baseline) - pseudo.count), 2)
    WriteTXT(baseline_norm, "Baseline/baseline_imputed_norm.txt")
  }


  # Run Net
  if("network" %in% tolower(do)){
    cat("Imputing data using network information\n")
    imputed$Network <- ImputeNetwork(log_masked_norm,
                                     network.path,
                                     cores,
                                     cluster.type,
                                     drop.exclude = drop.exclude,
                                     write.to.file = T,
                                     ...)
    net_norm <- round(scale * ( (2 ^ imputed$Network) - pseudo.count), 2)
    WriteTXT(net_norm, "Network/network_imputed_norm.txt")
  }


  # Run DrImpute
  if("drimpute" %in% tolower(do)){
    cat("Imputing data using DrImpute\n")
    imputed$DrImpute <- ImputeDrImpute(log_masked_norm)
  }

  # Combine results
  if("ensemble" %in% tolower(do)){
    cat("Combining imputation results into ensemble\n")
    imputed$Ensemble <- Combine(log_masked_norm, imputed, method.choice)
  }

  # Estimate true zeros

  if(any(grep(".txt", count_path))){
    infile <- "txt"
  } else if(any(grep(".rds", count_path))){
    infile <- "rds"
  } else if (any(grep(".csv", count_path))){
    infile <- "csv"
  }

  if(!is.null(true.zero.thr)){

    # error to run scImpute first

    # Get genelen if needed
    if (type == "TPM"){
      count_path <- paste0(strsplit(count_path,
                                    split = paste0("\\.", infile))[[1]][1],
                           "_red",
                           paste0(".", infile))
      infile <- "csv"
      # compute dropout probabilities according to scImpute
      droprob <- GetDropoutProbabilities(infile = infile, count_path = count_path,
                                         out_dir = "scImpute/", type = type,
                                         genelen = readRDS(paste0("outdir","genelength.rds")),
                                         drop_thre = true.zero.thr, data = utils::read.csv(count_path))
      WriteTXT(droprob, "dropout_probability.txt")

      # apply thresholds
      zerofiltered <- lapply(imputed, SetBiologicalZeros, drop_probs = droprob, thre = true.zero.thr,
                             was_zero = data == 0)

    } else {
      # compute dropout probabilities according to scImpute
      droprob <- GetDropoutProbabilities(infile = infile, count_path = count_path,
                                         out_dir = "scImpute/", type = type, genelen = NULL,
                                         drop_thre = true.zero.thr, data = data)
      WriteTXT(droprob, "dropout_probability.txt")

      # apply thresholds
      zerofiltered <- lapply(imputed, SetBiologicalZeros, drop_probs = droprob, thre = true.zero.thr,
                             was_zero = data == 0)
    }

    return(list("imputations" = imputed, "zerofiltered" = zerofiltered))

  } else{

    return(imputed)
  }

}
