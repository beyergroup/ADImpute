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
#' write = FALSE, train.ratio = .7, train.only = TRUE, mask.ratio = .1,
#' outdir = getwd(), scale = 1, pseudo.count = 1, labels = NULL,
#' cell.clusters = 2, drop_thre = NULL, type = "count", cores = 4,
#' network.coefficients = NULL, network.path = system.file("extdata",
#' "network.coefficients.zip", package="ADImpute"), network.imputation =
#' "iteration", transcript.length = NULL, drop.exclude = TRUE, bulk = NULL, ...)
#'
#' @description \code{EvaluateMethods} returns the best-performing imputation
#' method for each gene in the dataset
#'
#' @param data matrix; normalized counts, not logged (genes as rows and samples
#' as columns)
#' @param do character; choice of methods to be used for imputation. Currently
#' supported methods are \code{"Baseline"}, \code{"DrImpute"} and
#' \code{"Network"}. Not case-sensitive. Can include one or more methods. Non-
#' supported methods will be ignored.
#' @param write logical; write intermediary and imputed objects to files?
#' @param train.ratio numeric; ratio of samples to be used for training
#' @param train.only logical; if TRUE define only a training dataset, if
#' FALSE writes and returns both training and validation sets (defaults to TRUE)
#' @param mask.ratio numeric; ratio of samples to be masked per gene
#' @param outdir character; path to directory where output files are written.
#' Defaults to working directory
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
#' @param network.coefficients matrix; network coefficients. Please provide
#' either \code{network.coefficients} or \code{network.path}.
#' @param network.path character; path to .txt or .rds file with network
#' coefficients. Defaults to the zip file included in ADImpute installation
#' @param network.imputation character; either "iteration", for an iterative
#' solution, or "pseudoinv", to use Moore-Penrose pseudo-inversion as a
#' solution.
#' @param transcript.length matrix with at least 2 columns: "hgnc_symbol" and
#' "transcript_length"
#' @param drop.exclude logical; should zeros be discarded for the calculation
#' of genewise average expression levels? (defaults to TRUE)
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
#' @examples
#' # Normalize demo data
#' norm_data <- NormalizeRPM(ADImpute::demo_data)
#' method_choice <- EvaluateMethods(norm_data, do = c("Baseline","DrImpute"),
#' outdir = tempdir())
#' unlink(file.path(tempdir(),"ADImpute"), recursive = TRUE, force = TRUE)
#'
#' @seealso \code{\link{ImputeBaseline}},
#' \code{\link{ImputeDrImpute}},
#' \code{\link{ImputeNetwork}}
#'
#' @export
#'
EvaluateMethods <- function(data,
                            do = c("Baseline", "DrImpute", "Network"),
                            write = FALSE, train.ratio = .7,
                            train.only = TRUE, mask.ratio = .1,
                            outdir = getwd(), scale = 1, pseudo.count = 1,
                            labels = NULL, cell.clusters = 2, drop_thre = NULL,
                            type = "count", cores = 4,
                            network.coefficients = NULL,
                            network.path = system.file("extdata",
                              "network.coefficients.zip", package="ADImpute"),
                            network.imputation = "iteration",
                            transcript.length = NULL, drop.exclude = TRUE,
                            bulk = NULL, ...){

  # Check arguments
  if (is.null(transcript.length))
    transcript.length <- ADImpute::transcript_length
  Check <- CreateArgCheck(missing = list("data" = missing(data)),
                          match = list("type" = type),
                          acceptable = list("type" = c("TPM","count")))
  ArgumentCheck::finishArgCheck(Check)
  data <- DataCheck_Matrix(data)
  transcript.length <- DataCheck_TranscriptLength(transcript.length)

  if(write){
    savedir <- getwd(); dir.create(paste0(outdir,"/training"))
    setwd(paste0(outdir,"/training")); on.exit(setwd(savedir))
  }

  train_data <- CreateTrainData(data, train.ratio = train.ratio,
                                train.only = train.only,
                                mask = mask.ratio, write = write)

  train_imputed <- Impute(data = train_data$mask, do = do[do != "Ensemble"],
                          write = write, outdir = getwd(),
                          scale = scale, pseudo.count = pseudo.count,
                          # scImpute arguments below
                          count_path = NULL, labels = labels,
                          cell.clusters = cell.clusters, drop_thre = drop_thre,
                          type = type, transcript.length = transcript.length,
                          # add count_path argument for scImpute runs
                          bulk = bulk, # SCRABBLE argument
                          cores = cores,
                          network.coefficients = network.coefficients,
                          network.path = network.path,
                          drop.exclude = drop.exclude, ...)

  # Run optimum choice
  choice <- ChooseMethod(real = round(train_data$train, 2),
                         masked = round(train_data$mask, 2),
                         imputed = train_imputed, write.to.file = write)

  return(choice)
}


#' @title Dropout imputation using gene-specific best-performing methods
#'
#' @usage Impute(data, do = "Ensemble", write = FALSE,outdir = getwd(),
#' method.choice = NULL, scale = 1, pseudo.count = 1, count_path = NULL,
#' labels = NULL, cell.clusters = 2, drop_thre = NULL, type = "count",
#' transcript.length = NULL, cores = 4, network.coefficients = NULL,
#' network.path = system.file("extdata", "network.coefficients.zip", package =
#' "ADImpute"), network.imputation = "iteration", drop.exclude = TRUE, bulk =
#' NULL, true.zero.thr = NULL, prob.mat = NULL, ...)
#'
#' @description \code{Impute} performs dropout imputation based on the
#' performance results obtained in the training data, coupled to normalization
#' using \code{normalization.function}
#'
#' @param data matrix; raw counts (genes as rows and samples as columns)
#' @param do character; choice of methods to be used for imputation. Currently
#' supported methods are \code{"Baseline"}, \code{"DrImpute"},
#' \code{"Network"}, and \code{"Ensemble"}. Defaults to \code{"Ensemble"}.
#' Not case-sensitive. Can include one or more methods. Non-supported methods
#' will be ignored.
#' @param write logical; write intermediary and imputed objects to files?
#' @param outdir character; path to directory where output files are written.
#' Defaults to working directory
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
#' @param transcript.length matrix with at least 2 columns: "hgnc_symbol" and
#' "transcript_length"
#' @param cores integer; number of cores used for paralell computation
#' @param network.coefficients matrix; network coefficients. Please provide
#' either \code{network.coefficients} or \code{network.path}.
#' @param network.path character; path to .txt or .rds file with network
#' coefficients. Defaults to the zip file included in ADImpute installation
#' @param network.imputation character; either "iteration", for an iterative
#' solution, or "pseudoinv", to use Moore-Penrose pseudo-inversion as a
#' solution.
#' @param drop.exclude logical; should zeros be discarded for the calculation
#' of genewise average expression levels? (defaults to TRUE)
#' @param bulk vector of reference bulk RNA-seq, if available (average across
#' samples)
#' @param true.zero.thr if set to NULL (default), no true zero estimation is
#' performed. Set to numeric value between 0 and 1 for estimation. Value
#' corresponds to the threshold used to determine true zeros: if the probability
#' of dropout is lower than \code{true.zero.thr}, the imputed entries are set
#' to zero.
#' @param prob.mat matrix of the same size as data, filled with the dropout
#' probabilities for each gene in each cell
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
#'  \item \code{Ensemble}: is based on results on a training subset of the data
#'  at hand, indicating which method best predicts the expression of each gene.
#'  These results are supplied via \code{method.choice}. Applies the imputation
#'  results of the best performing method to the zero entries of each gene.
#' }
#' If \code{"Ensemble"} is included in \code{do}, \code{method.choice} has to
#' be provided (use output from \code{EvaluateMethods()}).
#' \code{Impute} creates a directory \code{imputation} containing the
#' imputation results of all methods in \code{do}.
#' If \code{true.zero.thr} is set, dropout probabilities are computed using
#' scImpute's framework. Expression values with dropout probabilities below
#' \code{true.zero.thr} will be set back to 0 if imputed, as they likely
#' correspond to true biological zeros (genes not expressed in cell) rather than
#' technical dropouts (genes expressed but not captured).
#'
#' @examples
#' # Normalize demo data
#' norm_data <- NormalizeRPM(demo_data)
#' # Impute with particular method(s)
#' imputed_data <- Impute(do = "Network", data = norm_data[,1:10],
#' network.coefficients = ADImpute::demo_net, cores = 2)
#' imputed_data <- Impute(do = "Network", data = norm_data[,1:10],
#' network.imputation = "pseudoinv", network.coefficients = ADImpute::demo_net,
#' cores = 2)
#' # Don't impute biological zeros
#' imputed_data <- Impute(do = "Baseline", data = norm_data, cores = 2,
#' true.zero.thr = .2)
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
                   write = FALSE, outdir = getwd(),
                   method.choice = NULL, scale = 1, pseudo.count = 1,
                   # scImpute arguments:
                   count_path = NULL, labels = NULL, cell.clusters = 2,
                   drop_thre = NULL, type = "count", transcript.length = NULL, #
                   cores = 4, network.coefficients = NULL,
                   network.path = system.file("extdata",
                                              "network.coefficients.zip",
                                              package="ADImpute"),
                   network.imputation = "iteration",
                   drop.exclude = TRUE,
                   bulk = NULL, # SCRABBLE argument
                   true.zero.thr = NULL, prob.mat = NULL, ...){

  # Check arguments
  if (is.null(transcript.length))
    transcript.length <- ADImpute::transcript_length
  if(is.null(method.choice) & ("ensemble" %in% tolower(do))){
    stop("Please provide method.choice for Ensemble imputation.
         Consider running EvaluateMethods()\n")
  }
  Check <- CreateArgCheck(missing = list("data" = missing(data)),
                          match = list("type" = type),
                          acceptable = list("type" = c("TPM","count")))
  ArgumentCheck::finishArgCheck(Check)
  data <- DataCheck_Matrix(data)
  transcript.length <- DataCheck_TranscriptLength(transcript.length)

  if(write){
    savedir <- getwd(); dir.create(paste0(outdir,"/imputation"))
    setwd(paste0(outdir,"/imputation")); on.exit(setwd(savedir))
  }

  if("ensemble" %in% tolower(do))
    do <- union(c(do, "Network"), unique(method.choice))

  imputed <- list()

  if("scimpute" %in% tolower(do)){
    cat("scImpute is not supported by default. Please check XXXX for details.\n")
    # # Get type of input masked file
    # if(grepl(".txt", count_path)){infile <- "txt"
    # } else if(grepl(".rds", count_path)){infile <- "rds"
    # } else if (grepl(".csv", count_path)){infile <- "csv"
    # } else{stop(paste0("Cannot recognize input count path to scImpute.",
    #                    "Please provide .txt, .rds or .csv file\n"))}
    #
    # imputed$scImpute <- ImputeScImpute(count_path,
    #                                    infile = infile,
    #                                    outfile = "rds",
    #                                    out_dir = "scImpute/",
    #                                    labeled = is.null(labels),
    #                                    Kcluster = cell.clusters,
    #                                    labels = labels,
    #                                    drop_thre = drop_thre,
    #                                    ncores = cores,
    #                                    type = type,
    #                                    transcript.length = transcript.length)
    # imputed$scImpute <- log2( (imputed$scImpute / scale) + pseudo.count)
  }

  if("saver" %in% tolower(do)){
    imputed$SAVER <- ImputeSAVER(data, cores, try.mean = FALSE,
                                 write.to.file = write)
    imputed$SAVER <- log2( (imputed$SAVER / scale) + pseudo.count)
  }

  if("scrabble" %in% tolower(do)){
    cat("scImpute is not supported by default. Please check XXXX for details.\n")
    # imputed$SCRABBLE <- ImputeSCRABBLE(data, bulk)
    # imputed$SCRABBLE <- log2( (imputed$SCRABBLE / scale) + pseudo.count)
  }

  log_masked_norm <- log2( (data / scale) + pseudo.count)

  if("baseline" %in% tolower(do))
    imputed$Baseline <- ImputeBaseline(log_masked_norm,
                                       drop.exclude = drop.exclude,
                                       write.to.file = write)

  if("network" %in% tolower(do))
    imputed$Network <- ImputeNetwork(log_masked_norm,
                                     network.coefficients,
                                     network.path,
                                     type = network.imputation,
                                     cores,
                                     drop.exclude = drop.exclude,
                                     write.to.file = write,
                                     ...)

  if("drimpute" %in% tolower(do))
    imputed$DrImpute <- ImputeDrImpute(log_masked_norm, write.to.file = write)

  if("ensemble" %in% tolower(do))
    imputed$Ensemble <- Combine(log_masked_norm, imputed, method.choice, write)

  # Estimate true zeros
  if(!is.null(true.zero.thr)){

    # the probability is not given - compute it
    if(is.null(prob.mat))
      prob.mat <- GetDropoutProbabilities(data = data, thre = true.zero.thr,
                                          cell.clusters = cell.clusters,
                                          labels = labels, type = type,
                                          ncores = cores,
                                          genelen = transcript.length)
    # use probability matrix
    cl <- parallel::makeCluster(cores)
    zerofiltered <- parallel::parLapply(cl, imputed, SetBiologicalZeros,
                                        drop_probs = prob.mat,
                                        thre = true.zero.thr,
                                        was_zero = data == 0)
    parallel::stopCluster(cl)

    return(list("imputations" = imputed, "zerofiltered" = zerofiltered))

  } else{

    return(imputed)
  }
}
