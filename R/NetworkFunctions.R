#' @title Data trimming
#'
#' @usage \code{ArrangeData(data, network.path = NULL)}
#'
#' @description \code{ArrangeData} finds common genes to the network and
#' provided data and limits both datasets to these
#'
#' @param data matrix with entries equal to zero to be imputed (genes as rows
#' and samples as columns)
#' @param network.coefficients matrix; coefficients of the gene regulatory
#' network (first column as intercept)
#'
#' @return list; data matrix, network coefficients matrix and intercept for
#' genes common between the data matrix and the network
#'
ArrangeData <- function(data,
                        network.coefficients = NULL){

  load(network.coefficients)

  O <- network.coefficients[,1]
  network_matrix <- network.coefficients[,-1]

  comm_targ <- intersect(rownames(network_matrix), rownames(data))
  comm_pred <- intersect(colnames(network_matrix), rownames(data))

  comm <- union(comm_targ, comm_pred)

  network_matrix <- round(network_matrix[comm_targ, comm_pred], 2)
  O              <- round(O[comm_targ], 2)
  data           <- data[comm, ]

  return(list("data"    = data,
              "network" = network_matrix,
              "O"       = O))
}


#' @title Data centering
#'
#' @usage \code{CenterData(data)}
#'
#' @description \code{CenterData} centers expression of each gene at 0
#'
#' @param data matrix of gene expression to be centered row-wise (genes as rows
#' and samples as columns)
#'
#' @return list; row-wise centers and centered data
#'
CenterData <- function(data, drop.exclude = F){

  cat("Centering expression of each gene at 0\n")

  if (drop.exclude){
    center <- apply(data, 1, function(x) mean(x[x != 0], na.rm = T))
    center[is.na(center)] <- 0
    center <- round(center, 2)

  } else{
    center <- round(apply(data, 1, function(x) mean(x, na.rm = T)), 2)
  }

  names(center) <- rownames(data)

  data <- data - center

  return(list("center" = center,
              "data"   = data))
}



#' @title Network-based parallel imputation
#'
#' @usage \code{ImputeNetParallel(dropout.matrix, arranged, cores = 4,
#' cluster.type = "SOCK", max.iter = 20)}
#'
#' @description \code{ImputeNetParallel} implements network-based imputation
#' in parallel
#'
#' @param dropout.matrix matrix, logical; dropout entries in the data matrix
#' (genes as rows and samples as columns)
#' @param arranged list; output of \code{\link{ArrangeData}}
#' @param cores integer; number of cores used for paralell computation
#' @param cluster.type character; either "SOCK" or "MPI"
#' @param max.iter numeric; maximum number of iterations for network
#' imputation. Set to -1 to remove limit (not recommended)
#'
#' @return matrix; imputation results incorporating network information
#'
ImputeNetParallel <- function(dropout.matrix,
                              arranged,
                              cores = 4,
                              cluster.type = "SOCK",
                              max.iter = 50){

  cluster <- snow::makeCluster(cores, type = cluster.type)

  dropouts <- intersect(rownames(dropout.matrix)[rowSums(dropout.matrix) != 0],
                        rownames(arranged$network)) # dropouts in at least one cell

  predictors <- colnames(arranged$network)[
    colSums(arranged$network[dropouts, ]) != 0] # all available predictors of the dropouts

  imp     <- arranged$centered
  network <- arranged$network[dropouts, predictors]
  O       <- arranged$O[dropouts]

  i <- 1
  repeat{

    # Limit number of iterations
    if (max.iter != -1){
      if (i > max.iter)

        break
    }

    # Print progress
    if ( (i %% 10) == 0)
      cat("Iteration", i, "\n")

    new <- round(O + snow::parMM(cluster, network, imp[predictors, ]), 2) # expression = intercept + (network coefficients * predictor expr.)

    # Check convergence
    if (any(new[dropout.matrix[dropouts, ]] != imp[dropouts, ][
      dropout.matrix[dropouts, ]])){

      imp[dropouts, ][dropout.matrix[dropouts, ]] <-
        new[dropout.matrix[dropouts, ]]
    }
    else{
      break
    }

    i <- i + 1
  }

  snow::stopCluster(cluster)

  cat("Network imputation complete\n")

  return(imp)
}
