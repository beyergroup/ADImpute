#' @title Network loading
#'
#' @usage \code{ReadNetwork(network.path)}
#'
#' @description \code{ReadNetwork} loads the matrix of network coefficients
#'
#' @param network.path character; path to .txt or .rds file with network
#' coefficients
#'
#' @return matrix with network coefficients and intercept in the first column
#'
ReadNetwork <- function(network.path){

  if (grepl(network.path, pattern = ".txt")){

    cat("Reading .txt file with network coefficients\n")

    m <- readLines(network.path)
    cnames <- m[1]
    cnames <- strsplit(cnames, split = "\t")[[1]]
    mm <- sapply(m[-1], function(x) strsplit(x, split = "\t")[[1]],
                 simplify = T, USE.NAMES = F)
    mm <- t(mm)
    rnames <- mm[,1]
    mm <- mm[,-1]
    storage.mode(mm) <- "numeric"
    rownames(mm) <- rnames
    colnames(mm) <- cnames

    return(mm)

  } else if (grepl(network.path, pattern = ".rds")){

    cat("Reading .rds file with network coefficients\n")

    load(network.path)

    return(network.coefficients)

  } else {

    stop("Please input txt or rds file")
  }

}

#' @title Data trimming
#'
#' @usage \code{ArrangeData(data, network.path = NULL)}
#'
#' @description \code{ArrangeData} finds common genes to the network and
#' provided data and limits both datasets to these
#'
#' @param data matrix with entries equal to zero to be imputed (genes as rows
#' and samples as columns)
#' @param network.path character; path to .txt or .rds file with network
#' coefficients
#'
#' @return list; data matrix, network coefficients matrix and intercept for
#' genes common between the data matrix and the network
#'
ArrangeData <- function(data,
                        network.path = NULL){

  network.coefficients <- ReadNetwork(network.path)

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
#' @usage \code{CenterData(data, drop.exclude = T)}
#'
#' @description \code{CenterData} centers expression of each gene at 0
#'
#' @param data matrix of gene expression to be centered row-wise (genes as rows
#' and samples as columns)
#' @param drop.exclude logical; should zeros be discarded for the calculation
#' of genewise average expression levels? (defaults to T)
#'
#' @return list; row-wise centers and centered data
#'
CenterData <- function(data, drop.exclude = T){

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
#' cluster.type = "SOCK", max.iter = 50)}
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
