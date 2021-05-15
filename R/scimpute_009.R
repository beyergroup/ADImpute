
# functions adapted from https://github.com/Vivianstats/scImpute version
# 0.0.9 see publication: https://www.nature.com/articles/s41467-018-03405-7

calculate_weight <- function(x, paramt) {
  pz1 <- paramt[1] * stats::dgamma(x, shape = paramt[2], rate = paramt[3])
  pz2 <- (1 - paramt[1]) * stats::dnorm(x, mean = paramt[4], sd = paramt[5])
  pz <- pz1/(pz1 + pz2)
  pz[pz1 == 0] <- 0
  return(cbind(pz, 1 - pz))
}


dmix <- function(x, pars) {
  pars[1] * stats::dgamma(x, shape = pars[2], rate = pars[3]) +
    (1 - pars[1]) * stats::dnorm(x, mean = pars[4], sd = pars[5])
}


find_hv_genes <- function(count, I, J) {
  count_nzero <- lapply(seq_len(I),
                        function(i) setdiff(count[i, ], log10(1.01)))
  mu <- vapply(count_nzero, mean, FUN.VALUE = 1)
  mu[is.na(mu)] <- 0
  sd <- vapply(count_nzero, stats::sd, FUN.VALUE = 1)
  sd[is.na(sd)] <- 0
  cv <- sd/mu
  cv[is.na(cv)] <- 0
  high_var_genes <- which(mu >= 1 & cv >= stats::quantile(cv, 0.25))
  if (length(high_var_genes) < 500) {
    high_var_genes <- seq_len(I)
  }
  count_hv <- count[high_var_genes, ]
  return(count_hv)
}


find_neighbors_labeled <- function(count_hv, J, Kcluster = NULL, cores,
                                   BPPARAM, cell_labels = NULL) {
  if (is.character(cell_labels)) {
    labels_uniq <- unique(cell_labels); labels_mth <- seq_along(labels_uniq)
    names(labels_mth) <- labels_uniq; clust <- labels_mth[cell_labels]
  } else { clust <- cell_labels }
  nclust <- length(unique(clust)); message("calculating cell distances ...")
  dist_list <- lapply(seq_len(nclust), function(ll) {
    cell_inds <- which(clust == ll)
    count_hv_sub <- count_hv[, cell_inds, drop = FALSE]
    if (length(cell_inds) < 1000) {
      var_thre <- 0.4; pca <- stats::prcomp(t(count_hv_sub))
      eigs <- (pca$sdev)^2; var_cum <- cumsum(eigs)/sum(eigs)
      if (max(var_cum) <= var_thre) { npc <- length(var_cum)
      } else { npc <- which.max(var_cum > var_thre) }
    } else { var_thre <- 0.6
    pca <- rsvd::rpca(t(count_hv_sub), k = 1000, center = TRUE,
                      scale = FALSE)
    eigs <- (pca$sdev)^2; var_cum <- cumsum(eigs)/sum(eigs)
    if (max(var_cum) <= var_thre) { npc <- length(var_cum)
    } else { npc <- which.max(var_cum > var_thre) }
    }
    if (npc < 3) { npc <- 3 }
    mat_pcs <- t(pca$x[, seq_len(npc)])
    dist_cells_list <- BiocParallel::bplapply(seq_along(cell_inds),
                                              function(id1) { d <- vapply(seq_len(id1), function(id2) {
                                                sse <- sum((mat_pcs[, id1] - mat_pcs[, id2])^2)
                                                sqrt(sse) }, FUN.VALUE = 1)
                                              return(c(d, rep(0, length(cell_inds) - id1))) }, BPPARAM = BPPARAM)
    dist_cells <- matrix(0, nrow = length(cell_inds),
                         ncol = length(cell_inds))
    for (cellid in seq_along(cell_inds)) {
      dist_cells[cellid, ] <- dist_cells_list[[cellid]]
    }; dist_cells <- dist_cells + t(dist_cells)
    return(dist_cells)})
  return(list(dist_list = dist_list, clust = clust))
}


find_neighbors_unlabeled <- function(count_hv, J, Kcluster = NULL, cores,
                                     BPPARAM, cell_labels = NULL) {

  message("dimension reduction ...")
  if (J < 5000) { var_thre <- 0.4; pca <- stats::prcomp(t(count_hv))
  eigs <- (pca$sdev)^2; var_cum <- cumsum(eigs)/sum(eigs)
  if (max(var_cum) <= var_thre) { npc <- length(var_cum)
  } else {
    npc <- which.max(var_cum > var_thre); npc <- max(npc, Kcluster) }
  } else { var_thre <- 0.6
  pca <- rsvd::rpca(t(count_hv), k = 1000, center = TRUE, scale = FALSE)
  eigs <- (pca$sdev)^2; var_cum <- cumsum(eigs)/sum(eigs)
  if (max(var_cum) <= var_thre) { npc <- length(var_cum)
  } else {
    npc <- which.max(var_cum > var_thre); npc <- max(npc, Kcluster) }
  }; if (npc < 3) { npc <- 3 }
  mat_pcs <- t(pca$x[, seq_len(npc)])  # columns are cells
  message("calculating cell distances ...")
  dist_cells_list <- BiocParallel::bplapply(seq_len(J), function(id1) {
    d <- vapply(seq_len(id1), function(id2) {
      sse <- sum((mat_pcs[, id1] - mat_pcs[, id2])^2)
      sqrt(sse) }, FUN.VALUE = 1)
    return(c(d, rep(0, J - id1))) }, BPPARAM = BPPARAM)
  dist_cells <- matrix(0, nrow = J, ncol = J)
  for (cellid in seq_len(J)) {
    dist_cells[cellid, ] <- dist_cells_list[[cellid]] }
  dist_cells <- dist_cells + t(dist_cells)
  min_dist <- vapply(seq_len(J), function(i) { min(dist_cells[i, -i]) },
                     FUN.VALUE = 1)
  iqr <- stats::quantile(min_dist, 0.75) - stats::quantile(min_dist, 0.25)
  outliers <- which(min_dist > 1.5 * iqr + stats::quantile(min_dist, 0.75))
  non_out <- setdiff(seq_len(J), outliers)
  spec_res <- kernlab::specc(t(mat_pcs[, non_out]), centers = Kcluster,
                             kernel = "rbfdot")
  message("cluster sizes:"); message(spec_res@size)
  nbs <- rep(NA, J); nbs[non_out] <- spec_res
  return(list(dist_cells = dist_cells, clust = nbs))
}


find_va_genes <- function(parslist, subcount) {
  point <- log10(1.01)
  valid_genes <- which((rowSums(subcount) > point * ncol(subcount)) &
                         stats::complete.cases(parslist))
  if (length(valid_genes) == 0)
    return(valid_genes)
  # find out genes that violate assumption
  mu <- parslist[, "mu"]
  sgene1 <- which(mu <= log10(1 + 1.01))

  dcheck1 <- stats::dgamma(mu + 1, shape = parslist[, "alpha"],
                           rate = parslist[, "beta"])
  dcheck2 <- stats::dnorm(mu + 1, mean = parslist[, "mu"], sd = parslist[,
                                                                         "sigma"])
  sgene3 <- which(dcheck1 >= dcheck2 & mu <= 1)
  sgene <- union(sgene1, sgene3)
  valid_genes <- setdiff(valid_genes, sgene)
  return(valid_genes)
}

### root-finding equation
fn <- function(alpha, target) {
  log(alpha) - digamma(alpha) - target
}


### estimate parameters in the mixture distribution
get_mix <- function(xdata, point) {
  inits <- rep(0, 5)
  inits[1] <- sum(xdata == point)/length(xdata)
  if (inits[1] == 0) {
    inits[1] <- 0.01
  }
  inits[2:3] <- c(0.5, 1)
  xdata_rm <- xdata[xdata > point]
  inits[4:5] <- c(mean(xdata_rm), stats::sd(xdata_rm))
  if (is.na(inits[5])) {
    inits[5] <- 0
  }
  paramt <- inits
  eps <- 10
  iter <- 0
  loglik_old <- 0

  while (eps > 0.5) {
    wt <- calculate_weight(xdata, paramt)
    paramt[1] <- sum(wt[, 1])/nrow(wt)
    paramt[4] <- sum(wt[, 2] * xdata)/sum(wt[, 2])
    paramt[5] <- sqrt(sum(wt[, 2] * (xdata - paramt[4])^2)/sum(wt[, 2]))
    paramt[2:3] <- update_gmm_pars(x = xdata, wt = wt[, 1])

    loglik <- sum(log10(dmix(xdata, paramt)))
    eps <- (loglik - loglik_old)^2
    loglik_old <- loglik
    iter <- iter + 1
    if (iter > 100) {
      break
    }
  }
  return(paramt)
}

get_mix_parameters <- function(count, point = log10(1.01), path,
                               cores, BPPARAM) {

  count <- as.matrix(count)
  null_genes <- which(abs(rowSums(count) - point * ncol(count)) < 1e-10)
  parslist <- BiocParallel::bplapply(seq_len(nrow(count)), function(ii) {
    if (ii%%2000 == 0) {
      gc()
      message(ii)
    }
    if (ii %in% null_genes) {
      return(rep(NA, 5))
    }
    xdata <- count[ii, ]
    paramt <- tryCatch(get_mix(xdata, point),
                       error = function(e) rep(NA, 5))
    return(paramt)
  }, BPPARAM = BPPARAM)
  parslist <- Reduce(rbind, parslist)
  colnames(parslist) <- c("rate", "alpha", "beta", "mu", "sigma")
  return(parslist)
}


imputation_model8 <- function(count, labeled = FALSE, point, drop_thre = 0.5,
                              Kcluster = 10, cores = BiocParallel::bpworkers(BPPARAM),
                              BPPARAM = BiocParallel::SnowParam(type = "SOCK")) {

  count <- as.matrix(count); I <- nrow(count)
  J <- ncol(count); count_imp <- count; count_hv <- find_hv_genes(count, I, J)
  message("searching candidate neighbors ... ")
  if (Kcluster == 1) { clust <- rep(1, J)
  if (J < 5000) { var_thre <- 0.4; pca <- stats::prcomp(t(count_hv))
  eigs <- (pca$sdev)^2; var_cum <- cumsum(eigs)/sum(eigs)
  if (max(var_cum) <= var_thre) {  npc <- length(var_cum)
  } else { npc <- which.max(var_cum > var_thre)
  npc <- max(npc, Kcluster) }
  } else { var_thre <- 0.6
  pca <- rsvd::rpca(t(count_hv), k = 1000, center = TRUE,
                    scale = FALSE)
  eigs <- (pca$sdev)^2; var_cum <- cumsum(eigs)/sum(eigs)
  if (max(var_cum) <= var_thre) { npc <- length(var_cum)
  } else {
    npc <- which.max(var_cum > var_thre); npc <- max(npc,Kcluster) }
  }
  if (npc < 3) { npc <- 3 }
  mat_pcs <- t(pca$x[, seq_len(npc)])  # columns are cells
  dist_cells_list <- BiocParallel::bplapply(seq_len(J), function(id1) {
    d <- vapply(seq_len(id1), function(id2) {
      sse <- sum((mat_pcs[, id1] - mat_pcs[, id2])^2)
      sqrt(sse) }, FUN.VALUE = 1); return(c(d, rep(0, J - id1))) },
    BPPARAM = BPPARAM)
  dist_cells <- matrix(0, nrow = J, ncol = J)
  for (cellid in seq_len(J)) {
    dist_cells[cellid, ] <- dist_cells_list[[cellid]] }
  dist_cells <- dist_cells + t(dist_cells)
  } else { message("inferring cell similarities ...")
    neighbors_res <- find_neighbors_unlabeled(count_hv = count_hv, J = J,
                                              Kcluster = Kcluster, cores = cores, BPPARAM = BPPARAM)
    dist_cells <- neighbors_res$dist_cells; clust <- neighbors_res$clust }
  return(MixtureModel(count, clust, cores, BPPARAM = BPPARAM, drop_thre))
}


imputation_wlabel_model8 <- function(count, labeled, cell_labels = NULL, point,
                                     drop_thre, Kcluster = NULL, cores = BiocParallel::bpworkers(BPPARAM),
                                     BPPARAM = BiocParallel::SnowParam(type = "SOCK")) {

  if (!(is.character(cell_labels) | is.numeric(cell_labels) |
        is.integer(cell_labels))) {
    stop("cell_labels should be a character or integer vector!")}

  count <- as.matrix(count); I <- nrow(count)
  J <- ncol(count); count_imp <- count
  count_hv <- find_hv_genes(count, I, J)
  message("searching candidate neighbors ... ")
  neighbors_res <- find_neighbors_labeled(count_hv = count_hv, J = J,
                                          cores = cores, BPPARAM = BPPARAM, cell_labels = cell_labels)
  dist_list <- neighbors_res$dist_list; clust <- neighbors_res$clust

  nclust <- sum(!is.na(unique(clust)))

  droprate <- list()
  for (cc in seq_len(nclust)) {
    message(paste("estimating dropout probability for type", cc, "..."))
    parslist <- get_mix_parameters(count = count[, which(clust == cc),
                                                 drop = FALSE], point = log10(1.01), cores = cores,
                                   BPPARAM = BPPARAM)
    cells <- which(clust == cc); if (length(cells) <= 1) { next }
    message("searching for valid genes ...")
    valid_genes <- find_va_genes(parslist, subcount = count[, cells])
    if (length(valid_genes) <= 10) { next }

    subcount <- count[valid_genes, cells, drop = FALSE]
    Ic <- length(valid_genes); Jc <- ncol(subcount)
    parslist <- parslist[valid_genes, , drop = FALSE]

    droprate[[cc]] <- t(vapply(seq_len(Ic), function(i) {
      wt <- calculate_weight(subcount[i, ], parslist[i, ])
      return(wt[, 1]) }, FUN.VALUE = rep(1, Jc)))
    dimnames(droprate[[cc]]) <- dimnames(subcount)
    mucheck <- sweep(subcount, MARGIN = 1, parslist[, "mu"], FUN = ">")
    droprate[[cc]][mucheck & droprate[[cc]] > drop_thre] <- 0
  }

  dropratem <- matrix(data = NA, nrow = nrow(count), ncol = ncol(count),
                      dimnames = dimnames(count))
  for (dmat in droprate) { dropratem[rownames(dmat), colnames(dmat)] <- dmat }

  return(dropratem)
}


MixtureModel <- function(count, clust, cores, BPPARAM, drop_thre) {

  nclust <- sum(!is.na(unique(clust)))

  droprate <- list()

  for (cc in seq_len(nclust)) {
    message(paste("estimating dropout probability for type", cc, "..."))
    parslist <- get_mix_parameters(count = count[, which(clust == cc),
                                                 drop = FALSE], point = log10(1.01), cores = cores,
                                   BPPARAM = BPPARAM)
    cells <- which(clust == cc)
    if (length(cells) <= 1) {
      next
    }
    message("searching for valid genes ...")
    valid_genes <- find_va_genes(parslist, subcount = count[, cells])
    if (length(valid_genes) <= 10) {
      message("Not enough valid genes found ...")
      next
    }

    subcount <- count[valid_genes, cells, drop = FALSE]
    Ic <- length(valid_genes)
    Jc <- ncol(subcount)
    parslist <- parslist[valid_genes, , drop = FALSE]

    droprate[[cc]] <- t(vapply(seq_len(Ic), function(i) {
      wt <- calculate_weight(subcount[i, ], parslist[i, ])
      return(wt[, 1])
    }, FUN.VALUE = rep(1, Jc)))
    dimnames(droprate[[cc]]) <- dimnames(subcount)
    mucheck <- sweep(subcount, MARGIN = 1, parslist[, "mu"], FUN = ">")
    droprate[[cc]][mucheck & droprate[[cc]] > drop_thre] <- 0
  }

  dropratem <- matrix(data = NA, nrow = nrow(count), ncol = ncol(count),
                      dimnames = dimnames(count))
  for (dmat in droprate) {
    dropratem[rownames(dmat), colnames(dmat)] <- dmat
  }

  return(dropratem)
}


read_count <- function(raw_count, type, genelen) {
  raw_count <- as.matrix(raw_count)
  message(paste("number of genes in raw count matrix", nrow(raw_count)))
  message(paste("number of cells in raw count matrix", ncol(raw_count)))

  if (type == "TPM") {
    if (length(genelen) != nrow(raw_count))
      stop("number of genes in 'genelen' and count matrix do not match! ")
    raw_count <- sweep(raw_count, 1, genelen, FUN = "*")
  }

  totalCounts_by_cell <- colSums(raw_count)
  totalCounts_by_cell[totalCounts_by_cell == 0] <- 1
  raw_count <- sweep(raw_count, MARGIN = 2, 10^6/totalCounts_by_cell,
                     FUN = "*")
  if (min(raw_count) < 0) {
    stop("smallest read count cannot be negative!")
  }
  count_lnorm <- log10(raw_count + 1.01)
  return(count_lnorm)
}


### update parameters in gamma distribution
update_gmm_pars <- function(x, wt){
  tp_s <- sum(wt)
  tp_t <- sum(wt * x)
  tp_u <- sum(wt * log(x))
  tp_v <- -tp_u/tp_s - log(tp_s/tp_t)
  if (tp_v <= 0) {
    alpha <- 20
  } else {
    alpha0 <- (3 - tp_v + sqrt((tp_v - 3)^2 + 24 * tp_v))/12/tp_v
    if (alpha0 >= 20) {
      alpha <- 20
    } else {
      alpha <- stats::uniroot(fn, c(0.9, 1.1) * alpha0, target = tp_v,
                              extendInt = "yes")$root
    }
  }
  ## need to solve log(x) - digamma(x) = tp_v We use this approximation to
  ## compute the initial value

  beta <- tp_s/tp_t * alpha
  return(c(alpha, beta))
}

