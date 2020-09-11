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


# code from scImpute ----------------------------------------------------------

read_count <-
  function (filetype, path, out_dir, type, genelen)
  {
    if(filetype == "csv") {
      raw_count = utils::read.csv(path, header = TRUE, row.names = 1)
    }else if(filetype == "txt") {
      raw_count = utils::read.table(path, header = TRUE, row.names = 1)
    }else if(filetype == "rds") {
      raw_count = readRDS(path)
    }else{
      print("filetype can be 'csv', 'txt', or 'rds'!")
      stop()
    }
    raw_count = as.matrix(raw_count)
    print(paste("number of genes in raw count matrix", nrow(raw_count)))
    print(paste("number of cells in raw count matrix", ncol(raw_count)))

    if(type == "TPM"){
      if(length(genelen) != nrow(raw_count))
        stop("number of genes in 'genelen' and count matrix do not match! ")
      raw_count = sweep(raw_count, 1, genelen, FUN = "*")
    }

    totalCounts_by_cell = colSums(raw_count)
    saveRDS(totalCounts_by_cell,
            file = paste0(out_dir, "totalCounts_by_cell.rds"))
    totalCounts_by_cell[totalCounts_by_cell == 0] = 1
    raw_count = sweep(raw_count, MARGIN = 2,
                      10^6/totalCounts_by_cell,
                      FUN = "*")
    if (min(raw_count) < 0) {
      stop("smallest read count cannot be negative!")
    }
    count_lnorm = log10(raw_count + 1.01)
    return(count_lnorm)
  }


find_va_genes = function(parslist, subcount){
  point = log10(1.01)
  valid_genes = which( (rowSums(subcount) > point * ncol(subcount)) &
                         stats::complete.cases(parslist) )
  if(length(valid_genes) == 0) return(valid_genes)
  # find out genes that violate assumption
  mu = parslist[, "mu"]
  sgene1 = which(mu <= log10(1+1.01))
  # sgene2 = which(mu <= log10(10+1.01) & mu - parslist[,5] > log10(1.01))

  dcheck1 = stats::dgamma(mu+1, shape = parslist[, "alpha"],
                          rate = parslist[, "beta"])
  dcheck2 = stats::dnorm(mu+1, mean = parslist[, "mu"],
                         sd = parslist[, "sigma"])
  sgene3 = which(dcheck1 >= dcheck2 & mu <= 1)
  sgene = union(sgene1, sgene3)
  valid_genes = setdiff(valid_genes, sgene)
  return(valid_genes)
}


calculate_weight <-
  function (x, paramt)
  {
    pz1 = paramt[1] * stats::dgamma(x, shape = paramt[2], rate = paramt[3])
    pz2 = (1 - paramt[1]) * stats::dnorm(x, mean = paramt[4], sd = paramt[5])
    pz = pz1/(pz1 + pz2)
    pz[pz1 == 0] = 0
    return(cbind(pz, 1 - pz))
  }

# -----------------------------------------------------------------------------


GetDropoutProbabilities <- function(infile, count_path, out_dir, type, genelen,
                                    drop_thre, data){

  count <- as.matrix(read_count(filetype = infile, path = count_path,
                                out_dir = out_dir, type = type,
                                genelen = genelen))

  clust <- readRDS(paste0(out_dir, "clust.rds"))

  drop_probs <- list()

  for(cc in seq_len(max(clust, na.rm = TRUE))){

    # code from scImpute ------------------------------------------------------
    cells = which(clust == cc)

    parslist = readRDS(paste0(out_dir, "pars", cc, ".rds"))

    print("searching for valid genes ...")
    valid_genes = find_va_genes(parslist, subcount = count[, cells])
    subcount = count[valid_genes, cells, drop = FALSE]
    parslist = parslist[valid_genes, , drop = FALSE]

    droprate = t(sapply(seq_along(valid_genes), function(i) {
      wt = calculate_weight(subcount[i, ], parslist[i, ])
      return(wt[, 1])
    }))
    mucheck = sweep(subcount, MARGIN = 1, parslist[, "mu"], FUN = ">")
    droprate[mucheck & droprate > drop_thre] = 0
    # -------------------------------------------------------------------------

    to.out <- droprate
    rownames(to.out) <- rownames(subcount)
    drop_probs[[cc]] <- to.out
  }

  droprob_mat <- matrix(data = NA, nrow = nrow(data), ncol = ncol(data),
                        dimnames = dimnames(data))
  for(r in drop_probs){
    droprob_mat[rownames(r),colnames(r)] <- r
  }
  saveRDS(droprob_mat, "scImpute_dropout_probabilities.rds")

  return(droprob_mat)
}


SetBiologicalZeros <- function(imputation, drop_probs, thre, was_zero){

  if(!all.equal(dim(imputation), dim(drop_probs))){
    # limit both to the set of common genes / cells
    common_genes <- intersect(rownames(imputation), rownames(drop_probs))
    common_cells <- intersect(colnames(imputation), colnames(drop_probs))
    imputation <- imputation[common_genes, common_cells]
    drop_probs <- drop_probs[common_genes, common_cells]
  }

  # do not impute entries with NA probability (give prob = 0)
  drop_probs[is.na(drop_probs)] <- 0

  # set prob < threshold to zero
  set_to_zero <- drop_probs < thre
  imputed_entries <- was_zero & (imputation != 0) # limit to imputed values
  imputation[set_to_zero & imputed_entries] <- 0

  return(imputation)
}
