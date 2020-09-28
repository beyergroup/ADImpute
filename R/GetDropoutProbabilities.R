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

#' @title Get dropout probabilities
#'
#' @description \code{GetDropoutProbabilities} computes dropout probabilities
#' (probability of being a dropout that should be imputed rather than a true
#' biological zero) using an adaptation of scImpute's approach
#'
#' @usage GetDropoutProbabilities(data, thre, cell.clusters, labels = NULL,
#' ncores)
#'
#' @param data matrix; original data before imputation
#' @param thre numeric; probability threshold to classify entries as biological
#' zeros
#' @param cell.clusters integer; number of cell subpopulations
#' @param labels character; vector specifying the cell type of each column of
#' \code{data}
#' @param ncores integer; number of cores used for paralell computation
#'
#' @details This function follows scImpute's model to distinguish between
#' true biological zeros and dropouts, and is based on adapted code from the
#' scImpute R package.
#'
#' @return matrix with same dimensions as \code{data} containing the dropout
#' probabilities for the corresponding entries
#'
GetDropoutProbabilities <- function(data, thre, cell.clusters = 2,
                                    labels = NULL, ncores){

  labeled <- !is.null(labels)

  if(!labeled){
    dropmat = imputation_model8(count = data, labeled = labeled,
                                point = log10(1.01), drop_thre = thre,
                                Kcluster = cell.clusters, ncores = ncores)
  }else{
    dropmat = imputation_wlabel_model8(count = data, labeled = labeled,
                                       cell_labels = labels, point = log10(1.01),
                                       drop_thre = thre, ncores = ncores,
                                       Kcluster = NULL)
  }

  return(dropmat)
}



#' @title Set biological zeros
#'
#' @description \code{SetBiologicalZeros} sets some of the entries back to zero
#' after dropout imputation, as they likely correspond to true biological zeros
#' (genes not expressed in given cell)
#'
#' @usage SetBiologicalZeros(imputation, drop_probs, thre = .2, was_zero)
#'
#' @param imputation matrix; imputed values
#' @param drop_probs matrix; dropout probabilities for each entry in
#' \code{imputation}. 0 means certain biological zero, while 1 means certain
#' dropout to be imputed
#' @param thre numeric; probability threshold to classify entries as biological
#' zeros
#' @param was_zero matrix; logical matrix: was the corresponding entry of
#' \code{imputation} originally a zero?
#'
#' @details Entries which were originally zero and have dropout probability
#' below \code{thre} are considered biological zeros and, if they were imputed,
#' are set back to 0.
#'
#' @return matrix containing likely biological zeros set back to 0.
#'
SetBiologicalZeros <- function(imputation, drop_probs, thre = .2, was_zero){

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
