#' A simulated dataset of read counts of 50 samples and 10,000 from a mixture of three cell types.
#'
#' The first 25 samples belong to one group, and the rest 25 samples belong to the other group.
#' 
#' Within the first 2,000 genes, the first and second cell types are differentially expressed  between two groups,
#' but the third cell type is not differentially expressed. The power is generally low to detect them.
#' The rest 8,000 genes have no cell type-specific differential expression.
#' 
#' The cellular frequencies roughly centers around 6:1:3 for the three cell types.
#' 
#' The dataset is generated using the config "n_50_DE_pattern_2_2_1_replicate_1".
#'
#' @format A list of:
#' \describe{
#'   \item{observed_read_count}{A matrix of 10,000 x 50 genes generated using negative binomial distribution.}
#'   \item{rho}{A 50 x 3 matrix of cellular frequencies.}
#'   \item{clinical_variables}{A matrix of one cell type-independent covariate, RIN.}
#'   \item{d}{Numeric. A vector of 50 sample read depths.}
#' }
"n50DE221rep1"