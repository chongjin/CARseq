#' CARseq package for cell type-specific differential analysis of count data
#' 
#' @references
#'
#' @author Chong Jin, Wei Sun
#' 
#' @docType package
#' @name CARseq-package
#' @aliases CARseq-package
#' @keywords package
#' @useDynLib CARseq
#' @importFrom Rcpp sourceCpp
NULL

#' log-likelihood function of negative binomial distribution 
#' whose mean is a sum of nonnegative terms with covariates 
#'
#' @param coefs a vector of c(beta, gamma, overdispersion), where beta is a length K vector of cell type-indepedent coefficients,
#'   gamma is a matrix of H x M dimension of cell type-specific non-negative coefficients, and overdispersion is a scalar.
#' @param cell_type_specific_variables a design matrix of n_B x K of cell type-indepedent variables. Can be NULL.
#' @param other_variables a design matrix of n_B x M of cell type-specific variables.
#' @param read_depth a vector of n_B sample-specific read depths.
#' @param cellular_proportions a matrix of n_B x H of cellular proportions.
#' @param counts A vector of n_B total read counts observed.
#' 
#' @import matrixStats
#' @import nloptr
#' 
#' @examples
#' library(CARseq)
#' set.seed(1234)
#' H = 4
#' n_B = 60
#' K = 1
#' M = 3
#' coefs = c(rep(0, K), rep(1, H*M), overdispersion=10)  # initial value
#' cell_type_specific_variables = matrix(0, nrow=n_B, ncol=M)
#' cell_type_specific_variables[1:(n_B/3), 1] = 1
#' cell_type_specific_variables[(n_B/3+1):(n_B*2/3), 2] = 1
#' cell_type_specific_variables[(n_B*2/3+1):n_B, 3] = 1
#' other_variables = matrix(runif(n_B * K, min = -1, max = 1), nrow=n_B, ncol=K)
#' read_depth = round(runif(n_B, min = 50, max = 100))
#' cellular_proportions = matrix(runif(n_B * H), nrow = n_B, ncol = H)
#' cellular_proportions = cellular_proportions / rowSums(cellular_proportions)
#' counts = round(runif(n_B, min = 100, max = 200))
#' loglik(coefs, cell_type_specific_variables, other_variables, read_depth, cellular_proportions, counts)
#' @export
loglik = function(coefs, cell_type_specific_variables, other_variables, read_depth, cellular_proportions, counts) {
  # obtain H, M, K, n_B
  n_B = length(counts)
  if (is.null(other_variables)) other_variables = matrix(ncol = 0, nrow = n_B)
  if (!is.matrix(cell_type_specific_variables)) stop("cell_type_specific_variables should be a matrix")
  if (!is.matrix(other_variables)) stop("other_variables should be a matrix")
  if (!is.matrix(cellular_proportions)) stop("cellular_proportions should be a matrix")
  if (nrow(cell_type_specific_variables) != n_B) stop("Number of rows in cell_type_specific_variables does not match the length of counts")
  if (nrow(other_variables) != n_B) stop("Number of rows in other_variables does not match the length of counts")
  if (nrow(cellular_proportions) != n_B) stop("Number of rows in cellular_proportions does not match the length of counts")
  M = ncol(cell_type_specific_variables)
  K = ncol(other_variables)
  H = ncol(cellular_proportions)
  if (length(coefs) != K + H * M + 1) stop("The length of coefs is not equal to K + H * M + 1")
  
  # split coefs into beta (cell type-independent) and gamma (cell type-specific)
  beta = coefs[seq_len(K)]
  gamma = matrix(coefs[seq_len(H * M) + K], nrow = H, ncol = M)
  overdispersion = coefs[length(coefs)]
  
  covariate_adjusted_read_depth = exp(other_variables %*% matrix(beta, ncol=1, nrow=K)) * read_depth  # n_B x 1
  mu = rep(NA, n_B)
  for (i in seq_len(n_B)) {
    # read depth, other effects, cell type-specific effects
    mu[i] = sum(cellular_proportions[i, ] * matrixStats::rowProds(gamma ^ rep(cell_type_specific_variables[i, ], each=H)))
  }
  mu = mu * covariate_adjusted_read_depth
  sum(lgamma(counts + overdispersion) - lgamma(overdispersion) - lgamma(counts + 1) +
        overdispersion*log(overdispersion) + counts*log(mu) -
        (counts + overdispersion) * log(overdispersion+mu))
}

#' gradient of log-likelihood function of negative binomial distribution 
#' whose mean is a sum of nonnegative terms with covariates 
#'
#' @param coefs a vector of c(beta, gamma, overdispersion), where beta is a length K vector of cell type-indepedent coefficients,
#'   gamma is a matrix of H x M dimension of cell type-specific non-negative coefficients, and overdispersion is a scalar.
#' @param cell_type_specific_variables a design matrix of n_B x K of cell type-indepedent variables. Can be NULL.
#' @param other_variables a design matrix of n_B x M of cell type-specific variables.
#' @param read_depth a vector of n_B sample-specific read depths.
#' @param cellular_proportions a matrix of n_B x H of cellular proportions.
#' @param counts A vector of n_B total read counts observed.
#' @examples
#' library(CARseq)
#' set.seed(1234)
#' H = 4
#' n_B = 60
#' K = 1
#' M = 3
#' coefs = c(rep(0, K), rep(1, H*M), overdispersion=10)  # initial value
#' cell_type_specific_variables = matrix(0, nrow=n_B, ncol=M)
#' cell_type_specific_variables[1:(n_B/3), 1] = 1
#' cell_type_specific_variables[(n_B/3+1):(n_B*2/3), 2] = 1
#' cell_type_specific_variables[(n_B*2/3+1):n_B, 3] = 1
#' other_variables = matrix(runif(n_B * K, min = -1, max = 1), nrow=n_B, ncol=K)
#' read_depth = round(runif(n_B, min = 50, max = 100))
#' cellular_proportions = matrix(runif(n_B * H), nrow = n_B, ncol = H)
#' cellular_proportions = cellular_proportions / rowSums(cellular_proportions)
#' counts = round(runif(n_B, min = 100, max = 200))
#' grad_loglik(coefs, cell_type_specific_variables, other_variables, read_depth, cellular_proportions, counts)
#' @export
grad_loglik = function(coefs, cell_type_specific_variables, other_variables, read_depth, cellular_proportions, counts) {
  # obtain H, M, K, n_B
  n_B = length(counts)
  if (is.null(other_variables)) other_variables = matrix(ncol = 0, nrow = n_B)
  if (!is.matrix(cell_type_specific_variables)) stop("cell_type_specific_variables should be a matrix")
  if (!is.matrix(other_variables)) stop("other_variables should be a matrix")
  if (!is.matrix(cellular_proportions)) stop("cellular_proportions should be a matrix")
  if (nrow(cell_type_specific_variables) != n_B) stop("Number of rows in cell_type_specific_variables does not match the length of counts")
  if (nrow(other_variables) != n_B) stop("Number of rows in other_variables does not match the length of counts")
  if (nrow(cellular_proportions) != n_B) stop("Number of rows in cellular_proportions does not match the length of counts")
  M = ncol(cell_type_specific_variables)
  K = ncol(other_variables)
  H = ncol(cellular_proportions)
  if (length(coefs) != K + H * M + 1) stop("The length of coefs is not equal to K + H * M + 1")
  # split coefs into beta (cell type-independent) and gamma (cell type-specific)
  beta = coefs[seq_len(K)]
  gamma = matrix(coefs[seq_len(H * M) + K], nrow = H, ncol = M)
  overdispersion = coefs[length(coefs)]
  
  g = rep(NA, K + H * M + 1)
  covariate_adjusted_read_depth = exp(other_variables %*% matrix(beta, ncol=1, nrow=K)) * read_depth  # n_B x 1
  mu = rep(NA, n_B)
  for (i in seq_len(n_B)) {
    # read depth, other effects, cell type-specific effects
    mu[i] = sum(cellular_proportions[i, ] * matrixStats::rowProds(gamma ^ rep(cell_type_specific_variables[i, ], each=H)))
  }
  mu = mu * covariate_adjusted_read_depth
  
  frac_difference = counts / mu - (counts + overdispersion) / (overdispersion + mu)
  
  # derivative wrt cell type-independent variables
  partial_l_partial_beta_i = matrix(NA, nrow = n_B, ncol = K)
  for (i in seq_len(n_B)) {
    partial_l_partial_beta_i[i, ] = frac_difference[i] * other_variables[i, ] * covariate_adjusted_read_depth[i] *
      sum(cellular_proportions[i, ] * matrixStats::rowProds(gamma ^ rep(cell_type_specific_variables[i, ], each=H)))
  }
  g[seq_len(K)] = colSums(partial_l_partial_beta_i)
  
  # derivative wrt cell type-specific variables
  partial_l_partial_gamma = matrix(NA, nrow = H, ncol = M)
  for (m in seq_len(M)) {
    partial_l_partial_gamma_i = matrix(NA, nrow = n_B, ncol = H)
    
    for (i in seq_len(n_B)) {
      cell_type_specific_variables_i = cell_type_specific_variables[i, ]
      cell_type_specific_variables_i[m] = cell_type_specific_variables_i[m] - 1
      partial_l_partial_gamma_i[i, ] = frac_difference[i] * covariate_adjusted_read_depth[i] *
        cellular_proportions[i, ] * cell_type_specific_variables[i, m] * 
        matrixStats::rowProds(gamma ^ rep(cell_type_specific_variables_i, each=H))
    }
    g[seq(H) + K + H * (m-1)] = colSums(partial_l_partial_gamma_i)
  }
  g[K + H * M + 1] = sum(digamma(counts + overdispersion) - digamma(overdispersion) +
                           log(overdispersion) + 1.0 - (counts + overdispersion) / (overdispersion + mu) -
                           log(overdispersion + mu))
  g
}
