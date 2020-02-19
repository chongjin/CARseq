#' Estimate the prior needed for LFC shrinkage
#' 
#' \code{estimate_shrinkage_prior}, for each cell type, matches the 95 percentile of absolute values of
#' all log fold changes (LFCs)
#' to the 97.5 percentile of a zero-centered normal distribution. It also matches the 95 percentile of 
#' absolute values of all cell type-specific expressions on log scale
#' to the 97.5 percentile of a zero-centered normal distribution.
#'
#' @param estimates A matrix of K (number of cell type-independent variables) + H (number of cell types) * M (number of groups) + 1 (overdispersion)
#' columns. The estimates need to be obtained without using LFC shrinkage.
#' @param K number of cell type-independent variables
#' @param H number of cell types
#' @param M number of group
#' @param prob a number to match \code{prob} quantile of absolute value of 
#' parameter estimates to \code{prob} confidence interval
#' 
#' @return a vector of length \eqn{K + H * (M + 1)}. The first K entries are 0,
#' the next H * M entries are based on 1/sigma^2 of length H of the normal prior 
#' for each cell type-specific LFC, and the last M entries are based on
#' 1/sigma^2 of length 1 of the normal prior for
#' the cell type-specific intercept terms. 
#' 
#' @export
estimate_shrinkage_prior = function(estimates, K, M, H, prob = 0.95) {
  lambda_minimum = 1e-6
  stopifnot(is.matrix(estimates))
  stopifnot(ncol(estimates) == K + M * H + 1)
  stopifnot(M >= 1)
  # The first H entries are 1/sigma^2 of the normal prior 
  # for each cell type-specific LFC
  LFCs = matrix(NA, nrow = nrow(estimates) * M * (M - 1) * 0.5, ncol = H)
  for (m1 in 1 + seq_len(M - 1)) {
    for (m2 in seq_len(m1 - 1)) {
      index = (m1 - 1) * (m1 - 2) * 0.5 + m2
      LFCs[seq_len(nrow(estimates)) + (index - 1) * nrow(estimates), seq_len(H)] = 
          estimates[, K + (m1 - 1) * H + seq_len(H)] - estimates[, K + (m2 - 1) * H + seq_len(H)]
    }
  }
  LFC_variance_parameter = estimate_variance_parameter(LFCs, prob)
  # the last entry is 1/sigma^2 of the normal prior for
  # the cell type-specific intercept terms
  mean_expression_variance_parameter = estimate_variance_parameter(as.numeric(estimates[, K + seq_len(M*H)]), prob)
  variance_parameter = matrix(c(rep(LFC_variance_parameter, M), rep(mean_expression_variance_parameter, H)),
                              nrow = H, ncol = M + 1)
  c(rep(lambda_minimum, K), 1.0 / variance_parameter)
}

#' Estimate variance parameter of zero-mean normal prior based on a vector of numbers
#' 
#' @param estimates a vector (or a matrix) of numbers from which to estimate normal prior parameters
#' @param prob a number to match \code{prob} quantile of \code{abs(x)} to \code{prob} confidence interval
#' of normal distribution
#' 
#' @return variance parameter sigma^2 (or variance parameters for each column of a matrix)
#' of zero-mean normal distribution
estimate_variance_parameter = function(estimates, prob = 0.95) {
  estimates = as.matrix(estimates)
  apply(estimates, 2, function(x) { 
    x = abs(x[is.finite(x)])
    q = quantile(x, prob)
    (q / qnorm(1 - 0.5 * (1 - prob)))^2
  })
}