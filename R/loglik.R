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
#' @import matrixStats
#' @import nloptr
#' @import zetadiv
#' @import MASS
#' @import bvls
NULL


#' fit a CARseq negative binomial model 
#' using iteratively weighted least squares / non-negative least squares
#' with backtracking line search.
#' 
#' \code{fit_model} fits a negative binomial distribution 
#' whose mean is a sum of nonnegative terms with covariates.
#' The overdispersion parameter is estimated by maximizing the 
#' adjusted profile log-likelihood.
#'
#' @param cell_type_specific_variables an array of n_B x H x K of cell type-independent variables.
#' @param other_variables a design matrix of n_B x M of cell type-specific variables. Can be NULL.
#' @param read_depth a vector of n_B sample-specific read depths. It it used as an offset term in the
#'        regression model. Alternatively, it can be 1, NA or NULL so that we don't have an offset term
#'        in the model, and log(read_depth) can be included as one of the cell type-independent variables.
#' @param cellular_proportions a matrix of n_B x H of cellular proportions.
#' @param counts A vector of n_B total read counts observed.
#' @param init a numeric vector of (K + H x M + 1) corresponding to the initial value of
#'             c(beta, gamma, overdispersion). Can be NULL.
#' @param fix_overdispersion logical (or numerical for a fixed overdispersion). 
#'        If \code{FALSE} (default), the overdispersion parameter will be estimated within the function.
#'        Otherwise, the overdispersion parameter can be 
#' @param number_of_resample numeric. The number of initial values we use in optimization. 
#'        The optimization algorithm generally is not sensitive to initial values, so we advise to leave it
#'        as the default value 1.
#' @param use_log_scale_algorithm logical. The default is FALSE, where the cell type-specific effects are
#'        fitted on a non-log scale and non-negative least squares is used (implemented in bvls).
#'        This is generally advisable than fitting them on a log scale, especially when the cell type-specific
#'        variables are binary (for example group indicators).
#' @param verbose logical. If \code{TRUE} (default), debugging information of negative log-likelihood of each
#'        iteration will be printed, grouped by each backtracking line search.
#'
#' @examples
#' library(CARseq)
#' set.seed(1234)
#' H = 4      # number of cell types
#' n_B = 60   # total number of subjects
#' K = 1      # number of cell type-independent covariates
#' M = 3      # number of cell type-specific covariates
#' cell_type_specific_variables = array(0, dim=c(n_B, H, M))
#' cell_type_specific_variables[1:(n_B/3), , 1] = 1
#' cell_type_specific_variables[(n_B/3+1):(n_B*2/3), , 2] = 1
#' cell_type_specific_variables[(n_B*2/3+1):n_B, , 3] = 1
#' other_variables = matrix(runif(n_B * K, min = -1, max = 1), nrow=n_B, ncol=K)
#' read_depth = round(runif(n_B, min = 50, max = 100))
#' cellular_proportions = matrix(runif(n_B * H), nrow = n_B, ncol = H)
#' cellular_proportions = cellular_proportions / rowSums(cellular_proportions)
#' counts = round(runif(n_B, min = 100, max = 200))
#' res = fit_model(cell_type_specific_variables, other_variables, read_depth, cellular_proportions, counts)
#' str(res)
#' @export
fit_model = function(cell_type_specific_variables, other_variables, read_depth, cellular_proportions, counts,
                     init = NULL, fix_overdispersion = FALSE, number_of_resample = 1, use_log_scale_algorithm = FALSE, verbose = TRUE) {
  
  # obtain H, M, K, n_B
  n_B = length(counts)
  if (is.null(other_variables)) other_variables = matrix(ncol = 0, nrow = n_B)
  if (!is.array(cell_type_specific_variables)) {
    stop("cell_type_specific_variables should be an array")
  } else if (length(dim(cell_type_specific_variables)) != 3) {
    stop("cell_type_specific_variables should be a 3-dimensional (n_B, H, M) array")
  }
  if (!is.matrix(other_variables)) stop("other_variables should be a matrix")
  if (!is.matrix(cellular_proportions)) stop("cellular_proportions should be a matrix")
  if (dim(cell_type_specific_variables)[1] != n_B) stop("Length of the 1st dimension in cell_type_specific_variables does not match the length of counts")
  if (nrow(other_variables) != n_B) stop("Number of rows in other_variables does not match the length of counts")
  if (nrow(cellular_proportions) != n_B) stop("Number of rows in cellular_proportions does not match the length of counts")
  if (is.null(read_depth) || identical(read_depth, NA) || length(read_depth) == 1) {
    read_depth = rep(1, n_B)
  }
  if (length(read_depth) != n_B) stop("Length of read_depth does not match the length of counts")
  M = dim(cell_type_specific_variables)[3]
  K = ncol(other_variables)
  H = ncol(cellular_proportions)
  if (identical(fix_overdispersion, TRUE) && is.null(init)) {
    stop("When fix_overdispersion is TRUE, it needs to be specified in init. Otherwise, you can directly specify fix_overdispersion as a numeric.")
  }
  coefs = init

  # We calculate which (cell type, cell type-specific variable) pairs are not included in the model
  # For future use.
  # is_active is a logical vector of length (K + H x M), recording which predictors are included in the model.
  is_reduced = apply(cell_type_specific_variables, c(2,3),
                     function(y) identical(y, rep(0, dim(cell_type_specific_variables)[1])))
  is_active = c(rep(TRUE, K), !is_reduced)
  
  # Upper limit and lower limit on entries on log scale
  log_lower = -30
  log_upper = 30  # some gene expression are extremely high
  # Upper limit and lower limit on entries of gamma
  if (use_log_scale_algorithm) {
    gamma_lower = log_lower
    gamma_upper = log_upper
  } else {
    gamma_lower = exp(log_lower)
    gamma_upper = exp(log_upper)
  }
  combined_lower = c(rep(-exp(log_upper), K), rep(gamma_lower, sum(is_active) - K))
  combined_upper = c(rep(exp(log_upper), K), rep(gamma_upper, sum(is_active) - K))
  
  # We need to sample many times to decide on an initial value.
  if (is.null(number_of_resample) || is.na(number_of_resample)) number_of_resample = 20
  # We only need number_of_resample different initial values and subsequent model fits.
  # However, we can try at most number_of_resample_max times.
  # We finally decide to keep them the same as the optimization algorithm is already good enough.
  number_of_resample_max = number_of_resample
  resample_size_list = rep(3 * (K + H * M), number_of_resample_max)    # An arbitrary number that should be larger than the number of parameters in the model
  resample_size_list[1] = n_B    # in the first instance, the initial value is based on all samples instead of a small subset of samples
  
  if (is.numeric(fix_overdispersion)) {
    overdispersion_theta = fix_overdispersion
    fix_overdispersion = TRUE
  } else if (fix_overdispersion) {
    overdispersion_theta = coefs[length(coefs)]
  }
  
  negloglik_all_time_best = NULL
  coefs_all_time_best = NULL
  number_of_successful_iter = 0
  
  # Need to generate multiple initial values. 
  # If we have initial value, we need to use it first.
  negloglik_curr_vector = rep(NA, number_of_resample_max)  # stores a vector of likelihood for debugging purposes
  for (iter in seq_len(number_of_resample_max)) {
    resample_size = resample_size_list[iter]
    
    # save the current likelihood for each optimization path with a new initial value
    negloglik_curr_vector[iter] = tryCatch({
      
      # Generate initial value through glm model.
      # In very few cases, the overdispersion parameter will be extremely large, and we need to do the sampling again
      # to regenerate initial values.
      max_init_iter = 10
      for (init_iter in 1:max_init_iter) {
      
        # sample from a matrix of cell type-specific covariates 
        array_inds_sampled_cell_type_specific = as.matrix(merge(cbind(1:n_B, sample(x = 1:H, size = n_B, replace = TRUE)), seq_len(M), by=NULL))
        variables_sampled_cell_type_specific = matrix(cell_type_specific_variables[array_inds_sampled_cell_type_specific], nrow=n_B, ncol=M)
        
        # Do we already have an overdispersion parameter yet?
        if (iter == 1 & (!fix_overdispersion)) {
          # Are there cell type-dependent variables?
          if (K == 0) {
            nbmodel = MASS::glm.nb(counts~offset(log(read_depth)) + 0 + variables_sampled_cell_type_specific,
                                   subset = sample(c(rep(TRUE, resample_size), rep(FALSE, n_B - resample_size))))
          } else {
            nbmodel = MASS::glm.nb(counts~offset(log(read_depth)) + 0 + other_variables + variables_sampled_cell_type_specific, 
                                   subset = sample(c(rep(TRUE, resample_size), rep(FALSE, n_B - resample_size))))
          }
          overdispersion_theta = nbmodel$theta
        } else {
          if (K == 0) {
            nbmodel = stats::glm(counts~offset(log(read_depth)) + 0 + variables_sampled_cell_type_specific,
                                   subset = sample(c(rep(TRUE, resample_size), rep(FALSE, n_B - resample_size))),
                          family = negative.binomial(overdispersion_theta))
          } else {
            nbmodel = stats::glm(counts~offset(log(read_depth)) + 0 + other_variables + variables_sampled_cell_type_specific, 
                                   subset = sample(c(rep(TRUE, resample_size), rep(FALSE, n_B - resample_size))),
                          family = negative.binomial(overdispersion_theta))
          }
        }
        
        if (overdispersion_theta < exp(10)) break
        if (init_iter == max_init_iter) stop("Cannot find a suitable initial value!")
      }
  
      # construct initial parameters
      if (is.null(init)) {
        # initialize using MASS::glm.nb
        if (use_log_scale_algorithm) {
          coefs = c(nbmodel$coefficients[seq_len(K)],                 # cell type-independent variables
                    rep(nbmodel$coefficients[seq_len(M) + K], each=H),   # cell type-specific variables
                    overdispersion = overdispersion_theta)
        } else {
          coefs = c(nbmodel$coefficients[seq_len(K)],                 # cell type-independent variables
                    rep(exp(nbmodel$coefficients[seq_len(M) + K]), each=H),   # cell type-specific variables
                    overdispersion = overdispersion_theta)
        }  
        coefs[is.na(coefs)] = 0
      } else if (length(coefs) != K + H * M + 1) {
        stop("The length of coefs is not equal to K + H * M + 1")
      } else if (!use_log_scale_algorithm) {
        # We need to take exp to make gamma in non-log scale within "coefs"
        coefs[seq_len(H * M) + K] = exp(coefs[seq_len(H * M) + K])
      }
      
      # TODO: use cooks.distance() and fitted() to process the outliers
      # https://stats.stackexchange.com/questions/11632/what-kind-of-residuals-and-cooks-distance-are-used-for-glm
      # At the current development stage, we did not calculate cooks distance from our fitted model.
      # Instead, we fitted a glm model to replace the outliers outside the function.
      
      # Convergence conditions of Iteratively Weighted Least Squares
      epsilon_convergence = 1e-4
      maxiter = 30
      # The negative log-likelihood needs to decrease each iteration of the loop.
      # Break the loop if it does not decrease a lot, or it increases:
      inner_iter = 0
      negloglik_prev = negloglik_curr = negloglik_best = Inf

      while (inner_iter == 0 ||
             (negloglik_prev - negloglik_curr > epsilon_convergence
              &&
              inner_iter < maxiter)) {
        
        # split coefs into beta (cell type-independent) and gamma (cell type-specific)
        beta = coefs[seq_len(K)]
        gamma = matrix(coefs[seq_len(H * M) + K], nrow = H, ncol = M)
        covariate_adjusted_read_depth = exp(other_variables %*% matrix(beta, ncol=1, nrow=K)) * read_depth  # n_B x 1
        mu_matrix = matrix(NA, n_B, H)
        for (i in seq_len(n_B)) {
          # read depth, other effects, cell type-specific effects
          if (use_log_scale_algorithm) {
            mu_matrix[i, ] = cellular_proportions[i, ] * exp(rowSums(gamma * cell_type_specific_variables[i, , ]))
          } else {
            mu_matrix[i, ] = cellular_proportions[i, ] * matrixStats::rowProds(gamma ^ cell_type_specific_variables[i, , ])
          }
        }
        mu_matrix = mu_matrix * as.numeric(covariate_adjusted_read_depth)
        mu = rowSums(mu_matrix)
        
        # current negative log-likelihood:
        negloglik_prev = negloglik_curr
        negloglik_curr = - sum(lgamma(counts + overdispersion_theta) - lgamma(overdispersion_theta) - lgamma(counts + 1) +
                                 overdispersion_theta*log(overdispersion_theta) + counts*log(mu) -
                                 (counts + overdispersion_theta) * log(overdispersion_theta+mu))
        if (!is.finite(negloglik_curr)) {
          # message("The optimization failed to converge.")
          break
        }
        
        # in weighted least squares:
        weights = overdispersion_theta / (mu * (mu+overdispersion_theta))
        if (use_log_scale_algorithm) {
          adjusted_design_matrix = cbind(other_variables * mu,
                                         matrix(cell_type_specific_variables, nrow=n_B) * as.numeric(mu_matrix))
        } else {
          adjusted_design_matrix = cbind(other_variables * mu,
                                         matrix(cell_type_specific_variables, nrow=n_B) / rep(gamma, each=n_B) * as.numeric(mu_matrix))
        }
        
        adjusted_response = adjusted_design_matrix %*% coefs[-length(coefs)] + (counts - mu)
        
        # score function
        score_curr = crossprod(adjusted_design_matrix[, is_active], weights * (counts - mu))
        
        # least squares
        if (use_log_scale_algorithm) {
          direction = - coefs[-length(coefs)][is_active] +
            tryCatch(stats::lsfit(adjusted_design_matrix[, is_active],
                                  adjusted_response,
                                  wt = weights,
                                  intercept = FALSE)$coef, 
                               warning = function(w) {coefs[-length(coefs)][is_active]},
                               error = function(e) {coefs[-length(coefs)][is_active]})
        } else {
          # bounded-variable least squares ensures that the active parameters in H*M cell type-specific ones are bounded
          direction = - coefs[-length(coefs)][is_active] +
            bvls::bvls(adjusted_design_matrix[, is_active] * sqrt(weights),
                       adjusted_response * sqrt(weights),
                       bl = combined_lower,
                       bu = combined_upper)$x
        }
        
        # backtracking line search (halving)
        c1 = 1e-4
        contraction = 0.5
        step_length = 1
        
        max_backtracking_iter = 10
        for (backtracking_iter in 1:max_backtracking_iter) {
          coefs_new = coefs
          coefs_new[-length(coefs_new)][is_active] = coefs[-length(coefs)][is_active] + step_length * direction
          
          # if the parameters are already out of the boundary, we half the step size until we are back inside the boundary:
          step_iter = 1
          max_step_iter = 20
          while (backtracking_iter == 1 &&
                 (any(combined_upper < coefs_new[-length(coefs_new)][is_active]) ||
                  any(combined_lower > coefs_new[-length(coefs_new)][is_active])) &&
                 step_iter < max_step_iter) {
            step_length = step_length * contraction
            coefs_new[-length(coefs_new)][is_active] = coefs[-length(coefs)][is_active] + step_length * direction
            step_iter = step_iter + 1
          }
          
          # Calculate the log-likelihood at the intended coefs_new:
          negloglik_new = negloglik(coefs = coefs_new,
                                  cell_type_specific_variables = cell_type_specific_variables,
                                  other_variables = other_variables,
                                  read_depth = read_depth,
                                  cellular_proportions = cellular_proportions,
                                  counts = counts,
                                  use_log_scale_params = use_log_scale_algorithm,
                                  verbose = FALSE)

          if (verbose) message(sprintf("Negative log-likelihood = %.4f at iteration %d",
                                       negloglik_new, inner_iter))
          if (!is.finite(negloglik_new)) {
            negloglik_new = Inf
            break
          }
          if (negloglik_new <= negloglik_curr + c1 * step_length * crossprod(-score_curr, direction))
            break
          step_length = step_length * contraction
        }
        if (verbose) cat("\n")
        coefs = coefs_new
        negloglik_curr = negloglik_new
        
        # save the best coefs and negative log-likelihood
        if (is.null(negloglik_best) || negloglik_curr < negloglik_best) {
          coefs_best = coefs
          negloglik_best = negloglik_curr
        }
        if (is.null(negloglik_all_time_best) || negloglik_curr < negloglik_all_time_best) {
          coefs_all_time_best = coefs
          negloglik_all_time_best = negloglik_curr
        }
        
        if (!is.finite(negloglik_curr)) next
        
        # Least squares can be fitted directly using QR decomposition
        # See https://math.stackexchange.com/questions/852327/efficient-algorithm-for-iteratively-reweighted-least-squares-problem
        # See https://doi.org/10.1093/nar/gks042: 
        #     "logdet is the sum of the logarithms of the diagonal elements of the Cholesky factor R."
        # qr_result = qr(sqrt(weights) * adjusted_design_matrix[, is_active])
        # coefs[-length(coefs)][is_active] = backsolve(qr.R(qr_result),
        #                                              crossprod(qr.Q(qr_result), sqrt(weights) * adjusted_response))
        
        inner_iter = inner_iter + 1
        
      }
      
      # throw a warning of maxiter reached
      if (inner_iter == maxiter) {
        message(paste("Max number of iterations", maxiter, "has been reached."))
      }
      
      number_of_successful_iter = number_of_successful_iter + 1
      
      if (number_of_successful_iter > number_of_resample) break
      
      negloglik_best
    # }, error = function(e) {warning(e); warning("Iteration has failed."); NA}, warning = function(w) {warning(w); NA})
    }, error = function(e) {warning(e); warning("Iteration has failed."); NA})
    
    # We will generate initial values from glm next iteration:
    init = NULL
    if (verbose) cat("\n")
    
  }
  
  # throw error if we do not have number_of_resample successful tries:
  # NOTE: should we add it back?
  # if (number_of_successful_iter != number_of_resample) stop("Optimization using different initial values has failed too many times!")
  
  # The best solution after using multiple initial values:
  coefs = coefs_all_time_best
  if (!use_log_scale_algorithm) {
    coefs[seq_len(H * M) + K] = log(coefs[seq_len(H * M) + K])
  }
    
  # update overdispersion parameter
  if (!fix_overdispersion) {
    
    # calculate likelihood without updated overdispersion to retrieve mu_matrix
    objective_old_overdispersion = negloglik(coefs = coefs,
                                             cell_type_specific_variables = cell_type_specific_variables,
                                             other_variables = other_variables,
                                             read_depth = read_depth,
                                             cellular_proportions = cellular_proportions,
                                             counts = counts,
                                             use_log_scale_params = TRUE,
                                             verbose = TRUE)
    # one dimensional optimization to obtain overdispersion parameter
    optimize_theta = stats::optim(
      par = log(coefs[length(coefs)]),
      negloglik_adjusted_profile,
      mu_matrix = attr(objective_old_overdispersion, "cell.type.specific.fitted.values"),
      coefs = coefs,
      cell_type_specific_variables = cell_type_specific_variables,
      other_variables = other_variables,
      read_depth = read_depth,
      cellular_proportions = cellular_proportions,
      counts = counts,
      is_active = is_active,
      use_log_scale_params = TRUE,
      method = "Brent",
      lower = log(coefs[length(coefs)]) - 1,
      upper = log(coefs[length(coefs)]) + 1
    )
    coefs[length(coefs)] = exp(optimize_theta$par)
  }
  
  # finally, a likelihood without using adjusted profile likelihood
  objective = negloglik(coefs = coefs,
                        cell_type_specific_variables = cell_type_specific_variables,
                        other_variables = other_variables,
                        read_depth = read_depth,
                        cellular_proportions = cellular_proportions,
                        counts = counts,
                        use_log_scale_params = TRUE,
                        verbose = TRUE)
  
  # SEE
  # First, we calculate the hessian matrix one last time:
  # split coefs into beta (cell type-independent) and gamma (cell type-specific)
  beta = coefs[seq_len(K)]
  gamma = matrix(coefs[seq_len(H * M) + K], nrow = H, ncol = M)
  covariate_adjusted_read_depth = exp(other_variables %*% matrix(beta, ncol=1, nrow=K)) * read_depth  # n_B x 1
  mu_matrix = matrix(NA, n_B, H)
  for (i in seq_len(n_B)) {
    # read depth, other effects, cell type-specific effects
    mu_matrix[i, ] = cellular_proportions[i, ] * exp(rowSums(gamma * cell_type_specific_variables[i, , ]))
  }
  mu_matrix = mu_matrix * as.numeric(covariate_adjusted_read_depth)
  mu = rowSums(mu_matrix)

  # in weighted least squares:
  weights = overdispersion_theta / (mu * (mu+overdispersion_theta))
  adjusted_design_matrix = cbind(other_variables * mu,
                                 matrix(cell_type_specific_variables, nrow=n_B) * as.numeric(mu_matrix))
  # Sometimes a column is all 0;
  # sometimes parameters are at the boundary.
  # Need to remove them from active covariates to ensure tXWX is invertible.
  is_active = is_active &
    !apply(adjusted_design_matrix, 2, function(x) isTRUE(all.equal(sum(abs(x)), 0))) &
    ((coefs[-length(coefs)] - log_lower)^2 >= .Machine$double.eps) &
    ((coefs[-length(coefs)] - log_upper)^2 >= .Machine$double.eps)
    
  tXWX = crossprod(adjusted_design_matrix[, is_active], weights * adjusted_design_matrix[, is_active])
  # SEE = sqrt(diag(solve(tXWX)))
  
  # build the matrix of coefficients
  coef_matrix = matrix(NA, nrow=length(coefs), ncol=2)
  colnames(coef_matrix) = c("Estimate", "Std. Error")
  rownames(coef_matrix) = names(coefs)
  coef_matrix[, 1] = coefs
  coef_matrix[, 2][-length(coefs)][is_active] = sqrt(diag(solve(tXWX)))
  
  # These coefficients are actually not in the (reduced) model are set to NAs:
  coef_matrix[K + which(is_reduced), ] = NA
  # The standard error of coeffcients at the boundary at set to NAs:
  coef_matrix[coef_matrix[, 1] == log_upper | coef_matrix[, 1] == log_lower, 2] = NA
  # we do not want to output boundaries, so we replace them with -Inf instead:
  coef_matrix[coef_matrix[, 1] == log_upper, 1] = Inf  # unlikely
  coef_matrix[coef_matrix[, 1] == log_lower, 1] = -Inf
  
  # build lfc matrix & z values
  lfc_matrix = NULL
  # calculate logFoldChange and lfcSE
  if (M >= 2) {
    lfc_matrix = matrix(NA, nrow=(M-1)*H, ncol=4)
    colnames(lfc_matrix) = c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    # assume we use group 1 of 1..M as the baseline:
    is_active_matrix_form = matrix(is_active[-seq_len(K)], nrow=H)
    add_NA_status = function(x, NA_status) {x[NA_status] = NA; x}
    for (m in 2:M) {
      NA_status = !(is_active_matrix_form[, m] & is_active_matrix_form[, 1])
      lfc_matrix_current_rows = (H * (m-2) + 1):(H * (m-1))
      lfc_matrix[lfc_matrix_current_rows, 1] = 
          coef_matrix[(K + H * (m-1) + 1):(K + H * m), 1] - coef_matrix[(K + 1):(K + H), 1]
      # contrast matrix
      R = matrix(0, nrow = H, ncol = K + H * M)
      for (h in seq_len(H)) {
        R[h, c(K + H * (m-1) + h, K + h)] = c(1, -1)
      }
      R = R[, is_active , drop=FALSE]
      # equivalent to sqrt(diag(R %*% solve(tXWX) %*% t(R)))
      lfc_matrix[lfc_matrix_current_rows, 2] = add_NA_status(
          sqrt(diag(R %*% solve(tXWX, t(R)))), NA_status)
    }
    # z value
    lfc_matrix[, 3] = lfc_matrix[, 1] / lfc_matrix[, 2]
    # Pr(>|z|)
    lfc_matrix[, 4] = pmin(1, 2 * pnorm(-abs(lfc_matrix[, 3])))
  }
  
  # Calculate Cook's distance
  # Williams, D. A. “Generalized Linear Model Diagnostics Using the Deviance and Single Case Deletions.”
  #     Applied Statistics 36, no. 2 (1987): 181. https://doi.org/10.2307/2347550.
  # C_i = 1/p * h_i/(1 - h_i) * r^2_{Pi}
  p = sum(is_active)
  sqrtWX = sqrt(weights) * adjusted_design_matrix[, is_active]
  hat_values = diag(sqrtWX %*% solve(tXWX, t(sqrtWX)))
  # variance of negative binomial is 1/weights
  # standardised Pearson residuals
  r_Pi = (counts - mu) / sqrt((1 - hat_values) / weights)
  cooks = 1/p * hat_values / (1 - hat_values) * r_Pi^2
  
  # Add human friendly names to coefficients
  # dimnames(cell_type_specific_variables)[[1]] = NULL
  if (is.null(dimnames(cell_type_specific_variables)[[2]])) {
    dimnames(cell_type_specific_variables)[[2]] = paste0("celltype", seq_len(H))
  }
  if (is.null(dimnames(cell_type_specific_variables)[[3]])) {
    dimnames(cell_type_specific_variables)[[3]] = paste0("variable", seq_len(M))
  }
  colnames_other_variables = NULL
  if (!is.null(other_variables) && ncol(other_variables) > 0) {
    if (is.null(colnames(other_variables))) {
      colnames_other_variables = paste0("other_variables", seq_len(K))
    } else {
      colnames_other_variables = colnames(other_variables)
    }
  }
  
  # Cell type-specific effects as interaction terms
  rownames(coef_matrix)[-length(coefs)] = 
      c(colnames_other_variables, 
        paste(rep(dimnames(cell_type_specific_variables)[[3]], each=H),
              dimnames(cell_type_specific_variables)[[2]],
              sep=":"))
  if (M >= 2) {
    rownames(lfc_matrix) = 
        sprintf("%s_vs_%s:%s",
                rep(dimnames(cell_type_specific_variables)[[3]][-1], each=H),
                dimnames(cell_type_specific_variables)[[3]][1],
                dimnames(cell_type_specific_variables)[[2]])
  }
  
  # return a list of coefficients and negative log-likelihood
  list(coefficients = coef_matrix,
       lfc = lfc_matrix,
       value = as.numeric(objective),
       cooks.distance = cooks,
       fitted.values = mu
       # ,
       # cell.type.specific.fitted.values = attr(objective, "cell.type.specific.fitted.values"),
       # negloglik_curr_vector = negloglik_curr_vector
       )
}






#' negative log-likelihood function of a CARseq negative binomial model
#' 
#' \code{negloglik} provides a negative log-likelihood function 
#' of negative binomial distribution 
#' whose mean is a sum of nonnegative terms with covariates 
#'
#' @param coefs a vector of c(beta, gamma, overdispersion), where beta is a length K vector of cell type-indepedent coefficients,
#'   gamma is a matrix of H x M dimension of cell type-specific non-negative coefficients, and overdispersion is a scalar.
#' @param cell_type_specific_variables an array of n_B x H x K of cell type-independent variables.
#' @param other_variables a design matrix of n_B x M of cell type-specific variables. Can be NULL.
#' @param read_depth a vector of n_B sample-specific read depths.
#' @param cellular_proportions a matrix of n_B x H of cellular proportions.
#' @param counts A vector of n_B total read counts observed.
#' @param use_log_scale_params logical. The default is FALSE, where the cell type-specific effects are
#'        parametrized on a non-log scale. Otherwise, they are parametrized on a log scale.
#' @param verbose logical. If \code{TRUE}, a matrix fitted cell type-specific expression will be attached.
#' 
#' @examples
#' library(CARseq)
#' set.seed(1234)
#' H = 4
#' n_B = 60
#' K = 1
#' M = 3
#' overdispersion_theta = 10
#' coefs = c(rep(0, K), rep(1, H*M), overdispersion_theta)  # initial value
#' cell_type_specific_variables = array(0, dim=c(n_B, H, M))
#' cell_type_specific_variables[1:(n_B/3), , 1] = 1
#' cell_type_specific_variables[(n_B/3+1):(n_B*2/3), , 2] = 1
#' cell_type_specific_variables[(n_B*2/3+1):n_B, , 3] = 1
#' other_variables = matrix(runif(n_B * K, min = -1, max = 1), nrow=n_B, ncol=K)
#' read_depth = round(runif(n_B, min = 50, max = 100))
#' cellular_proportions = matrix(runif(n_B * H), nrow = n_B, ncol = H)
#' cellular_proportions = cellular_proportions / rowSums(cellular_proportions)
#' counts = round(runif(n_B, min = 100, max = 200))
#' CARseq:::negloglik(coefs, cell_type_specific_variables, other_variables, read_depth, cellular_proportions, counts)
negloglik = function(coefs, cell_type_specific_variables, other_variables, read_depth, cellular_proportions, counts,
                     use_log_scale_params = FALSE, verbose = FALSE) {
  # obtain H, M, K, n_B
  n_B = length(counts)
  if (is.null(other_variables)) other_variables = matrix(ncol = 0, nrow = n_B)
  if (!is.array(cell_type_specific_variables)) {
    stop("cell_type_specific_variables should be an array")
  } else if (length(dim(cell_type_specific_variables)) != 3) {
    stop("cell_type_specific_variables should be a 3-dimensional (n_B, H, M) array")
  }
  if (!is.matrix(other_variables)) stop("other_variables should be a matrix")
  if (!is.matrix(cellular_proportions)) stop("cellular_proportions should be a matrix")
  if (dim(cell_type_specific_variables)[1] != n_B) stop("Length of the 1st dimension in cell_type_specific_variables does not match the length of counts")
  if (nrow(other_variables) != n_B) stop("Number of rows in other_variables does not match the length of counts")
  if (nrow(cellular_proportions) != n_B) stop("Number of rows in cellular_proportions does not match the length of counts")
  M = dim(cell_type_specific_variables)[3]
  K = ncol(other_variables)
  H = ncol(cellular_proportions)
  if (length(coefs) != K + H * M + 1) stop("The length of coefs is not equal to K + H * M + 1")
  
  # split coefs into beta (cell type-independent) and gamma (cell type-specific)
  beta = coefs[seq_len(K)]
  gamma = matrix(coefs[seq_len(H * M) + K], nrow = H, ncol = M)
  overdispersion = coefs[length(coefs)]
  
  covariate_adjusted_read_depth = exp(other_variables %*% matrix(beta, ncol=1, nrow=K)) * read_depth  # n_B x 1
  
  mu_matrix = matrix(NA, n_B, H)
  if (use_log_scale_params) {
    for (i in seq_len(n_B)) {
      # read depth, other effects, cell type-specific effects
      mu_matrix[i, ] = cellular_proportions[i, ] * exp(rowSums(gamma * cell_type_specific_variables[i, , ]))
    }
  } else {
    for (i in seq_len(n_B)) {
      mu_matrix[i, ] = cellular_proportions[i, ] * matrixStats::rowProds(gamma ^ cell_type_specific_variables[i, , ])
    }
  }
  mu_matrix = mu_matrix * as.numeric(covariate_adjusted_read_depth)
  mu = rowSums(mu_matrix)
  objective = - sum(lgamma(counts + overdispersion) - lgamma(overdispersion) - lgamma(counts + 1) +
                      overdispersion*log(overdispersion) + counts*log(mu) -
                      (counts + overdispersion) * log(overdispersion+mu))
  if (verbose)
    attr(objective, "cell.type.specific.fitted.values") = mu_matrix
  objective
}


#' negative adjusted profile log-likelihood function of a CARseq negative binomial model
#' 
#' \code{negloglik_adjusted_profile} provides a negative
#' adjusted profile log-likelihood function 
#' of negative binomial distribution 
#' whose mean is a sum of nonnegative terms with covariates.
#' This is used to estimate the overdispersion parameter.
#'
#' @param logtheta log(overdispersion) as a scalar. Overdispersion is parametrized as \code{theta} in \link[MASS]{glm.nb}.
#' @param mu_matrix cell type-specific fitted values from CARseq negative binomial model.
#' @param coefs a vector of c(beta, gamma, overdispersion), where beta is a length K vector of cell type-indepedent coefficients,
#'   gamma is a matrix of H x M dimension of cell type-specific non-negative coefficients, and overdispersion is a scalar.
#'   Note: only gamma is used in the evaluation here.
#' @param cell_type_specific_variables an array of n_B x H x K of cell type-independent variables.
#' @param other_variables a design matrix of n_B x M of cell type-specific variables. Can be NULL.
#' @param read_depth a vector of n_B sample-specific read depths.
#' @param cellular_proportions a matrix of n_B x H of cellular proportions.
#' @param counts A vector of n_B total read counts observed.
#' @param is_active A logical vector recording which predictors are included in the model.
#' @param use_log_scale_params logical. The default is FALSE, where the cell type-specific effects are
#'        parametrized on a non-log scale. Otherwise, they are parametrized on a log scale.
#' 
#' @examples
#' library(CARseq)
#' set.seed(1234)
#' H = 4
#' n_B = 60
#' K = 1
#' M = 3
#' overdispersion_theta = 10
#' coefs = c(rep(0, K), rep(1, H*M), overdispersion_theta)  # initial value
#' cell_type_specific_variables = array(0, dim=c(n_B, H, M))
#' cell_type_specific_variables[1:(n_B/3), , 1] = 1
#' cell_type_specific_variables[(n_B/3+1):(n_B*2/3), , 2] = 1
#' cell_type_specific_variables[(n_B*2/3+1):n_B, , 3] = 1
#' other_variables = matrix(runif(n_B * K, min = -1, max = 1), nrow=n_B, ncol=K)
#' read_depth = round(runif(n_B, min = 50, max = 100))
#' cellular_proportions = matrix(runif(n_B * H), nrow = n_B, ncol = H)
#' cellular_proportions = cellular_proportions / rowSums(cellular_proportions)
#' counts = round(runif(n_B, min = 100, max = 200))
#' is_active = rep(TRUE, K + H*M)
#' res = CARseq:::negloglik(coefs, cell_type_specific_variables, other_variables, read_depth, cellular_proportions, counts, verbose=TRUE)
#' CARseq:::negloglik_adjusted_profile(log(overdispersion_theta), attr(res, "cell.type.specific.fitted.values"), coefs, cell_type_specific_variables, other_variables, read_depth, cellular_proportions, counts, is_active)
#' @references Davis J. McCarthy, Yunshun Chen, Gordon K. Smyth, Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation, \emph{Nucleic Acids Research}, Volume 40, Issue 10, 1 May 2012, Pages 4288-4297, https://doi.org/10.1093/nar/gks042
negloglik_adjusted_profile = function(logtheta,
                                      mu_matrix,
                                      coefs,
                                      cell_type_specific_variables,
                                      other_variables,
                                      read_depth,
                                      cellular_proportions,
                                      counts, 
                                      is_active = NULL,
                                      use_log_scale_params = FALSE) {

  overdispersion = exp(logtheta)
  mu = rowSums(mu_matrix)
  
  # obtain H, M, K, n_B
  n_B = length(counts)
  if (is.null(other_variables)) other_variables = matrix(ncol = 0, nrow = n_B)
  if (!is.array(cell_type_specific_variables)) {
    stop("cell_type_specific_variables should be an array")
  } else if (length(dim(cell_type_specific_variables)) != 3) {
    stop("cell_type_specific_variables should be a 3-dimensional (n_B, H, M) array")
  }
  if (!is.matrix(other_variables)) stop("other_variables should be a matrix")
  if (!is.matrix(cellular_proportions)) stop("cellular_proportions should be a matrix")
  if (dim(cell_type_specific_variables)[1] != n_B) stop("Length of the 1st dimension in cell_type_specific_variables does not match the length of counts")
  if (nrow(other_variables) != n_B) stop("Number of rows in other_variables does not match the length of counts")
  if (nrow(cellular_proportions) != n_B) stop("Number of rows in cellular_proportions does not match the length of counts")
  M = dim(cell_type_specific_variables)[3]
  K = ncol(other_variables)
  H = ncol(cellular_proportions)
  if (length(coefs) != K + H * M + 1) stop("The length of coefs is not equal to K + H * M + 1")
  
  # split coefs into beta (cell type-independent) and gamma (cell type-specific)
  beta = coefs[seq_len(K)]   # not used later
  gamma = matrix(coefs[seq_len(H * M) + K], nrow = H, ncol = M)
  
  # -1/2 logdet
  weights = overdispersion / (mu * (mu+overdispersion))
  if (use_log_scale_params) {
    adjusted_design_matrix = cbind(other_variables * mu,
        matrix(cell_type_specific_variables, nrow=n_B) * as.numeric(mu_matrix))
  } else {
    adjusted_design_matrix = cbind(other_variables * mu,
        matrix(cell_type_specific_variables, nrow=n_B) / rep(gamma, each=n_B) * as.numeric(mu_matrix))
  }
  qr_result = qr(sqrt(weights) * adjusted_design_matrix[, is_active])
  logdet = sum(log(abs(diag(qr.R(qr_result)))))
  
  - sum(lgamma(counts + overdispersion) - lgamma(overdispersion) - lgamma(counts + 1) +
          overdispersion*log(overdispersion) + counts*log(mu) -
          (counts + overdispersion) * log(overdispersion+mu)) +
    0.5 * logdet
}


#' Gradient of negative log-likelihood of a negative binomial model
#' 
#' \code{grad_negloglik} provides gradient of negative log-likelihood function
#'  of negative binomial distribution 
#' whose mean is a sum of nonnegative terms with covariates.
#'
#' The cell type-specific effects are parametrized on a non-log scale.
#'
#' @param coefs a vector of c(beta, gamma, overdispersion), where beta is a length K vector of cell type-indepedent coefficients,
#'   gamma is a matrix of H x M dimension of cell type-specific non-negative coefficients, and overdispersion is a scalar.
#' @param cell_type_specific_variables a design matrix of n_B x K of cell type-independent variables.
#' @param other_variables a design matrix of n_B x M of cell type-specific variables. Can be NULL.
#' @param read_depth a vector of n_B sample-specific read depths.
#' @param cellular_proportions a matrix of n_B x H of cellular proportions.
#' @param counts A vector of n_B total read counts observed.
#' 
#' @examples
#' library(CARseq)
#' set.seed(1234)
#' H = 4
#' n_B = 60
#' K = 1
#' M = 3
#' coefs = c(rep(0, K), rep(1, H*M), overdispersion=10)  # initial value
#' cell_type_specific_variables = array(0, dim=c(n_B, H, M))
#' cell_type_specific_variables[1:(n_B/3), , 1] = 1
#' cell_type_specific_variables[(n_B/3+1):(n_B*2/3), , 2] = 1
#' cell_type_specific_variables[(n_B*2/3+1):n_B, , 3] = 1
#' other_variables = matrix(runif(n_B * K, min = -1, max = 1), nrow=n_B, ncol=K)
#' read_depth = round(runif(n_B, min = 50, max = 100))
#' cellular_proportions = matrix(runif(n_B * H), nrow = n_B, ncol = H)
#' cellular_proportions = cellular_proportions / rowSums(cellular_proportions)
#' counts = round(runif(n_B, min = 100, max = 200))
#' CARseq:::grad_negloglik(coefs, cell_type_specific_variables, other_variables, read_depth, cellular_proportions, counts)
grad_negloglik = function(coefs, cell_type_specific_variables, other_variables, read_depth, cellular_proportions, counts) {
  # obtain H, M, K, n_B
  n_B = length(counts)
  if (is.null(other_variables)) other_variables = matrix(ncol = 0, nrow = n_B)
  if (!is.array(cell_type_specific_variables)) {
    stop("cell_type_specific_variables should be an array")
  } else if (length(dim(cell_type_specific_variables)) != 3) {
    stop("cell_type_specific_variables should be a 3-dimensional (n_B, H, M) array")
  }
  if (!is.matrix(other_variables)) stop("other_variables should be a matrix")
  if (!is.matrix(cellular_proportions)) stop("cellular_proportions should be a matrix")
  if (dim(cell_type_specific_variables)[1] != n_B) stop("Length of the 1st dimension in cell_type_specific_variables does not match the length of counts")
  if (nrow(other_variables) != n_B) stop("Number of rows in other_variables does not match the length of counts")
  if (nrow(cellular_proportions) != n_B) stop("Number of rows in cellular_proportions does not match the length of counts")
  M = dim(cell_type_specific_variables)[3]
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
  # TODO: The following loop takes half of total time. Need to use Rcpp to make it faster:
  for (i in seq_len(n_B)) {
    # read depth, other effects, cell type-specific effects
    mu[i] = sum(cellular_proportions[i, ] * matrixStats::rowProds(gamma ^ cell_type_specific_variables[i, , ]))
  }
  mu = mu * covariate_adjusted_read_depth
  
  frac_difference = counts / mu - (counts + overdispersion) / (overdispersion + mu)
  
  # derivative wrt cell type-independent variables
  partial_l_partial_beta_i = matrix(NA, nrow = n_B, ncol = K)
  for (i in seq_len(n_B)) {
    partial_l_partial_beta_i[i, ] = frac_difference[i] * other_variables[i, ] * covariate_adjusted_read_depth[i] *
      sum(cellular_proportions[i, ] * matrixStats::rowProds(gamma ^ cell_type_specific_variables[i, , ]))
  }
  g[seq_len(K)] = colSums(partial_l_partial_beta_i)
  
  # derivative wrt cell type-specific variables
  partial_l_partial_gamma = matrix(NA, nrow = H, ncol = M)
  for (m in seq_len(M)) {
    partial_l_partial_gamma_i = matrix(NA, nrow = n_B, ncol = H)
    
    for (i in seq_len(n_B)) {
      cell_type_specific_variables_i = cell_type_specific_variables[i, , ]
      cell_type_specific_variables_i[, m] = cell_type_specific_variables_i[, m] - 1
      partial_l_partial_gamma_i[i, ] = frac_difference[i] * covariate_adjusted_read_depth[i] *
        cellular_proportions[i, ] * cell_type_specific_variables[i, , m] * 
        matrixStats::rowProds(gamma ^ cell_type_specific_variables_i)
    }
    g[seq(H) + K + H * (m-1)] = colSums(partial_l_partial_gamma_i)
  }
  g[K + H * M + 1] = sum(digamma(counts + overdispersion) - digamma(overdispersion) +
                           log(overdispersion) + 1.0 - (counts + overdispersion) / (overdispersion + mu) -
                           log(overdispersion + mu))
  - g
}



#' negative log-likelihood function of a CARseq negative binomial model
#' 
#' \code{negloglik_legacy} provides a negative log-likelihood function 
#' of negative binomial distribution 
#' whose mean is a sum of nonnegative terms with covariates 
#'
#' @param coefs a vector of c(beta, gamma, overdispersion), where beta is a length K vector of cell type-indepedent coefficients,
#'   gamma is a matrix of H x M dimension of cell type-specific non-negative coefficients, and overdispersion is a scalar.
#' @param cell_type_specific_variables an array of n_B x H x K of cell type-independent variables.
#' @param other_variables a design matrix of n_B x M of cell type-specific variables. Can be NULL.
#' @param read_depth a vector of n_B sample-specific read depths.
#' @param cellular_proportions a matrix of n_B x H of cellular proportions.
#' @param counts A vector of n_B total read counts observed.
#' 
#' @examples
#' library(CARseq)
#' set.seed(1234)
#' H = 4
#' n_B = 60
#' K = 1
#' M = 3
#' overdispersion_theta = 10
#' coefs = c(rep(0, K), rep(1, H*M), overdispersion_theta)  # initial value
#' cell_type_specific_variables = array(0, dim=c(n_B, H, M))
#' cell_type_specific_variables[1:(n_B/3), , 1] = 1
#' cell_type_specific_variables[(n_B/3+1):(n_B*2/3), , 2] = 1
#' cell_type_specific_variables[(n_B*2/3+1):n_B, , 3] = 1
#' other_variables = matrix(runif(n_B * K, min = -1, max = 1), nrow=n_B, ncol=K)
#' read_depth = round(runif(n_B, min = 50, max = 100))
#' cellular_proportions = matrix(runif(n_B * H), nrow = n_B, ncol = H)
#' cellular_proportions = cellular_proportions / rowSums(cellular_proportions)
#' counts = round(runif(n_B, min = 100, max = 200))
#' CARseq:::negloglik_legacy(coefs, cell_type_specific_variables, other_variables, read_depth, cellular_proportions, counts)
#' @export
negloglik_legacy = function(coefs, cell_type_specific_variables, other_variables, read_depth, cellular_proportions, counts) {
  # obtain H, M, K, n_B
  n_B = length(counts)
  if (is.null(other_variables)) other_variables = matrix(ncol = 0, nrow = n_B)
  if (!is.array(cell_type_specific_variables)) {
    stop("cell_type_specific_variables should be an array")
  } else if (length(dim(cell_type_specific_variables)) != 3) {
    stop("cell_type_specific_variables should be a 3-dimensional (n_B, H, M) array")
  }
  if (!is.matrix(other_variables)) stop("other_variables should be a matrix")
  if (!is.matrix(cellular_proportions)) stop("cellular_proportions should be a matrix")
  if (dim(cell_type_specific_variables)[1] != n_B) stop("Length of the 1st dimension in cell_type_specific_variables does not match the length of counts")
  if (nrow(other_variables) != n_B) stop("Number of rows in other_variables does not match the length of counts")
  if (nrow(cellular_proportions) != n_B) stop("Number of rows in cellular_proportions does not match the length of counts")
  M = dim(cell_type_specific_variables)[3]
  K = ncol(other_variables)
  H = ncol(cellular_proportions)
  if (length(coefs) != K + H * M + 1) stop("The length of coefs is not equal to K + H * M + 1")
  
  # split coefs into beta (cell type-independent) and gamma (cell type-specific)
  beta = coefs[seq_len(K)]
  gamma = matrix(coefs[seq_len(H * M) + K], nrow = H, ncol = M)
  overdispersion = coefs[length(coefs)]
  
  covariate_adjusted_read_depth = exp(other_variables %*% matrix(beta, ncol=1, nrow=K)) * read_depth  # n_B x 1

  mu_matrix = matrix(NA, n_B, H)
  for (i in seq_len(n_B)) {
    # read depth, other effects, cell type-specific effects
    mu_matrix[i, ] = cellular_proportions[i, ] * matrixStats::rowProds(gamma ^ cell_type_specific_variables[i, , ])
  }
  mu_matrix = mu_matrix * as.numeric(covariate_adjusted_read_depth)
  mu = rowSums(mu_matrix)
  objective = - sum(lgamma(counts + overdispersion) - lgamma(overdispersion) - lgamma(counts + 1) +
        overdispersion*log(overdispersion) + counts*log(mu) -
        (counts + overdispersion) * log(overdispersion+mu))
  attr(objective, "cell.type.specific.fitted.values") = mu_matrix
  objective
}


#' fit a CARseq negative binomial model
#' 
#' \code{fit_model_legacy} fits a negative binomial distribution 
#' whose mean is a sum of nonnegative terms with covariates.
#' The overdispersion parameter is estimated by maximizing the 
#' adjusted profile log-likelihood.
#' 
#' This is an old version where the result is highly dependent on initial values.
#' Please use \code{\link{fit_model}} instead. The arguments are a little different.
#'
#' @param init a numeric vector of (K + H x M + 1) corresponding to the initial value of
#'             c(beta, gamma, overdispersion). Can be NULL.
#' @param cell_type_specific_variables an array of n_B x H x K of cell type-independent variables.
#' @param other_variables a design matrix of n_B x M of cell type-specific variables. Can be NULL.
#' @param read_depth a vector of n_B sample-specific read depths.
#' @param cellular_proportions a matrix of n_B x H of cellular proportions.
#' @param counts A vector of n_B total read counts observed.
#' @param is_active A logical vector of length (K + H x M), recording which predictors are included in the model.
#' @param fix_overdispersion logical. If \code{TRUE}, the overdispersion parameter will not be updated.
#' @param number_of_resample numeric. 
#'
#' @examples
#' library(CARseq)
#' set.seed(1234)
#' H = 4
#' n_B = 60
#' K = 1
#' M = 3
#' cell_type_specific_variables = array(0, dim=c(n_B, H, M))
#' cell_type_specific_variables[1:(n_B/3), , 1] = 1
#' cell_type_specific_variables[(n_B/3+1):(n_B*2/3), , 2] = 1
#' cell_type_specific_variables[(n_B*2/3+1):n_B, , 3] = 1
#' other_variables = matrix(runif(n_B * K, min = -1, max = 1), nrow=n_B, ncol=K)
#' read_depth = round(runif(n_B, min = 50, max = 100))
#' cellular_proportions = matrix(runif(n_B * H), nrow = n_B, ncol = H)
#' cellular_proportions = cellular_proportions / rowSums(cellular_proportions)
#' counts = round(runif(n_B, min = 100, max = 200))
#' is_active = rep(TRUE, K + H*M)
#' res = CARseq:::fit_model_legacy(init = NULL, cell_type_specific_variables, other_variables, read_depth, cellular_proportions, counts, is_active)
#' str(res)
fit_model_legacy = function(init = NULL, cell_type_specific_variables, other_variables, read_depth, cellular_proportions, counts,
                     is_active = NULL, fix_overdispersion = FALSE, number_of_resample = 20) {
  
  # obtain H, M, K, n_B
  n_B = length(counts)
  if (is.null(other_variables)) other_variables = matrix(ncol = 0, nrow = n_B)
  if (!is.array(cell_type_specific_variables)) {
    stop("cell_type_specific_variables should be an array")
  } else if (length(dim(cell_type_specific_variables)) != 3) {
    stop("cell_type_specific_variables should be a 3-dimensional (n_B, H, M) array")
  }
  if (!is.matrix(other_variables)) stop("other_variables should be a matrix")
  if (!is.matrix(cellular_proportions)) stop("cellular_proportions should be a matrix")
  if (dim(cell_type_specific_variables)[1] != n_B) stop("Length of the 1st dimension in cell_type_specific_variables does not match the length of counts")
  if (nrow(other_variables) != n_B) stop("Number of rows in other_variables does not match the length of counts")
  if (nrow(cellular_proportions) != n_B) stop("Number of rows in cellular_proportions does not match the length of counts")
  M = dim(cell_type_specific_variables)[3]
  K = ncol(other_variables)
  H = ncol(cellular_proportions)
  coefs = init
  if (is.null(is_active)) is_active = rep(TRUE, K + H * M)
  
  # We need to sample many times to decide on an initial value.
  if (is.null(number_of_resample) || is.na(number_of_resample)) number_of_resample = 20
  # We only need number_of_resample different initial values and subsequent model fits.
  # However, we can try at most number_of_resample_max times.
  number_of_resample_max = number_of_resample * 2
  resample_size_list = rep(3 * (K + H * M), number_of_resample_max)    # An arbitrary number that should be larger than the number of parameters in the model
  resample_size_list[1] = n_B    # in the first instance, the initial value is based on all samples instead of a small subset of samples
  
  if (fix_overdispersion) init_overdispersion = coefs[length(coefs)]
  
  negloglik_best = NULL
  coefs_best = NULL
  number_of_successful_iter = 0
  
  negloglik_curr_vector = rep(NA, number_of_resample_max)  # stores a vector of likelihood for debugging purposes
  for (iter in seq_len(number_of_resample_max)) {
    resample_size = resample_size_list[iter]
    
    # save the current likelihood for each optimization path with a new initial value
    negloglik_curr_vector[iter] =  tryCatch({
    
      # # generate initial value using glm and a subset of samples. Similar to "bootstrap root search" in the literature.
      # if (K == 0) {
      #   nbmodel = MASS::glm.nb(counts~offset(log(read_depth)),
      #                          subset = sample(c(rep(TRUE, resample_size), rep(FALSE, n_B - resample_size))))
      # } else {
      #   nbmodel = MASS::glm.nb(counts~offset(log(read_depth))+other_variables, 
      #                          subset = sample(c(rep(TRUE, resample_size), rep(FALSE, n_B - resample_size))))
      # }
      # 
      # if (is.null(init)) {
      #   # initialize using MASS::glm.nb
      #   coefs = c(nbmodel$coefficients[-1],                         # cell type-independent variables
      #             rep(exp(nbmodel$coefficients[1]), H * M),         # cell type-specific variables
      #             overdispersion = nbmodel$theta)
      #   coefs[is.na(coefs)] = 0
      
      # generate initial value using glm and a subset of samples. Similar to "bootstrap root search" in the literature.
      
      array_inds_sampled_by_cell_type = as.matrix(merge(cbind(1:n_B, sample(x = 1:H, size = n_B, replace = TRUE)), seq_len(M), by=NULL))
      variables_sampled_cell_type_specific = matrix(cell_type_specific_variables[array_inds_sampled_by_cell_type], nrow=n_B, ncol=M)
      if (K == 0) {
        nbmodel = MASS::glm.nb(counts~offset(log(read_depth)) + 0 + variables_sampled_cell_type_specific,
                               subset = sample(c(rep(TRUE, resample_size), rep(FALSE, n_B - resample_size))))
      } else {
        nbmodel = MASS::glm.nb(counts~offset(log(read_depth)) + 0 + other_variables + variables_sampled_cell_type_specific, 
                               subset = sample(c(rep(TRUE, resample_size), rep(FALSE, n_B - resample_size))))
      }
     
      if (is.null(init)) {
        # initialize using MASS::glm.nb
        coefs = c(nbmodel$coefficients[seq_len(K)],                     # cell type-independent variables
                  rep(exp(nbmodel$coefficients[-seq_len(K)]), each=H),   # cell type-specific variables
                  overdispersion = nbmodel$theta)
        coefs[is.na(coefs)] = 0
        
        # TODO: use cooks.distance() and fitted() to process the outliers
      } else if (length(coefs) != K + H * M + 1) {
        stop("The length of coefs is not equal to K + H * M + 1")
      }
      
      if (fix_overdispersion) {
        coefs[length(coefs)] = init_overdispersion
      }
      
      # some practical bounds to make the method work numerically:
      maxcoef_log = 50
      maxmu = 1e6
      minmu = 1e-6
      
      epsilon_convergence = 1e-3
      maxiter = 20
      
      overdispersion_theta = coefs[length(coefs)]
      
      # variables to control the loop:
      # the negative log-likelihood needs to decrease.
      # break the loop if it does not decrease a lot, or it increases:
      inner_iter = 0
      negloglik_prev = negloglik_curr = Inf
      coefs_prev = coefs
      
      while (inner_iter == 0 ||
           (abs(negloglik_prev - negloglik_curr) > epsilon_convergence
            &&
            inner_iter < maxiter)) {
      
        # split coefs into beta (cell type-independent) and gamma (cell type-specific)
        beta = coefs[seq_len(K)]
        gamma = matrix(coefs[seq_len(H * M) + K], nrow = H, ncol = M)
        covariate_adjusted_read_depth = exp(other_variables %*% matrix(beta, ncol=1, nrow=K)) * read_depth  # n_B x 1
        mu_matrix = matrix(NA, n_B, H)
        for (i in seq_len(n_B)) {
          # read depth, other effects, cell type-specific effects
          mu_matrix[i, ] = cellular_proportions[i, ] * matrixStats::rowProds(gamma ^ cell_type_specific_variables[i, , ])
        }
        mu_matrix = mu_matrix * as.numeric(covariate_adjusted_read_depth)
        mu = rowSums(mu_matrix)
        
        # current negative log-likelihood:
        negloglik_prev = negloglik_curr
        negloglik_curr = - sum(lgamma(counts + overdispersion_theta) - lgamma(overdispersion_theta) - lgamma(counts + 1) +
                                 overdispersion_theta*log(overdispersion_theta) + counts*log(mu) -
                                 (counts + overdispersion_theta) * log(overdispersion_theta+mu))
        if (!is.finite(negloglik_curr)) break
        
        # save the best coefs and negative log-likelihood
        if (is.null(negloglik_best) || negloglik_curr < negloglik_best) {
          coefs_best = coefs
          negloglik_best = negloglik_curr   # the matching coefs should be coefs_prev, but it would not matter anyway when converged
        }
        
        # print(negloglik_curr)
        
        if (!is.finite(negloglik_curr)) next
        
        # need to roll back if we find the negative-loglikelihood starts to increase:
        coefs_prev = coefs
        
        # in weighted least squares:
        weights = overdispersion_theta / (mu * (mu+overdispersion_theta))
        adjusted_design_matrix = cbind(other_variables * mu,
                                       matrix(cell_type_specific_variables, nrow=n_B) / rep(gamma, each=n_B) * as.numeric(mu_matrix))
        
        adjusted_response = adjusted_design_matrix %*% coefs[-length(coefs)] + (counts - mu)

        # bounded-variable least squares ensures that the active parameters in H*M cell type-specific ones are bounded
        coefs[-length(coefs)][is_active] = bvls::bvls(adjusted_design_matrix[, is_active] * sqrt(weights),
                                                      adjusted_response * sqrt(weights),
                                                      bl = c(rep(-maxcoef_log, K), rep(minmu, sum(is_active) - K)),
                                                      bu = c(rep( maxcoef_log, K), rep(maxmu, sum(is_active) - K)))$x
        
        # newfit = stats::lsfit(adjusted_design_matrix[, is_active], adjusted_response, wt=weights, intercept=FALSE)$coef
        # ind_min = K + which.min(newfit[K + seq_len(sum(is_active) - K)])
        # # ind_min = K + which.min(newfit[K + which(is_active[K+seq(H*M)]]))
        # if (newfit[ind_min] < minmu) {
        #   lambda = (coefs[is_active][ind_min] - minmu) / (coefs[is_active][ind_min] - newfit[ind_min])
        #   coefs[-length(coefs)][is_active] = coefs[-length(coefs)][is_active] * (1 - lambda) + newfit * lambda
        # } else {
        #   coefs[-length(coefs)][is_active] = newfit
        # }
        
        # lsfit cannot guarantee non-negativity:
        # coefs[-length(coefs)][is_active] = stats::lsfit(adjusted_design_matrix[, is_active], adjusted_response, wt=weights, intercept=FALSE)$coef
        
        # Least squares can be fitted directly using QR decomposition
        # See https://math.stackexchange.com/questions/852327/efficient-algorithm-for-iteratively-reweighted-least-squares-problem
        # See https://doi.org/10.1093/nar/gks042: 
        #     "logdet is the sum of the logarithms of the diagonal elements of the Cholesky factor R."
        # qr_result = qr(sqrt(weights) * adjusted_design_matrix[, is_active])
        # coefs[-length(coefs)][is_active] = backsolve(qr.R(qr_result),
        #                                              crossprod(qr.Q(qr_result), sqrt(weights) * adjusted_response))
        
        inner_iter = inner_iter + 1
      
      }
      
      # throw a warning of maxiter reached
      if (inner_iter == maxiter) {
        message("Max number of iterations has been reached.")
      }
      
      number_of_successful_iter = number_of_successful_iter + 1
      
      if (number_of_successful_iter == number_of_resample) break
      
      negloglik_curr
    
    }, error = function(e) {warning("Iteration has failed.")}, warning = function(w) {warning(w)})
    
    # We will generate initial values from glm next iteration:
    init = NULL
    cat("\n")
    
  }
  
  # throw error if we do not have number_of_resample successful tries:
  # TODO: need to add back
  # if (number_of_successful_iter != number_of_resample) stop("Optimization using different initial values has failed too many times!")
  
  # The best solution after using multiple initial values:
  coefs = coefs_best

  # update overdispersion parameter
  if (!fix_overdispersion) {
    
    # likelihood without updated overdispersion
    objective_old_overdispersion = negloglik(coefs = coefs,
                                             cell_type_specific_variables = cell_type_specific_variables,
                                             other_variables = other_variables,
                                             read_depth = read_depth,
                                             cellular_proportions = cellular_proportions,
                                             counts = counts,
                                             verbose = TRUE)
    
    optimize_theta = stats::optim(
      par = log(coefs[length(coefs)]),
      negloglik_adjusted_profile,
      mu_matrix = attr(objective_old_overdispersion, "cell.type.specific.fitted.values"),
      coefs = coefs,
      cell_type_specific_variables = cell_type_specific_variables,
      other_variables = other_variables,
      read_depth = read_depth,
      cellular_proportions = cellular_proportions,
      counts = counts,
      is_active = is_active,
      method = "L-BFGS-B",
      lower = -5,
      upper = 10
      )
    coefs[length(coefs)] = exp(optimize_theta$par)
  }
  
  # finally, a likelihood without using adjusted profile likelihood
  objective = negloglik(coefs = coefs,
                        cell_type_specific_variables = cell_type_specific_variables,
                        other_variables = other_variables,
                        read_depth = read_depth,
                        cellular_proportions = cellular_proportions,
                        counts = counts,
                        verbose = TRUE)
  
  # SEE
  # tXWX = crossprod(adjusted_design_matrix[, is_active], weights * adjusted_design_matrix[, is_active])
  # rbind(estimates=coefs[-length(coefs)][is_active], SEE=sqrt(diag(solve(tXWX))))
  
  list(par = coefs, 
       value = as.numeric(objective),
       cell.type.specific.fitted.values = attr(objective, "cell.type.specific.fitted.values"),
       negloglik_curr_vector = negloglik_curr_vector)
}