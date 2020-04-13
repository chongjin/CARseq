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
#' @import foreach
#' @import doMC
#' @import matrixStats
#' @import zetadiv
#' @import MASS
#' @import bvls
#' @import nloptr
NULL


#' Conduct a CARseq test
#'
#' @param count_matrix A matrix of G x n total read counts observed.
#' @param cellular_proportions A matrix of n x H of cellular proportions.
#' @param groups A vector of length n indicating groups that we would like to test. Will be coerced to be factors.
#' @param formula A formula of an intercept term "1" and other cell type-indepedent variables. Can be NULL.
#' @param data A data frame containing the cell type-indepedent variables specified in \code{formula}. Can be NULL.
#' @param read_depth A vector of n sample-specific read depths. It it used as an offset term in the
#'        regression model. Alternatively, it can be 1, NA or NULL so that we don't have an offset term
#'        in the model, and log(read_depth) can be included as one of the cell type-independent variables.
#' @param shrunken_lfc Logical. If \code{TRUE} (default), provide shrunken log fold change for cell type-specific variables.
#' @param cores Numeric. Number of cores to use in \code{parallel::makePSOCKcluster}.
#'        Note that for faster execution of the package, OPENBLAS or MKL library is recommended.
#'        When OPENBLAS is used, add environment variable \code{Sys.setenv(OPENBLAS_NUM_THREADS=1)}.
#'        When MKL is used,  add environment variables \code{Sys.setenv(MKL_NUM_THREADS=1)} and \code{Sys.setenv(MKL_THREADING_LAYER="GNU")}.
#' @param fix_overdispersion Logical or numeric. In general, when the sample size is sufficiently large (for example 10 samples per degree of freedom),
#' fix_overdispersion should be FALSE so that the overdispersion parameter is re-estimated in the reduced model.
#' However, when sample size is smaller, overdispersion parameter is hard to estimate and it is might be advantageous to
#' always use overdispersion parameter estimated under the full model, which is similar to how DESeq2 LRT works.
#' Alternatively, only for experimental purposes, 
#' a list of overdispersion parameters equal to the number of samples,
#' parametrized so that the variance of the negative binomial distribution is \code{mean + mean^2/overdispersion},
#' can be provided so that pre-computed overdispersion parameters can be used.
#'
#' @return
#' Returns a list mostly of matrices. Note that the matrices with "shrunken" in their names are only available when \code{shrunken_lfc} is \code{TRUE}: 
#' \describe{
#'   \item{p}{A matrix of G x H p-values}
#'   \item{padj}{A matrix of G x H p-values adjusted using Benjamini & Hochberg (1995).}
#'   \item{shrunken_lfc}{A matrix of shrunken log fold change between cell type-specific effects.}
#'   \item{shrunken_lfcSE}{A matrix of standard errors of shrunken log fold change of cell type-specific effects between different groups.}
#'   \item{shrunken_coefficients}{A matrix of shrunken coefficient estimates.}
#'   \item{shrunken_coefficientsSE}{A matrix of standard errors of shrunken coefficient estimates.}
#'   \item{lfc}{A matrix of log fold change between cell type-specific effects.}
#'   \item{lfcSE}{A matrix of standard errors of log fold change of cell type-specific effects between between different groups of cell type-specific effects.}
#'   \item{coefficients}{A matrix of MLE coefficient estimates.}
#'   \item{coefficientsSE}{A matrix of standard errors of coefficient estimates.}
#'   \item{overdispersion}{A matrix of overdispersion parameters. The overdispersion parameter
#'                         is parametrized so that the variance of the negative binomial distribution is \code{mean + mean^2/overdispersion}.}
#'   \item{elapsed_time}{The time elapsed.}
#' }
#' 
#' @examples
#' data("n50DE221rep1")
#' # As an example, run only the first 10 genes among the all 10,000 genes to save time:
#' res = run_CARseq(count_matrix = n50DE221rep1$observed_read_count[1:10,],
#'                  cellular_proportions = n50DE221rep1$rho,
#'                  groups = gl(2, 25),
#'                  formula = ~ RIN,
#'                  data = n50DE221rep1$clinical_variables,
#'                  read_depth = n50DE221rep1$d,
#'                  shrunken_lfc = TRUE,
#'                  cores = 1,
#'                  fix_overdispersion = FALSE
#' )
#'
#' @export
run_CARseq = function(
    count_matrix, cellular_proportions, groups, formula = NULL, data = NULL,
    read_depth = 1, 
    shrunken_lfc = TRUE,
    cores = 1, fix_overdispersion = FALSE) {
  
  elapsed = Sys.time()
  if (!is.matrix(count_matrix)) count_matrix = as.matrix(count_matrix)
  if (!is.matrix(cellular_proportions)) cellular_proportions = as.matrix(cellular_proportions)
  n = ncol(count_matrix)
  G = nrow(count_matrix)
  stopifnot(length(groups) == n)
  stopifnot(nrow(data) == n)
  stopifnot(nrow(cellular_proportions) == n)
  stopifnot(is.logical(fix_overdispersion) || (is.numeric(fix_overdispersion) && length(fix_overdispersion) == G))
  # "formula" contains covariates that are cell type-independent.
  # If an intercept term is detected, it will be automatically removed with a warning.
  # construct the matrix of "other_variables" from "formula" and "data":
  other_variables = NULL
  K = 0
  if (!is.null(formula)) {
    stopifnot(identical(class(formula), "formula"))
    other_variables = model.matrix(formula, data = as.data.frame(data))
    # remove the intercept term, which is covered in "groups"
    other_variables = other_variables[, apply(other_variables, 2, function(x) any(x != 1)), drop = FALSE]
    K = ncol(other_variables)
  }
  H = ncol(cellular_proportions)
  
  # construct the array of "cell_type_specific_variables"
  # (an array of n x H x M of cell type-independent variables) from "groups"
  # specify design matrix in full model
  groups = as.factor(groups)
  M = length(levels(groups))
  stopifnot(M >= 2)
  cell_type_specific_variables_full = array(0, dim=c(n, H, M))
  for (m in seq_len(M)) {
    cell_type_specific_variables_full[groups == levels(groups)[m], , m] = 1
  }
  if (is.null(colnames(count_matrix))) colnames(count_matrix) = paste0("sample", seq_len(n))
  if (is.null(colnames(cellular_proportions))) colnames(cellular_proportions) = paste0("celltype", seq_len(H))
  dimnames(cell_type_specific_variables_full) = list(colnames(count_matrix), colnames(cellular_proportions), levels(groups))
  
  # allocate cores for parallel computation
  `%dopar%` = foreach::`%do%`
  blocksize = 1
  job_schedule_list = list(seq_len(G))
  cl = NULL
  if (cores >= 2) {
    cl = parallel::makeCluster(cores, outfile="")  # FORK is not suitable for Windows or GUIs
    doParallel::registerDoParallel(cl)
    # doMC::registerDoMC(cores)
    `%dopar%` = foreach::`%dopar%`
    # foreach has too much computational overhead.
    # We would like to keep each computational blocksize sufficiently long.
    #  We set blocksize to 1 if number of genes < 100,
    #     set blocksize to 10 if 100 < number of genes < 10,000,
    # and set blocksize to 100 if number of genes >= 10,000.
    if (G >= 100 && G < 10000) {
      blocksize = 10
    } else if (G >= 10000) {
      blocksize = 100
    }
    job_schedule_list = split(seq_len(G), (seq_len(G) - 1) %/% blocksize)
  }
  
  # 1st pass: calculate LR statistics
  result_list = foreach::foreach(job=seq_along(job_schedule_list), .packages=c("MASS", "CARseq"), .combine=c) %dopar% {
    fit_model_list = list()
    for (j in job_schedule_list[[job]]) {
      fit_model_list[[j + 1 - job_schedule_list[[job]][1]]] = tryCatch({
        pvalues = rep(NA, H)
        overdispersion = FALSE
        if (is.numeric(fix_overdispersion)) overdispersion = fix_overdispersion[j]
        # call CARseq::fit_model to obtain estimates and negative log-likelihood
        res_optim_full = CARseq::fit_model(
          cell_type_specific_variables = cell_type_specific_variables_full,
          other_variables = other_variables,
          read_depth = read_depth,
          cellular_proportions = cellular_proportions,
          counts = count_matrix[j, ],
          init = NULL,
          fix_overdispersion = overdispersion,
          number_of_resample = 1,
          use_log_scale_algorithm = FALSE,
          lambda = 0,
          verbose = FALSE)
        
        if (is.logical(fix_overdispersion) && fix_overdispersion) overdispersion = res_optim_full$coefficients[nrow(res_optim_full$coefficients), 1]
        
        for (h in seq_len(H)) {
          # specify design matrix in reduced model
          cell_type_specific_variables_reduced = array(0, dim=c(n, H, M))
          for (m in seq_len(M)) {
            cell_type_specific_variables_reduced[groups == levels(groups)[m], , m] = 1
          }
          # For the cell type specified by indices:
          cell_type_specific_variables_reduced[, h, 1] = 1
          cell_type_specific_variables_reduced[, h, -1] = 0
          dimnames(cell_type_specific_variables_reduced) = list(colnames(count_matrix), colnames(cellular_proportions), levels(groups))
  
          # fit reduced model
          res_optim = CARseq::fit_model(
            cell_type_specific_variables = cell_type_specific_variables_reduced,
            other_variables = other_variables,
            read_depth = read_depth,
            cellular_proportions = cellular_proportions,
            counts = count_matrix[j, ],
            init = NULL,
            fix_overdispersion = overdispersion,
            number_of_resample = 1,
            use_log_scale_algorithm = FALSE,
            lambda = 0,
            verbose = FALSE)
  
          pvalues[h] = pchisq(2*(res_optim$value - res_optim_full$value), df = M - 1, lower.tail = FALSE)
        }
        list(pvalues        = pvalues,
             coefficients   = res_optim_full$coefficients,
             lfc            = res_optim_full$lfc
        )
      }, error = function(e) {
        NA
      }, finally = {
      })
    }
    if (j %% (10 * blocksize) == 0 || j == G) {
      cat(sprintf("Gene %d of %d has been processed at %s [1st pass]\n", j, G, Sys.time()))
    }
    fit_model_list
    
  }

  first_j_not_NA = {
    j = 1
    while (identical(NA, result_list[[j]])) j = j + 1
    j
  }
  coefficient_names = rownames(result_list[[first_j_not_NA]]$coefficients)
  coefficients_mat = do.call(rbind, lapply(result_list, function(x) {
    if (identical(NA, x)) {
      rep(NA, length(coefficient_names))
    } else {
      x$coefficients[,1]
    }
  }))
  coefficientsSE_mat = do.call(rbind, lapply(result_list, function(x) {
    if (identical(NA, x)) {
      rep(NA, length(coefficient_names))
    } else {
      x$coefficients[,2]
    }
  }))
  colnames(coefficients_mat) = colnames(coefficientsSE_mat) = coefficient_names
  rownames(coefficients_mat) = rownames(coefficientsSE_mat) = rownames(count_matrix)[seq_len(G)]
  pval_mat = do.call(rbind, lapply(result_list, function(x) if (identical(NA, x)) rep(NA, H) else x$pvalues))
  colnames(pval_mat) = colnames(cellular_proportions)
  rownames(pval_mat) = rownames(count_matrix)[seq_len(G)]
  # adjusted p values -- Benjamini-Hochberg is used here. Can also use qvalue::qvalue as an alternative:
  pval_adj_mat = pval_mat
  pval_adj_mat[] = apply(pval_adj_mat, 2, function(x) stats::p.adjust(x, method = "BH"))
  # MLE estimates and standard error estimators for log fold change of cell type-specific covariates:
  lfc_mat = do.call(rbind, lapply(result_list, function(x) {if (identical(NA, x)) rep(NA, H) else x$lfc[, "Estimate"]}))
  lfcSE_mat = do.call(rbind, lapply(result_list, function(x) {if (identical(NA, x)) rep(NA, H) else x$lfc[, "Std. Error"]}))
  rownames(lfc_mat) = rownames(lfcSE_mat) = rownames(count_matrix)[seq_len(G)]
  
  overdispersion_mat = coefficients_mat[, ncol(coefficients_mat), drop=FALSE]
  
  # 2nd pass run with prior when we need shrunken LFC
  if (shrunken_lfc) {
    
    lambda_empirical_bayes = estimate_shrinkage_prior(estimates=coefficients_mat[,seq_len(K+M*H),drop=FALSE], K=K, M=M, H=H)
    
    shrunken_result_list = foreach::foreach(job=seq_along(job_schedule_list), .packages=c("MASS", "CARseq"), .combine=c) %dopar% {
      fit_model_list = list()
      for (j in job_schedule_list[[job]]) {
        fit_model_list[[j + 1 - job_schedule_list[[job]][1]]] = tryCatch({
          pvalues = rep(NA, H)
          # call CARseq::fit_model to obtain estimates and negative log-likelihood
          res_optim_full = CARseq::fit_model(
            cell_type_specific_variables = cell_type_specific_variables_full,
            other_variables = other_variables,
            read_depth = read_depth,
            cellular_proportions = cellular_proportions,
            counts = count_matrix[j, ],
            init = NULL,
            fix_overdispersion = FALSE,
            number_of_resample = 1,
            use_log_scale_algorithm = TRUE,
            lambda = lambda_empirical_bayes,
            verbose = FALSE)
          
          list(coefficients   = res_optim_full$coefficients,
               lfc            = res_optim_full$lfc
          )
          
        }, error = function(e) {
          NA
        }, finally = {
        })
      }
      if (j %% (10 * blocksize) == 0 || j == G) {
        cat(sprintf("Gene %d of %d has been processed at %s [2nd pass for shrunken LFC]\n", j, G, Sys.time()))
      }
      fit_model_list
      
    }
    
    shrunken_coefficient_names = rownames(shrunken_result_list[[first_j_not_NA]]$coefficients)
    shrunken_coefficients_mat = do.call(rbind, lapply(shrunken_result_list, function(x) {
      if (identical(NA, x)) {
        rep(NA, length(shrunken_coefficient_names))
      } else {
        x$coefficients[,1]
      }
    }))
    shrunken_coefficientsSE_mat = do.call(rbind, lapply(shrunken_result_list, function(x) {
      if (identical(NA, x)) {
        rep(NA, length(shrunken_coefficient_names))
      } else {
        x$coefficients[,2]
      }
    }))
    colnames(shrunken_coefficients_mat) = colnames(shrunken_coefficientsSE_mat) = shrunken_coefficient_names
    rownames(shrunken_coefficients_mat) = rownames(shrunken_coefficientsSE_mat) = rownames(count_matrix)[seq_len(G)]
    shrunken_lfc_mat = do.call(rbind, lapply(shrunken_result_list, function(x) {if (identical(NA, x)) rep(NA, H) else x$lfc[, "Estimate"]}))
    shrunken_lfcSE_mat = do.call(rbind, lapply(shrunken_result_list, function(x) {if (identical(NA, x)) rep(NA, H) else x$lfc[, "Std. Error"]}))
    rownames(shrunken_lfc_mat) = rownames(shrunken_lfcSE_mat) = rownames(count_matrix)[seq_len(G)]
  }
  
  if (cores >= 2) parallel::stopCluster(cl)
  
  elapsed = Sys.time() - elapsed
  
  if (shrunken_lfc) {
    list(p = pval_mat,
         padj = pval_adj_mat,
         shrunken_lfc = shrunken_lfc_mat,
         shrunken_lfcSE = shrunken_lfcSE_mat,
         shrunken_coefficients = shrunken_coefficients_mat[, -ncol(shrunken_coefficients_mat), drop=FALSE],
         shrunken_coefficientsSE = shrunken_coefficientsSE_mat[, -ncol(shrunken_coefficients_mat), drop=FALSE],
         lfc = lfc_mat,
         lfcSE = lfcSE_mat,
         coefficients = coefficients_mat[, -ncol(coefficients_mat), drop=FALSE],
         coefficientsSE = coefficientsSE_mat[, -ncol(coefficients_mat), drop=FALSE],
         overdispersion = overdispersion_mat,
         lambda = lambda_empirical_bayes,
         elapsed_time = elapsed)
  } else {
    list(p = pval_mat,
         padj = pval_adj_mat,
         lfc = lfc_mat,
         lfcSE = lfcSE_mat,
         coefficients = coefficients_mat[, -ncol(coefficients_mat), drop=FALSE],
         coefficientsSE = coefficientsSE_mat[, -ncol(coefficients_mat), drop=FALSE],
         overdispersion = overdispersion_mat,
         elapsed_time = elapsed)
  }
}

#' fit a CARseq negative binomial model 
#' using iteratively weighted least squares / non-negative least squares
#' with backtracking line search.
#' 
#' \code{fit_model} fits a negative binomial distribution 
#' whose mean is a sum of nonnegative terms with covariates.
#' The overdispersion parameter is estimated by maximizing the 
#' adjusted profile log-likelihood.
#' 
#' Please use the wrapper \link{run_CARseq}. This function is not intended for general use.
#'
#' @param cell_type_specific_variables an array of n_B x H x K of cell type-independent variables.
#' @param other_variables a design matrix of n_B x M of cell type-specific variables. Can be NULL.
#' @param read_depth a vector of n_B sample-specific read depths. It it used as an offset term in the
#'        regression model. Alternatively, it can be 1, NA or NULL so that we don't have an offset term
#'        in the model, and log(read_depth) can be included as one of the cell type-independent variables.
#' @param cellular_proportions a matrix of n_B x H of cellular proportions.
#' @param counts A vector of n_B total read counts observed.
#' @param verbose logical. If \code{TRUE} (default), display information of negative log-likelihood of each iteration.
#' @param init ("expert" argument) a numeric vector of (K + H x M + 1) corresponding to the initial value of
#'             c(beta, gamma, overdispersion). Can be NULL. For gamma, log scale is assumed when
#'             \code{use_log_scale_algorithm} is TRUE, and non-log scale is assumed when \code{use_log_scale_algorithm} is FALSE.
#' @param fix_overdispersion ("expert" argument) logical (or numerical for a fixed overdispersion). 
#'        If \code{FALSE} (default), the overdispersion parameter will be estimated within the function.
#'        Otherwise, the overdispersion parameter can be 
#' @param number_of_resample ("expert" argument) numeric. The number of initial values we use in optimization. 
#'        The optimization algorithm generally is not sensitive to initial values, so we advise to leave it
#'        as the default value 1.
#' @param use_log_scale_algorithm ("expert" argument) logical. The default is FALSE, where the cell type-specific effects are
#'        fitted on a non-log scale and non-negative least squares is used (implemented in bvls).
#'        This is generally advisable than fitting them on a log scale, especially when the cell type-specific
#'        variables are binary (for example group indicators).
#' @param lambda ("expert" argument) numerical. \eqn{1/\sigma^2} in Gaussian prior. Only takes effect when 
#'        \code{use_log_scale_algorithm} is \code{TRUE}. Defaults to 0. It can also be a vector of length 
#'        \eqn{K + H * (M + 1)}.
#' @param return_coefficients_only ("expert" argument) logical. If \code{FALSE} (default), lfc, log-likelihood,
#'        Cook's distance and fitted values 
#'        will be provided together with coefficient estimates. Otherwise, only coefficient estimates are returned.
#' @param allow_negative_expression ("expert" argument) logical. If \code{FALSE} (default), the cell type-specific expression
#'        cannot be negative and will be returned in a log scale, along with LFCs. If \code{TRUE}, the cell type-specific expression
#'        will be on a non-log scale and could possibly be negative, while LFCs will not be returned.
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
                     verbose = TRUE,
                     init = NULL,
                     fix_overdispersion = FALSE,
                     number_of_resample = 1, 
                     use_log_scale_algorithm = FALSE,
                     lambda = 1e-6,
                     return_coefficients_only = FALSE,
                     allow_negative_expression = FALSE) {
  
  # cannot calculate LFC when allow_negative_expression = TRUE since log of negative expression is undefined:
  if (allow_negative_expression) {
    return_coefficients_only = TRUE
    use_log_scale_algorithm = FALSE
  }
  
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
  lambda_minimum = 1e-6
  if (!use_log_scale_algorithm) {
    lambda = 0
    lambda_minimum = 0
  }
  # Check whether design matrix is normal or expanded by checking the singularity of the cell type-specific design matrix.
  # A whole column being all 0s is removed before checking the singularity; such columns indicate reduced model.
  # NOTE: an alternative way is to let users specify this argument themselves 
  is_design_matrix_expanded = any(apply(cell_type_specific_variables, 2, function(x) kappa(x[, colSums(x^2) != 0, drop=FALSE])) > 1e9)
  if (is_design_matrix_expanded && all(lambda <= lambda_minimum)) {
    stop("The cell type-specific matrix is singular. Please specify a stronger prior and let use_log_scale_algorithm = TRUE.")
  }

  # We calculate which (cell type, cell type-specific variable) pairs are not included in the model
  # For future use.
  # is_active is a logical vector of length (K + H x M), recording which predictors are included in the model.
  is_reduced = apply(cell_type_specific_variables, c(2,3), function(y) all(y == 0))
  is_active = c(rep(TRUE, K), !is_reduced)
  if (!all(is_active) && use_log_scale_algorithm && sum(lambda) > lambda_minimum) {
    stop("Cannot add prior to a parameters in a reduced model. For Wald test, please fit the full model only. For LR test, please set prior to 0 and set use_log_scale_algorithm = FALSE.")
  }
  
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
  combined_upper = c(rep(exp(log_upper), K), rep(gamma_upper, sum(is_active) - K))
  if (allow_negative_expression) {
    combined_lower = -combined_upper
  } else {
    combined_lower = c(rep(-exp(log_upper), K), rep(gamma_lower, sum(is_active) - K))
  }
  
  # We need to sample many times to decide on an initial value.
  if (is.null(number_of_resample) || is.na(number_of_resample)) number_of_resample = 20
  # We only need number_of_resample different initial values and subsequent model fits.
  # However, we can try at most number_of_resample_max times.
  # We finally decide to keep them the same as the optimization algorithm is already good enough.
  number_of_resample_max = number_of_resample
  resample_size_list = rep(3 * (K + H * M), number_of_resample_max)    # An arbitrary number that should be larger than the number of parameters in the model
  resample_size_list[1] = n_B    # in the first instance, the initial value is based on all samples instead of a small subset of samples
  
  # Ridge regression penalty
  # lambda_vector = rep(lambda, K + H * M)
  if (length(lambda) == 1) {
    lambda_vector = c(rep(lambda_minimum, K), rep(lambda, H * M))
    for (h in seq_len(H)) {
      for (m in seq_len(M)) {
        if (all(1 == cell_type_specific_variables[ , h, m])) {
          lambda_vector[K + H * (m-1) + h] = lambda_minimum   # a less informative prior for cell type-specific expression
        }
      }
    }
    lambda_vector_expanded = c(lambda_vector, rep(lambda_minimum, H))
  } else {
    lambda_vector_expanded = lambda
    lambda_vector = lambda[seq_len(K + H * M)]
  }
  
  if (is.numeric(fix_overdispersion)) {
    overdispersion_theta = fix_overdispersion
    fix_overdispersion = TRUE
  } else if (fix_overdispersion) {
    overdispersion_theta = coefs[length(coefs)]
  }
  
  iter_overdispersion_update = 0
  overdispersion_has_converged = FALSE
  # we will alternate between updating mean model parameters and overdispersion parameter:
  while (!overdispersion_has_converged) {

    overdispersion_has_converged = TRUE  # it will be flipped to FALSE if later we detect no convergence
    negloglik_all_time_best = NULL
    coefs_all_time_best = NULL
    number_of_successful_iter = 0
    
    # Need to generate multiple initial values. 
    # If we have initial value, we need to use it first.
    negloglik_curr_vector = rep(NA, number_of_resample_max)  # stores a vector of likelihood for debugging purposes
    for (iter in seq_len(number_of_resample_max)) {
      resample_size = resample_size_list[iter]
      
      # save the current likelihood for each optimization path with a new initial value
      negloglik_curr_vector[iter] = suppressWarnings(tryCatch({
        
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
            # Sometimes glm.nb cannot converge. Since the estimates here are only used as initial values, this is not a big issue and is ignored:
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
                            family = MASS::negative.binomial(overdispersion_theta))
            } else {
              nbmodel = stats::glm(counts~offset(log(read_depth)) + 0 + other_variables + variables_sampled_cell_type_specific, 
                                     subset = sample(c(rep(TRUE, resample_size), rep(FALSE, n_B - resample_size))),
                            family = MASS::negative.binomial(overdispersion_theta))
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
        } # else if (!use_log_scale_algorithm) {
        #   # We need to take exp to make gamma in non-log scale within "coefs"
        #   coefs[seq_len(H * M) + K] = exp(coefs[seq_len(H * M) + K])
        # }
        
        # TODO: use cooks.distance() and fitted() to process the outliers
        # https://stats.stackexchange.com/questions/11632/what-kind-of-residuals-and-cooks-distance-are-used-for-glm
        # At the current development stage, we did not calculate cooks distance from our fitted model.
        # Instead, we fitted a glm model to replace the outliers outside the function.
        
        # Convergence conditions of Iteratively Weighted Least Squares
        epsilon_convergence = 1e-6
        maxiter = 100
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
          if (use_log_scale_algorithm) {
            # read depth, other effects, cell type-specific effects
            mu_matrix[] = 0
            for (m in seq_len(M)) {
              mu_matrix = mu_matrix + (matrix(rep(gamma[, m], each = n_B), nrow = n_B) * cell_type_specific_variables[, , m])
            }
            mu_matrix = exp(mu_matrix)
          } else {
            mu_matrix[] = 1
            for (m in seq_len(M)) {
              mu_matrix = mu_matrix * (matrix(rep(gamma[, m], each = n_B), nrow = n_B) ^ cell_type_specific_variables[, , m])
            }
          }
          mu_matrix = mu_matrix * cellular_proportions
          mu_matrix = mu_matrix * as.numeric(covariate_adjusted_read_depth)
          mu = rowSums(mu_matrix)
          
          # current negative log-likelihood:
          negloglik_prev = negloglik_curr
          negloglik_curr = - sum(lgamma(counts + overdispersion_theta) - lgamma(overdispersion_theta) - lgamma(counts + 1) +
                                   overdispersion_theta*log(overdispersion_theta) + counts*log(mu) -
                                   (counts + overdispersion_theta) * log(overdispersion_theta+mu)) +
                                sum(lambda_vector * coefs[-length(coefs)]^2 * 0.5)
                                
          if (!is.finite(negloglik_curr)) {
            if (verbose) message("The optimization failed to converge.")
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
            # direction = - coefs[-length(coefs)][is_active] +
            #   tryCatch(stats::lsfit(adjusted_design_matrix[, is_active],
            #                         adjusted_response,
            #                         wt = weights,
            #                         intercept = FALSE)$coef, 
            #                      warning = function(w) {coefs[-length(coefs)][is_active]},
            #                      error = function(e) {coefs[-length(coefs)][is_active]})
            # direction = - coefs[-length(coefs)][is_active] +
            #   tryCatch(as.numeric(glmnet::glmnet(adjusted_design_matrix[, is_active],
            #                           adjusted_response,
            #                           weights = weights,
            #                           alpha = 0,
            #                           lambda = lambda,
            #                           intercept = FALSE)$beta),
            #            warning = function(w) {coefs[-length(coefs)][is_active]},
            #            error = function(e) {coefs[-length(coefs)][is_active]})
            direction = - coefs[-length(coefs)][is_active] +
              tryCatch({
                tXWX_lambda = crossprod(adjusted_design_matrix[, is_active], weights * adjusted_design_matrix[, is_active])
                diag(tXWX_lambda) = diag(tXWX_lambda) + lambda_vector[is_active]
                L = chol(tXWX_lambda)
                w = backsolve(L, t(adjusted_design_matrix[, is_active]) %*% (weights * adjusted_response), upper.tri = TRUE, transpose = TRUE)
                forwardsolve(L, w, upper.tri = TRUE, transpose = FALSE)
              },
              warning = function(w) {coefs[-length(coefs)][is_active]},
              error = function(e) {coefs[-length(coefs)][is_active]})
            
          } else {
            # bounded-variable least squares ensures that the active parameters in H*M cell type-specific ones are bounded
            # direction = - coefs[-length(coefs)][is_active] +
            #   bvls::bvls(adjusted_design_matrix[, is_active] * sqrt(weights),
            #              adjusted_response * sqrt(weights),
            #              bl = combined_lower,
            #              bu = combined_upper)$x
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
          
          # Define a function that can decide whether a coefficient estimate is out of the boundary
          # If "allow_negative_expression" is TRUE,
          # While we no longer need to check whether cell type-specific coefficient estimates are all non-negative,
          # we need to check whether the predicted expression of each sample is positive or not.
          is_out_of_boundary = function(coefs, is_active, lower, upper,
                                        allow_negative_expression = FALSE, cell_type_specific_variables = NULL) {
            
            is_out = any(upper < coefs[-length(coefs)][is_active])
            
            if (allow_negative_expression) {
              # If we allow negative cell type-specific estimates, we would
              # need to check whether the predicted expression of each sample is positive or not.
              stopifnot(!is.null(cell_type_specific_variables))
              n_B = dim(cell_type_specific_variables)[1]
              H   = dim(cell_type_specific_variables)[2]
              M   = dim(cell_type_specific_variables)[3]
              gamma = matrix(coefs[seq_len(H * M) + K], nrow = H, ncol = M)
              mu_matrix[] = matrix(1, n_B, H)
              for (m in seq_len(M)) {
                mu_matrix = mu_matrix * (matrix(rep(gamma[, m], each = n_B), nrow = n_B) ^ cell_type_specific_variables[, , m])
              }
              mu_matrix = mu_matrix * cellular_proportions
              mu = rowSums(mu_matrix)
              is_positive = all(mu > 0)
              is_out = is_out || !is_positive
            } else {
              is_out = is_out || any(lower > coefs[-length(coefs)][is_active])
            }
            
            is_out
          }
          
          for (backtracking_iter in 1:max_backtracking_iter) {
            coefs_new = coefs
            coefs_new[-length(coefs_new)][is_active] = coefs[-length(coefs)][is_active] + step_length * direction
            
            # if the parameters are already out of the boundary, we half the step size until we are back inside the boundary:
            step_iter = 1
            max_step_iter = 20
            while (backtracking_iter == 1 &&
                   is_out_of_boundary(coefs_new, is_active, combined_lower, combined_upper, allow_negative_expression, cell_type_specific_variables) &&
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
                                    verbose = FALSE) +
                            sum(lambda_vector * coefs_new[-length(coefs_new)]^2 * 0.5)
  
            if (verbose) message(sprintf("Negative log-likelihood = %.4f at iteration %d.%d",
                                         negloglik_new, inner_iter, backtracking_iter))
            if (!is.finite(negloglik_new)) {
              negloglik_new = Inf
              break
            }
            if (negloglik_new <= negloglik_curr + c1 * step_length * crossprod(-score_curr, direction))
              break
            step_length = step_length * contraction
          }
          if (verbose) cat("\n")
          if (backtracking_iter != max_backtracking_iter) {
            coefs = coefs_new
            negloglik_curr = negloglik_new
          }
          
          # save the best coefs and negative log-likelihood
          if (is.null(negloglik_best) || negloglik_curr < negloglik_best) {
            coefs_best = coefs
            negloglik_best = negloglik_curr
          }
          if (is.null(negloglik_all_time_best) || negloglik_curr < negloglik_all_time_best) {
            coefs_all_time_best = coefs
            negloglik_all_time_best = negloglik_curr
          }
          
          if (backtracking_iter == max_backtracking_iter) break
          
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
        if (inner_iter == maxiter && verbose) {
          message(paste("Max number of iterations", maxiter, "has been reached."))
        }
        
        number_of_successful_iter = number_of_successful_iter + 1
        
        if (number_of_successful_iter > number_of_resample) break
        
        negloglik_best
      # }, error = function(e) {warning(e); warning("Iteration has failed."); NA}, warning = function(w) {warning(w); NA})
      }, error = function(e) {warning(e); warning("Iteration has failed."); NA}))
      
      # We will generate initial values from glm next iteration:
      init = NULL
      if (verbose) cat("\n")
      
    }
    
    # throw error if we do not have number_of_resample successful tries:
    # NOTE: should we add it back?
    # if (number_of_successful_iter != number_of_resample) stop("Optimization using different initial values has failed too many times!")
    
    # The best solution after using multiple initial values:
    coefs = coefs_all_time_best
    if (use_log_scale_algorithm) {
      coefs[seq_len(H * M) + K] = exp(coefs[seq_len(H * M) + K])
    }
      
    # update overdispersion parameter
    # NOTE: when adding normal prior to beta, need to use standard (not expanded) design matrix.
    if (!fix_overdispersion) {
      
      # calculate likelihood without updated overdispersion to retrieve mu_matrix
      objective_old_overdispersion = negloglik(coefs = coefs,
                                               cell_type_specific_variables = cell_type_specific_variables,
                                               other_variables = other_variables,
                                               read_depth = read_depth,
                                               cellular_proportions = cellular_proportions,
                                               counts = counts,
                                               use_log_scale_params = FALSE,
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
        use_log_scale_params = FALSE,
        method = "Brent",
        lower = -5,
        upper = 10
        # lower = log(coefs[length(coefs)]) - 1,
        # upper = log(coefs[length(coefs)]) + 1
      )
      
      # When overdispersion parameter is update, we only re-fit mean regression parameters when it is a major update:
      if (abs(optimize_theta$par - log(coefs[length(coefs)])) > 0.1) {
        # we will continue to update mean regression parameters:
        overdispersion_has_converged = FALSE
        iter_overdispersion_update = 1 + iter_overdispersion_update
        number_of_resample_max = number_of_resample = 1
        if (use_log_scale_algorithm) {
          coefs[seq_len(H * M) + K] = log(coefs[seq_len(H * M) + K])
        }
        init = coefs
      }
      if (iter_overdispersion_update >= 10) {
        stop("Overdispersion parameter failed to converge!")
      }
      
      coefs[length(coefs)] = exp(optimize_theta$par)
    }
  }
  
  # finally, a likelihood without using adjusted profile likelihood
  objective = negloglik(coefs = coefs,
                        cell_type_specific_variables = cell_type_specific_variables,
                        other_variables = other_variables,
                        read_depth = read_depth,
                        cellular_proportions = cellular_proportions,
                        counts = counts,
                        use_log_scale_params = FALSE,
                        verbose = TRUE)
  
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
  
  
  # A "lite" version of to return results without reporting cooks distance, LFC change, etc.
  if (return_coefficients_only || allow_negative_expression) {
    if (use_log_scale_algorithm) {
      coefs[seq_len(H * M) + K] = log(coefs[seq_len(H * M) + K])
    }
    coef_matrix = matrix(coefs, nrow=length(coefs), ncol=1)
    colnames(coef_matrix) = "Estimate"
    rownames(coef_matrix) =
      c(colnames_other_variables, 
        paste(rep(dimnames(cell_type_specific_variables)[[3]], each=H),
              dimnames(cell_type_specific_variables)[[2]],
              sep=":"),
        "overdispersion")
    
    if (!allow_negative_expression) {
      return(list(coefficients = coef_matrix,
                  value = as.numeric(objective)))
    } else {
      # TODO: The "diff_matrix" Wald test does not actually work properly due to the fact that quadratic approximation
      # of likelihood surface at the MLE cannot work well. It is currently still included for debugging purposes.
      # There is no justification to believe the statistics and p-values reported here.
      
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
      # the overdispersion used in obtaining standard errors (Wald test) are still the MLE estimate
      # overdispersion_theta = coefs[length(coefs)]  # use the overdispersion parameter estimated using profile likelihood
      weights = overdispersion_theta / (mu * (mu+overdispersion_theta))
      adjusted_design_matrix = cbind(other_variables * mu,
                                     matrix(cell_type_specific_variables, nrow=n_B) * as.numeric(mu_matrix))
      tXWX = crossprod(adjusted_design_matrix[, is_active], weights * adjusted_design_matrix[, is_active])
      L = chol(tXWX)
      w = backsolve(L, tXWX, upper.tri = TRUE, transpose = TRUE)
      fit = forwardsolve(L, w, upper.tri = TRUE, transpose = FALSE)
      fit = t(fit)
      w = backsolve(L, fit, upper.tri = TRUE, transpose = TRUE)
      solve_tXWX = forwardsolve(L, w, upper.tri = TRUE, transpose = FALSE)
      # build diff matrix & z values
      diff_matrix = NULL
      # These coefficients are actually not in the (reduced) model are set to NAs:
      # calculate logFoldChange and diffSE
      if (M >= 2) {
        diff_matrix = matrix(NA, nrow=(M-1)*H, ncol=4)
        colnames(diff_matrix) = c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
        # assume we use group 1 of 1..M as the baseline:
        is_active_matrix_form = matrix(is_active[seq_len(H*M)+K], nrow=H)
        add_NA_status = function(x, NA_status) {x[NA_status] = NA; x}
        for (m in 2:M) {
          NA_status = !(is_active_matrix_form[, m] & is_active_matrix_form[, 1])
          diff_matrix_current_rows = (H * (m-2) + 1):(H * (m-1))
          diff_matrix[diff_matrix_current_rows, 1] = 
            add_NA_status(coef_matrix[(K + H * (m-1) + 1):(K + H * m), 1] - coef_matrix[(K + 1):(K + H), 1], NA_status)
          # contrast matrix
          R = matrix(0, nrow = H, ncol = K + H * M)
          for (h in seq_len(H)) {
            R[h, c(K + H * (m-1) + h, K + h)] = c(1, -1)
          }
          R = R[, is_active , drop=FALSE]
          # equivalent to sqrt(diag(R %*% solve(tXWX) %*% t(R)))
          diff_matrix[diff_matrix_current_rows, 2] = add_NA_status(
            sqrt(diag(R %*% solve_tXWX %*% t(R))), NA_status)
        }
        # z value
        diff_matrix[, 3] = diff_matrix[, 1] / diff_matrix[, 2]
        # Pr(>|z|)
        diff_matrix[, 4] = pmin(1, 2 * pnorm(-abs(diff_matrix[, 3])))
      }
      coef_matrix[K + which(is_reduced), ] = NA

      return(list(coefficients = coef_matrix,
                  diff_matrix = diff_matrix,
                  value = as.numeric(objective)))
    }
  }
  
  if (!allow_negative_expression) {
    
    coefs[seq_len(H * M) + K] = log(coefs[seq_len(H * M) + K])
    
    # fit model using expanded matrices
    # NOTE: "fit_model" is a recursive call (executed only once since return_coefficients_only==TRUE). We might restructure it in the future:
    if (use_log_scale_algorithm && sum(lambda) > lambda_minimum && !is_design_matrix_expanded) {
      
      # add cell type-specific intercepts to make the "lambda" prior symmetric on all levels
      is_design_matrix_expanded = TRUE
      add_intercept_to_dimnames = function(x) {x[[3]] = c(x[[3]], "intercept"); x}
      cell_type_specific_variables_expanded = array(data = cell_type_specific_variables,
                                                    dim = dim(cell_type_specific_variables) + c(0,0,1),
                                                    dimnames = add_intercept_to_dimnames(dimnames(cell_type_specific_variables)))
      cell_type_specific_variables_expanded[, , dim(cell_type_specific_variables_expanded)[3]] = 1
      
      coefs_init = c(coefs[seq_len(K)],                                   # cell type-independent variables
                rep(0, (M + 1) * H),                                      # cell type-specific variables
                overdispersion = overdispersion_theta)
      coefs = fit_model(cell_type_specific_variables = cell_type_specific_variables_expanded,
                        other_variables = other_variables,
                        read_depth = read_depth,
                        cellular_proportions = cellular_proportions,
                        counts = counts,
                        init = coefs_init,
                        fix_overdispersion = TRUE,
                        number_of_resample = number_of_resample, 
                        use_log_scale_algorithm = TRUE,
                        lambda = lambda,
                        verbose = verbose, 
                        return_coefficients_only = TRUE,
                        allow_negative_expression = FALSE)$coefficients[, "Estimate"]
      
      # use expanded design matrix
      cell_type_specific_variables_standard = cell_type_specific_variables
      cell_type_specific_variables = cell_type_specific_variables_expanded
      M = M + 1
      lambda_vector = lambda_vector_expanded
    }
    
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
    # the overdispersion used in obtaining standard errors (Wald test) are still the MLE estimate
    # overdispersion_theta = coefs[length(coefs)]  # use the overdispersion parameter estimated using profile likelihood
    weights = overdispersion_theta / (mu * (mu+overdispersion_theta))
    adjusted_design_matrix = cbind(other_variables * mu,
                                   matrix(cell_type_specific_variables, nrow=n_B) * as.numeric(mu_matrix))
    # Sometimes a column is all 0;
    # sometimes parameters are at the boundary.
    # Need to remove them from active covariates to ensure tXWX is invertible.
    if (use_log_scale_algorithm && sum(lambda) > lambda_minimum) {
      # with sufficient shrinkage, all of the parameters should be "active", i.e. not on the boundary
      is_active = c(rep(TRUE, K), !is_reduced, rep(TRUE, H)) #&   # M is M + 1 here
        #(coefs[seq_len(K + H * M)] - log_lower >= .Machine$double.eps^0.5) &
        #(coefs[seq_len(K + H * M)] - log_upper <= -.Machine$double.eps^0.5)
    } else {
      is_active = is_active &
        !apply(adjusted_design_matrix, 2, function(x) isTRUE(all.equal(sum(abs(x)), 0))) &
        (coefs[seq_len(K + H * M)] - log_lower >= .Machine$double.eps^0.5) &
        (coefs[seq_len(K + H * M)] - log_upper <= -.Machine$double.eps^0.5)
    }
      
    tXWX = crossprod(adjusted_design_matrix[, is_active], weights * adjusted_design_matrix[, is_active])
    tXWX_lambda = tXWX + diag(lambda_vector[is_active])
    # solve_tXWX_lambda = solve(tXWX_lambda)
    # solve_tXWX_lambda = solve_tXWX_lambda %*% tXWX %*% solve_tXWX_lambda
    L = chol(tXWX_lambda)
    w = backsolve(L, tXWX, upper.tri = TRUE, transpose = TRUE)
    ridge_fit = forwardsolve(L, w, upper.tri = TRUE, transpose = FALSE)
    ridge_fit = t(ridge_fit)
    w = backsolve(L, ridge_fit, upper.tri = TRUE, transpose = TRUE)
    solve_tXWX_lambda = forwardsolve(L, w, upper.tri = TRUE, transpose = FALSE)
    
    # SEE = sqrt(diag(solve(tXWX)))
    
    # build the matrix of coefficients
    coef_matrix = matrix(NA, nrow=length(coefs), ncol=2)
    colnames(coef_matrix) = c("Estimate", "Std. Error")
    rownames(coef_matrix) = names(coefs)
    coef_matrix[, 1] = coefs
    coef_matrix[, 2][-length(coefs)][is_active] = sqrt(diag(solve_tXWX_lambda))
    
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
      is_active_matrix_form = matrix(is_active[seq_len(H*M)+K], nrow=H)
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
            sqrt(diag(R %*% solve_tXWX_lambda %*% t(R))), NA_status)
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
    hat_values = diag(sqrtWX %*% solve_tXWX_lambda %*% t(sqrtWX))
    # variance of negative binomial is 1/weights
    # standardised Pearson residuals
    r_Pi = (counts - mu) / sqrt((1 - hat_values) / weights)
    cooks = 1/p * hat_values / (1 - hat_values) * r_Pi^2
    
    # use standard design matrix
    if (use_log_scale_algorithm && is_design_matrix_expanded) {
      cell_type_specific_variables = cell_type_specific_variables_standard
      M = M - 1
    }
    
    # Name cell type-specific effects as interaction terms
    if (use_log_scale_algorithm && is_design_matrix_expanded) {
      rownames(coef_matrix)[-length(coefs)] = 
          c(colnames_other_variables, 
            paste(rep(c(dimnames(cell_type_specific_variables)[[3]], "intercept"), each=H),
                  dimnames(cell_type_specific_variables)[[2]],
                  sep=":"))
      if (M >= 2) {
        lfc_matrix = lfc_matrix[seq_len(H * (M - 1)), ]
      }
    } else {
      rownames(coef_matrix)[-length(coefs)] = 
        c(colnames_other_variables, 
          paste(rep(dimnames(cell_type_specific_variables)[[3]], each=H),
                dimnames(cell_type_specific_variables)[[2]],
                sep=":"))
    }
    
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
  # if (use_log_scale_params) {
  #   for (i in seq_len(n_B)) {
  #     # read depth, other effects, cell type-specific effects
  #     mu_matrix[i, ] = cellular_proportions[i, ] * exp(rowSums(gamma * cell_type_specific_variables[i, , ]))
  #   }
  # } else {
  #   for (i in seq_len(n_B)) {
  #     mu_matrix[i, ] = cellular_proportions[i, ] * matrixStats::rowProds(gamma ^ cell_type_specific_variables[i, , ])
  #   }
  # }
  if (use_log_scale_params) {
    # read depth, other effects, cell type-specific effects
    mu_matrix[] = 0
    for (m in seq_len(M)) {
      mu_matrix = mu_matrix + (matrix(rep(gamma[, m], each = n_B), nrow = n_B) * cell_type_specific_variables[, , m])
    }
    mu_matrix = exp(mu_matrix)
  } else {
    mu_matrix[] = 1
    for (m in seq_len(M)) {
      mu_matrix = mu_matrix * (matrix(rep(gamma[, m], each = n_B), nrow = n_B) ^ cell_type_specific_variables[, , m])
    }
  }
  mu_matrix = mu_matrix * cellular_proportions
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
#' whose mean is a sum of nonnegative terms with covariates.
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
      
      epsilon_convergence = 1e-4
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
        warning(sprintf("Max number of iterations %d has been reached.", maxiter))
      }
      
      number_of_successful_iter = number_of_successful_iter + 1
      
      if (number_of_successful_iter == number_of_resample) break
      
      negloglik_curr
    
    }, error = function(e) {stop("Iteration has failed.")}, warning = function(w) {stop(w)})
    
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



